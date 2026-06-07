"""Structure parsing + PEAK — the only structure work the pipeline needs.

In this pipeline the fancy interface metrics (iLIS/LIS/cLIS/ipSAE/actifpTM) all come
from vendored lis.py. The one thing lis.py does NOT output is PEAK
( = 1 − min-interchain-PAE/30 ), which needs chain boundaries from the structure.
So this module = a dependency-light PDB/mmCIF parser (representative atom Cβ, Cα for
Gly) + a PAE loader + PEAK. (Standalone actifpTM lives in the repo-root actifptm.py,
used by the legacy collectors and as the r=0.92 cross-check against lis.py.)
"""
import glob
import numpy as np


def _rep_atoms_pdb(path):
    best, order = {}, []
    for ln in open(path):
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        atom = ln[12:16].strip(); resn = ln[17:20].strip()
        key = (ln[21], ln[22:26].strip(), ln[26])
        want = 2 if (atom == "CB" and resn != "GLY") else (1 if atom == "CA" else 0)
        if not want:
            continue
        try:
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        except ValueError:
            continue
        if key not in best:
            best[key] = (want, xyz); order.append(key)
        elif want > best[key][0]:
            best[key] = (want, xyz)
    return [k[0] for k in order], np.array([best[k][1] for k in order], float)


def _rep_atoms_cif(path):
    header, rows, in_loop = [], [], False
    for ln in open(path):
        s = ln.strip()
        if s.startswith("_atom_site."):
            header.append(s); in_loop = True; continue
        if in_loop and (s.startswith("_") or s.startswith("loop_") or s == "#"):
            if rows:
                break
            continue
        if in_loop and header and s and not s.startswith("_"):
            rows.append(s.split())
    idx = {n.split(".")[1]: i for i, n in enumerate(header)}
    def c(*names):
        for n in names:
            if n in idx:
                return idx[n]
    ci = dict(atom=c("label_atom_id", "auth_atom_id"), comp=c("label_comp_id", "auth_comp_id"),
              ch=c("auth_asym_id", "label_asym_id"), seq=c("auth_seq_id", "label_seq_id"),
              x=c("Cartn_x"), y=c("Cartn_y"), z=c("Cartn_z"), model=c("pdbx_PDB_model_num"))
    best, order = {}, []
    for r in rows:
        if ci["model"] is not None and r[ci["model"]] not in ("1", "."):
            continue
        atom = r[ci["atom"]].strip('"'); resn = r[ci["comp"]]
        key = (r[ci["ch"]], r[ci["seq"]])
        want = 2 if (atom == "CB" and resn != "GLY") else (1 if atom == "CA" else 0)
        if not want:
            continue
        xyz = (float(r[ci["x"]]), float(r[ci["y"]]), float(r[ci["z"]]))
        if key not in best:
            best[key] = (want, xyz); order.append(key)
        elif want > best[key][0]:
            best[key] = (want, xyz)
    return [k[0] for k in order], np.array([best[k][1] for k in order], float)


def parse_structure(path):
    """Return (asym_id [N int chain index in order of appearance], coords [N,3])."""
    chains, coords = (_rep_atoms_cif(path) if path.endswith((".cif", ".mmcif"))
                      else _rep_atoms_pdb(path))
    seen, asym = {}, []
    for ch in chains:
        if ch not in seen:
            seen[ch] = len(seen)
        asym.append(seen[ch])
    return np.array(asym), coords


def load_pae(path):
    d = np.load(path)
    return d[d.files[0]]


def find_structure(pdir, model):
    for pat in ("*_model_%d.cif" % model, "*_model_%d.pdb" % model):
        h = sorted(glob.glob("%s/%s" % (pdir, pat)))
        if h:
            return h[0]


def peak_per_chainpair(pae, asym_id, cutoff=30.0):
    """{(chr_i, chr_j): PEAK} for each interchain pair. PEAK = max(0, 1 - minPAE/cutoff)."""
    chains = sorted(set(asym_id.tolist()))
    bnd = {c: (np.where(asym_id == c)[0][0], np.where(asym_id == c)[0][-1]) for c in chains}
    out = {}
    for a in range(len(chains)):
        for b in range(a + 1, len(chains)):
            si, ei = bnd[chains[a]]; sj, ej = bnd[chains[b]]
            blk = np.concatenate([pae[si:ei + 1, sj:ej + 1].ravel(),
                                  pae[sj:ej + 1, si:ei + 1].ravel()])
            out[(chr(65 + a), chr(65 + b))] = round(max(0.0, 1 - float(blk.min()) / cutoff), 4)
    return out
