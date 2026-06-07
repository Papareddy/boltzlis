"""collect_all -- one call to compute EVERY interface metric for a Boltz out dir.

Runs vendored AFM-LIS lis.py (iLIS/LIS/cLIS/LIA/ipSAE/actifpTM/ipTM, parallel -w),
adds PEAK ( = 1 - min_interchain_PAE/30 ) per chain pair, and aggregates the per-model
rows to a ranked per-pair table (mean + max over the N diffusion models).

CLI:  python -m boltzlis.collect <out_dir> -o results.tsv [-w 8] [--rank iLIS_max]
"""
import os, sys, glob, argparse, subprocess
import numpy as np, pandas as pd
from . import structure as st

HERE = os.path.dirname(os.path.abspath(__file__))
LIS_PY = os.path.join(HERE, "lis.py")
# metrics to aggregate (mean + max across models); others are kept from the best model
AGG = ["iLIS", "LIS", "cLIS", "iLIA", "LIA", "cLIA", "ipSAE", "actifpTM", "ipTM",
       "pTM", "PEAK", "pLDDT_i", "pLDDT_j"]


def peak_table(out_dir):
    """PEAK per (name, model, chain_i, chain_j) from PAE npz + structure.
    Uses os.walk (robust on gpfs / symlinks) rather than a recursive glob."""
    rows = []
    for root, _, files in os.walk(out_dir):
        for fn in files:
            if not (fn.startswith("pae_") and "_model_" in fn and fn.endswith(".npz")):
                continue
            pf = os.path.join(root, fn); name = os.path.basename(root)
            try:
                m = int(fn.rsplit("_model_", 1)[1].split(".")[0])
            except Exception:
                continue
            sf = st.find_structure(root, m)
            if sf is None:
                continue
            try:
                pae = st.load_pae(pf); asym, _ = st.parse_structure(sf)
            except Exception:
                continue
            if pae.shape[0] != len(asym):
                continue
            for (ci, cj), peak in st.peak_per_chainpair(pae, asym).items():
                rows.append(dict(name=name, model=m, chain_i=ci, chain_j=cj, PEAK=peak))
    return pd.DataFrame(rows)


def run_lis(out_dir, csv_path, python_exe=None, lis_py=LIS_PY, workers=8):
    python_exe = python_exe or sys.executable
    d = os.path.dirname(os.path.abspath(csv_path)) or "."
    os.makedirs(d, exist_ok=True)
    cmd = [python_exe, lis_py, out_dir, "-w", str(workers),
           "-o", os.path.basename(csv_path), "-d", d]
    subprocess.run(cmd, check=True)
    return pd.read_csv(csv_path)


def collect_all(out_dir, out_tsv, python_exe=None, lis_py=LIS_PY, workers=8,
                rank_by="iLIS_max"):
    lis_csv = os.path.splitext(out_tsv)[0] + ".permodel.csv"
    lis = run_lis(out_dir, lis_csv, python_exe, lis_py, workers)
    peak = peak_table(out_dir)
    pk = {(r["name"], int(r["model"]), r["chain_i"], r["chain_j"]): r["PEAK"]
          for _, r in peak.iterrows()}
    lis["PEAK"] = [pk.get((r["name"], int(r["model"]), r["chain_i"], r["chain_j"]), np.nan)
                   for _, r in lis.iterrows()]
    sys.stderr.write("[collect] peak rows=%d, PEAK matched %d/%d lis rows\n"
                     % (len(peak), int(lis["PEAK"].notna().sum()), len(lis)))
    # aggregate per pair (name + chain pair) across models
    keys = ["name", "chain_i", "chain_j"]
    have = [c for c in AGG if c in lis.columns]
    g = lis.groupby(keys)
    agg = g[have].agg(["mean", "max"])
    agg.columns = ["%s_%s" % (m, s) for m, s in agg.columns]
    agg["n_models"] = g.size()
    agg = agg.reset_index()
    if rank_by not in agg.columns:
        rank_by = "PEAK_max" if "PEAK_max" in agg.columns else agg.columns[3]
    agg = agg.sort_values(rank_by, ascending=False)
    agg.to_csv(out_tsv, sep="\t", index=False)
    lis.to_csv(lis_csv, index=False)
    return agg


def main():
    ap = argparse.ArgumentParser(description="Compute all interface metrics for a Boltz out dir")
    ap.add_argument("out_dir", help="Boltz output dir (contains <pair>/boltz_results_*/predictions/...)")
    ap.add_argument("-o", "--output", default="metrics.tsv")
    ap.add_argument("-w", "--workers", type=int, default=8)
    ap.add_argument("--python", default=None, help="python exe with numpy+scipy for lis.py")
    ap.add_argument("--lis", default=LIS_PY)
    ap.add_argument("--rank", default="iLIS_max")
    a = ap.parse_args()
    agg = collect_all(a.out_dir, a.output, a.python, a.lis, a.workers, a.rank)
    print("wrote %s  (%d pairs)" % (a.output, len(agg)))
    cols = [c for c in ("name", "n_models", "iLIS_max", "PEAK_max", "actifpTM_max", "ipSAE_max") if c in agg.columns]
    print(agg[cols].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
