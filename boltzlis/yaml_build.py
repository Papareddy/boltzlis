"""Build Boltz-2 input YAMLs from resolved chain-sets.

Two modes:
  grid    : every A-id x every B-id -> one 2-chain YAML per pair (the screen case).
            "single or multiple" on each side just changes how many pairs you get.
  complex : ONE YAML containing all A-ids + all B-ids as separate chains
            (multi-chain assembly).
Chain IDs in the YAML are A, B, C, ... in order.
"""
import os, string

CHAIN_IDS = list(string.ascii_uppercase)


def _yaml(chains):
    """chains: list of (chain_id, sequence) -> Boltz v1 YAML text."""
    out = ["version: 1", "sequences:"]
    for cid, seq in chains:
        out += ["  - protein:", "      id: %s" % cid, "      sequence: %s" % seq]
    return "\n".join(out) + "\n"


def safe(name):
    return "".join(c if (c.isalnum() or c in "-._") else "_" for c in name)


def build_grid(a_set, b_set, outdir):
    """a_set,b_set: list of (label, seq). Writes one 2-chain YAML per (a,b) pair.
    Returns list of (name, yaml_path)."""
    os.makedirs(outdir, exist_ok=True)
    jobs = []
    for al, aseq in a_set:
        for bl, bseq in b_set:
            name = "%s__%s" % (safe(al), safe(bl))
            p = os.path.join(outdir, name + ".yaml")
            open(p, "w").write(_yaml([("A", aseq), ("B", bseq)]))
            jobs.append((name, p))
    return jobs


def build_complex(a_set, b_set, name, outdir):
    """Single multi-chain YAML: all A then all B as chains A,B,C,...
    Returns [(name, yaml_path)]."""
    os.makedirs(outdir, exist_ok=True)
    chains = []
    for i, (_, seq) in enumerate(list(a_set) + list(b_set)):
        chains.append((CHAIN_IDS[i], seq))
    p = os.path.join(outdir, safe(name) + ".yaml")
    open(p, "w").write(_yaml(chains))
    return [(safe(name), p)]
