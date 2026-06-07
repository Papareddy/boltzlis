"""Local unit tests -- no cluster/GPU/network needed.
Run: python -m pytest tests/  (or python tests/test_metrics.py)"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from boltzlis import structure, yaml_build, fetch


def _toy_complex(n_a=10, n_b=8):
    """Two chains; place chain B residues near the first few of A to form an interface."""
    rng = np.random.RandomState(0)
    a = rng.rand(n_a, 3) * 50
    b = a[:n_b] + rng.rand(n_b, 3) * 2.0      # B sits ~2A from first n_b of A -> contacts
    coords = np.vstack([a, b])
    asym = np.array([0] * n_a + [1] * n_b)
    return asym, coords


def test_peak_per_chainpair():
    asym, coords = _toy_complex()
    N = len(asym)
    pae = np.full((N, N), 30.0); pae[0, N-1] = pae[N-1, 0] = 3.0
    pk = structure.peak_per_chainpair(pae, asym)
    assert abs(pk[('A','B')] - 0.9) < 1e-6, "min inter PAE 3 -> PEAK 0.9"


def test_peak_definition():
    asym, coords = _toy_complex()
    N = len(asym)
    pae = np.full((N, N), 30.0)
    pae[0, N - 1] = 3.0; pae[N - 1, 0] = 3.0      # one good inter-chain pair
    bnd = {c: (np.where(asym == c)[0][[0, -1]]) for c in (0, 1)}
    (si, ei), (sj, ej) = bnd[0], bnd[1]
    blk = np.concatenate([pae[si:ei + 1, sj:ej + 1].ravel(), pae[sj:ej + 1, si:ei + 1].ravel()])
    peak = max(0.0, 1 - float(blk.min()) / 30)
    assert abs(peak - 0.9) < 1e-9, "min inter PAE 3 -> PEAK 0.9"


def test_yaml_grid_and_complex(tmp=None):
    import tempfile
    d = tempfile.mkdtemp()
    A = [("UFL1", "MKMKMKMK"), ("UFC1", "ACDEFGHIK")]
    B = [("DDRGK1", "WYWYWYWY")]
    jobs = yaml_build.build_grid(A, B, os.path.join(d, "g"))
    assert len(jobs) == 2
    txt = open(jobs[0][1]).read()
    assert "id: A" in txt and "id: B" in txt and "version: 1" in txt
    cj = yaml_build.build_complex(A, B, "tri", os.path.join(d, "c"))
    ct = open(cj[0][1]).read()
    assert ct.count("protein:") == 3, "complex = 3 chains (2 A + 1 B)"


def test_fetch_parsers():
    label, seq = fetch.resolve_one("UFL1=MKMKMKMKMKMKMKMKMKMKMK")  # inline seq
    assert label == "UFL1" and seq.startswith("MK")
    assert fetch.TAIR_LOCUS.match("AT1G01010") and not fetch.TAIR_LOCUS.match("P69905")
    assert fetch.UNIPROT_ACC.match("P69905") and fetch.UNIPROT_ACC.match("Q9FPE7")


if __name__ == "__main__":
    test_peak_per_chainpair(); test_peak_definition()
    test_yaml_grid_and_complex(); test_fetch_parsers()
    print("all tests passed")
