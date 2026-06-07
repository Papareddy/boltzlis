"""boltzlis -- Boltz-2 protein-protein interface pipeline with LIS/actifpTM/PEAK metrics.

  fetch       protein IDs (UniProt acc / TAIR locus / raw seq) -> sequences
  yaml_build  chain-sets -> Boltz YAML(s) (grid screen or single multi-chain complex)
  collect     run AFM-LIS lis.py + PEAK -> ranked per-pair metric table
  structure   PDB/mmCIF parse + PEAK (chain split); actifpTM now comes from lis.py
  lis         vendored AFM-LIS CLI (github.com/flyark/AFM-LIS)
"""
__version__ = "0.1.0"
