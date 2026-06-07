#!/bin/bash
# Example: screen UFMylation E3/E2 against the ER adaptor + C53, from a local FASTA.
set -euo pipefail
cd "$(dirname "$0")/.."
python pipeline.py --name UFM_screen \
    --chainA UFL1,UFC1 \
    --chainB DDRGK1,CDK5RAP3 \
    --fasta examples/ufm_machinery.fasta \
    --outdir runs/UFM_screen \
    --diffusion-samples 5
echo
echo "Built YAMLs + submit.sbatch. After the cluster run finishes:"
echo "  python pipeline.py --collect runs/UFM_screen/out -o runs/UFM_screen/metrics.tsv --config config.yaml"
