# boltzlis — Boltz-2 interface screening with LIS / actifpTM / PEAK

Give **chain A** and **chain B** as protein IDs (each one *or many*), get Boltz-2
input YAMLs + a ready-to-submit SLURM array, and — after the run — a single ranked
table with **every interface metric**:

| metric | what it is |
|---|---|
| **iLIS** = √(LIS·cLIS) | AFM-LIS primary score; interaction call ≈ **iLIS ≥ 0.223** |
| LIS / cLIS | local interaction score (all / contact-restricted) |
| LIA / cLIA | local interaction *area* (interface size) |
| **actifpTM** | interface-restricted ipTM (flank-robust) |
| ipSAE | Dunbrack interface score |
| **PEAK** | 1 − min(inter-chain PAE)/30 |
| ipTM / pTM / pLDDT | native Boltz confidence |

Metrics come from vendored [AFM-LIS](https://github.com/flyark/AFM-LIS) `lis.py`
(parallel `-w`) plus a PEAK pass, aggregated to mean+max over the diffusion models.

> **What this is.** `boltzlis` is a thin **orchestration wrapper** — built by
> **Ranjith Papareddy** — that runs **Boltz-2** for structure prediction and computes
> *published* interface-confidence metrics (**LIS / iLIS**, **actifpTM**, **ipSAE**, **PEAK**).
> The modelling and the metrics are other people's science; boltzlis just automates the
> boring part end-to-end — IDs → sequences → folds → all metrics → one ranked table, across a
> cluster. **If you use it, please cite the underlying methods** (see
> [Credits & citation](#credits--citation)).

---

## Step-by-step guide (local + Helix / bwForCluster)

**Who runs what, where** — the work splits across three contexts:

| step | runs on | needs |
|---|---|---|
| 1. build inputs (resolve IDs → sequences → YAMLs) | **your laptop** (or any internet machine) | internet (UniProt), Python |
| 2. Boltz folds (GPU) | **Helix GPU node** | Boltz install + GPU |
| 3. collect metrics | **Helix** (or laptop, if outputs copied) | numpy + scipy + pandas |
| 4. plot / inspect | laptop | matplotlib |

> A coding agent (e.g. Claude Code) can do all four end-to-end: give it this repo +
> SSH access to Helix and it can install, build, submit, collect, and pull results.

### 0. One-time install — laptop
```bash
git clone https://github.com/Papareddy/boltzlis.git
cd boltzlis
python -m venv .venv && source .venv/bin/activate     # or: conda create -n boltzlis python=3.11
pip install -r requirements.txt                        # numpy scipy pandas matplotlib
python tests/test_metrics.py                            # -> "all tests passed"
```

### 1. One-time setup — Helix
```bash
ssh helix                                               # your bwForCluster login
# (a) a workspace to hold runs (60-day, extendable):
ws_allocate boltzlis 60        # -> prints a path; call it $WS
export WS=$(ws_find boltzlis)
# (b) clone the repo on Helix too (needed for collect + lis.py):
git clone https://github.com/Papareddy/boltzlis.git $WS/boltzlis_repo
# (c) a Boltz GPU env (one-time; ~minutes). Conda recommended:
conda create -y -n boltz python=3.11 && conda activate boltz
pip install boltz                                       # Boltz-2
#     model weights download to --cache on first run (needs internet once;
#     run one tiny job from a login/internet node to pre-fill the cache).
# (d) an analysis env for collect/lis (numpy+scipy+pandas). If your boltz env
#     already has them, reuse it; else: conda create -y -n analysis numpy scipy pandas
# (e) write your config (gitignored; never commit real paths):
cd $WS/boltzlis_repo
cp config.example.yaml config.yaml
#     edit config.yaml ->  env: <path to boltz env>   cache: <cache dir>
#                          python: <analysis python with numpy/scipy/pandas>
#                          partition/gpumem/time/mem/cpus: your GPU partition
```

### 2. Build inputs — laptop
```bash
python pipeline.py --name UFM_screen \
    --chainA UFL1,UFC1 \
    --chainB DDRGK1,CDK5RAP3 \
    --fasta examples/ufm_machinery.fasta \
    --outdir runs/UFM_screen --diffusion-samples 5
# writes runs/UFM_screen/{yaml/*.yaml, submit.sbatch, yaml_list.txt, RUNBOOK.md}
```
- IDs = a **name in your `--fasta`** file (above), a **UniProt accession** (`P61960`),
  a **TAIR locus** (`AT1G01010`), or **`LABEL=SEQUENCE`**. `LABEL=` is just a display name.
  Bundling a local FASTA keeps your targets off any web lookup.
- `--chainA` / `--chainB` are comma-separated lists → "single or multiple" scales the
  grid: `1×1`, `N×1`, or `N×M` pairwise 2-chain folds.
- One multi-chain assembly instead of a grid: add `--mode complex`.
- No internet / non-UniProt sequences: pass `--fasta my_seqs.fa` and use those names as IDs.

### 3. Sync + submit — Helix
```bash
rsync -av runs/UFM_screen  helix:$WS/runs/                # copy the built run dir over
ssh helix
cd $WS/boltzlis_repo
N=$(wc -l < $WS/runs/UFM_screen/yaml_list.txt)
sbatch --array=1-$N $WS/runs/UFM_screen/submit.sbatch     # one GPU fold per pair
squeue --me                                              # watch; logs in runs/.../logs/
```
> **MSA note (important):** Boltz needs an MSA. The template uses `--use_msa_server`
> (in `config.yaml: boltz_extra`), which auto-generates via the ColabFold server and
> **needs internet on the run node**. If your GPU nodes are offline, either run from an
> internet-capable node, or precompute MSAs and add an `msa:` path per YAML and drop
> `--use_msa_server` from `boltz_extra`.

### 4. Collect every metric — Helix (one command)
```bash
# after the array finishes:
$(grep '^python:' config.yaml | awk '{print $2}') -m boltzlis.collect \
    $WS/runs/UFM_screen/out -o $WS/runs/UFM_screen/metrics.tsv -w 8
# (or simply, from the repo dir with the analysis env active:)
python -m boltzlis.collect $WS/runs/UFM_screen/out -o $WS/runs/UFM_screen/metrics.tsv -w 8
```
`metrics.tsv` = one row per pair, ranked by `iLIS_max`, with **mean+max of every metric**
(iLIS/LIS/cLIS/LIA/ipSAE/actifpTM/ipTM/pTM/PEAK/pLDDT).

### 5. Pull + inspect — laptop
```bash
rsync -av helix:$WS/runs/UFM_screen/metrics.tsv  runs/UFM_screen/
column -t -s$'\t' runs/UFM_screen/metrics.tsv | less -S
```
Rule of thumb: a confident interaction passes **iLIS ≥ 0.22 and PEAK ≥ 0.7**.

## Layout
```
pipeline.py            main CLI (build / submit-recipe / collect)
boltzlis/
  fetch.py             ID (UniProt acc | TAIR locus | seq) -> sequence (+cache)
  yaml_build.py        chain-sets -> Boltz YAML(s)  (grid | complex)
  collect.py           lis.py + PEAK -> ranked per-pair metric table   [python -m boltzlis.collect]
  structure.py         PDB/mmCIF parser + PEAK chain split (actifpTM comes from lis.py, not here)
  lis.py               vendored AFM-LIS (iLIS/LIS/cLIS/LIA/ipSAE/actifpTM)
slurm/boltz.sbatch.tmpl  templated GPU array job
config.example.yaml    cluster paths/resources
tests/                 local unit tests (no cluster/GPU needed)
```

## Notes
- Sequence fetch runs **locally** (compute nodes have no internet); Boltz GPU folds
  run on the cluster; metric collection runs anywhere with numpy+scipy.
- iLIS 0.223 is AF-Multimer/Y2H-calibrated — treat as approximate on Boltz.
- `collect.py` ingests Boltz output dirs (`<pair>/boltz_results_*/predictions/...`)
  as-is via lis.py's native Boltz adapter.

## Credits & citation

`boltzlis` (this wrapper) is by **Ranjith Papareddy**, MIT-licensed. It does **not** introduce
a new method — it orchestrates and scores the tools below. **If `boltzlis` is useful in your
work, cite the underlying methods** (and a link to this repo is appreciated):

| component | what it does here | cite |
|---|---|---|
| **Boltz-2** | structure prediction (the folds) | Passaro *et al.* 2025, bioRxiv [10.1101/2025.06.14.659707](https://doi.org/10.1101/2025.06.14.659707); Boltz-1: Wohlwend *et al.* 2024 [10.1101/2024.11.19.624167](https://doi.org/10.1101/2024.11.19.624167) |
| **iLIS / LIS / cLIS / LIA** | local interaction scores (vendored `lis.py`) | iLIS: Kim *et al.* 2026 [10.64898/2026.04.14.718529](https://doi.org/10.64898/2026.04.14.718529); LIS+LIA: Kim *et al.* 2024, bioRxiv [10.1101/2024.02.19.580970](https://doi.org/10.1101/2024.02.19.580970) — repo: [flyark/AFM-LIS](https://github.com/flyark/AFM-LIS) |
| **actifpTM** | interface-restricted ipTM | Varga, Ovchinnikov & Schueler-Furman 2025, *Bioinformatics* [10.1093/bioinformatics/btaf107](https://doi.org/10.1093/bioinformatics/btaf107) |
| **ipSAE** | aligned-error interface score | Dunbrack 2025, bioRxiv [10.1101/2025.02.10.637595](https://doi.org/10.1101/2025.02.10.637595) |
| **ipTM / pTM** | native confidence (AlphaFold-Multimer metric, reported by Boltz) | Evans *et al.* 2021 [10.1101/2021.10.04.463034](https://doi.org/10.1101/2021.10.04.463034) |

`PEAK` ( = 1 − min-interchain-PAE/30 ) is a convenience scalar defined in this repo; no separate citation.

<details><summary>BibTeX</summary>

> DOIs and author-years are verified from source. The two Kim *et al.* (LIS/iLIS) **titles**
> below are abbreviated — confirm the exact title/author list at each DOI before citing.

```bibtex
@article{passaro2025boltz2, title={Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction}, author={Passaro, Saro and Corso, Gabriele and Wohlwend, Jeremy and others}, journal={bioRxiv}, year={2025}, doi={10.1101/2025.06.14.659707}}
@article{wohlwend2024boltz1, title={Boltz-1: Democratizing Biomolecular Interaction Modeling}, author={Wohlwend, Jeremy and Corso, Gabriele and Passaro, Saro and others}, journal={bioRxiv}, year={2024}, doi={10.1101/2024.11.19.624167}}
@article{kim2026ilis, title={Integrated Local Interaction Score (iLIS)}, author={Kim, Ah-Ram and others}, year={2026}, doi={10.64898/2026.04.14.718529}}
@article{kim2024lis, title={Enhancing Protein-Protein Interaction Prediction with Local Interaction Score from AlphaFold-Multimer}, author={Kim, Ah-Ram and others}, journal={bioRxiv}, year={2024}, doi={10.1101/2024.02.19.580970}}
@article{varga2025actifptm, title={actifpTM: a refined confidence metric of AlphaFold2 predictions involving flexible regions}, author={Varga, Julia K. and Ovchinnikov, Sergey and Schueler-Furman, Ora}, journal={Bioinformatics}, year={2025}, doi={10.1093/bioinformatics/btaf107}}
@article{dunbrack2025ipsae, title={Res ipSAE loquunt: What's wrong with AlphaFold's ipTM score and how to fix it}, author={Dunbrack, Roland L.}, journal={bioRxiv}, year={2025}, doi={10.1101/2025.02.10.637595}}
```
</details>
