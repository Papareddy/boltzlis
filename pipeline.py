#!/usr/bin/env python3
"""boltzlis main pipeline -- give chain A and chain B (each one or many protein IDs),
get Boltz YAMLs + a ready-to-submit SLURM array + the post-run metric-collect command.

Examples
--------
# screen: every A x every B as 2-chain folds (IDs from a local FASTA, or UniProt/TAIR)
python pipeline.py --name UFM_screen \
    --chainA UFL1,UFC1 --chainB DDRGK1,CDK5RAP3 \
    --fasta examples/ufm_machinery.fasta --outdir runs/UFM_screen

# single multi-chain complex (all A + all B chains in one structure)
python pipeline.py --name UFM_complex --mode complex \
    --chainA UFM1,UBA5 --chainB UFC1 --fasta examples/ufm_machinery.fasta \
    --outdir runs/UFM_complex

# just (re)collect metrics for a finished Boltz out dir
python pipeline.py --collect runs/UFM_screen/out -o runs/UFM_screen/metrics.tsv
"""
import os, sys, argparse, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from boltzlis import fetch, yaml_build
from boltzlis.collect import collect_all

HERE = os.path.dirname(os.path.abspath(__file__))
TMPL = os.path.join(HERE, "slurm", "boltz.sbatch.tmpl")


def load_config(path):
    cfg = dict(env="/path/to/boltz-env", cache="/path/to/boltz-cache",
               python="python3", partition="gpu-single", gpumem="80G",
               time="12:00:00", mem="120g", cpus="12", workdir="$PWD",
               boltz_extra="--use_msa_server --no_kernels --override")
    if path and os.path.exists(path):
        import re
        for ln in open(path):
            ln = ln.split("#", 1)[0].strip()
            if ":" in ln:
                k, _, v = ln.partition(":")
                cfg[k.strip()] = v.strip().strip('"').strip("'")
    return cfg


def render_sbatch(jobs, outdir, name, cfg, diffusion, output_format="cif"):
    listfile = os.path.join(outdir, "yaml_list.txt")
    open(listfile, "w").write("\n".join(p for _, p in jobs) + "\n")
    tmpl = open(TMPL).read()
    txt = tmpl.format(JOBNAME=name[:14], PARTITION=cfg["partition"], GPUMEM=cfg["gpumem"],
                      TIME=cfg["time"], MEM=cfg["mem"], CPUS=cfg["cpus"], NJOBS=len(jobs),
                      ENV=cfg["env"], CACHE=cfg["cache"], OUTDIR=os.path.abspath(outdir),
                      LISTFILE=os.path.abspath(listfile), DIFFUSION=diffusion,
                      OUTFMT=output_format,
                      BOLTZ_EXTRA=cfg.get("boltz_extra", "--use_msa_server --no_kernels --override"))
    sb = os.path.join(outdir, "submit.sbatch")
    open(sb, "w").write(txt)
    return sb, listfile


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", default="run")
    ap.add_argument("--chainA", help="comma-sep IDs (UniProt acc / TAIR locus / LABEL=SEQ)")
    ap.add_argument("--chainB", help="comma-sep IDs")
    ap.add_argument("--outdir", default="runs/run")
    ap.add_argument("--mode", choices=["grid", "complex"], default="grid",
                    help="grid: A x B pairwise 2-chain folds (screen). complex: one multi-chain fold.")
    ap.add_argument("--diffusion-samples", type=int, default=5)
    ap.add_argument("--fasta", action="append", default=[], help="FASTA file(s) overriding ID lookup")
    ap.add_argument("--cache", default="seq_cache", help="sequence cache dir")
    ap.add_argument("--config", default=os.path.join(HERE, "config.yaml"))
    ap.add_argument("--output-format", choices=["cif", "pdb"], default="cif")
    ap.add_argument("--collect", metavar="OUT_DIR", help="skip build; just collect metrics for this Boltz out dir")
    ap.add_argument("-o", "--output", default=None, help="metrics TSV (with --collect)")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()

    if a.collect:
        out = a.output or os.path.join(os.path.dirname(a.collect.rstrip("/")) or ".", "metrics.tsv")
        cfg = load_config(a.config)
        agg = collect_all(a.collect, out, python_exe=cfg.get("python"), workers=a.workers)
        print("wrote %s (%d pairs)" % (out, len(agg)))
        return

    if not (a.chainA and a.chainB):
        ap.error("--chainA and --chainB are required (unless --collect)")
    overrides = {}
    for f in a.fasta:
        overrides.update(fetch.load_fasta(f))
    os.makedirs(a.outdir, exist_ok=True)
    A = fetch.resolve_chain(a.chainA, a.cache, overrides)
    B = fetch.resolve_chain(a.chainB, a.cache, overrides)
    print("chain A: %s" % ", ".join("%s(%daa)" % (l, len(s)) for l, s in A))
    print("chain B: %s" % ", ".join("%s(%daa)" % (l, len(s)) for l, s in B))

    ydir = os.path.join(a.outdir, "yaml")
    if a.mode == "grid":
        jobs = yaml_build.build_grid(A, B, ydir)
    else:
        jobs = yaml_build.build_complex(A, B, a.name, ydir)
    print("built %d YAML(s) -> %s" % (len(jobs), ydir))

    cfg = load_config(a.config)
    sb, listfile = render_sbatch(jobs, a.outdir, a.name, cfg, a.diffusion_samples, a.output_format)
    runbook = os.path.join(a.outdir, "RUNBOOK.md")
    open(runbook, "w").write(
        "# %s run\n\n%d fold(s), mode=%s, diffusion_samples=%d.\n\n"
        "## On the cluster\n```bash\n"
        "# 1. sync this run dir to the workspace, then:\n"
        "sbatch --array=1-%d %s\n"
        "# 2. when the array finishes, collect ALL metrics:\n"
        "python %s/pipeline.py --collect %s/out -o %s/metrics.tsv --config <config.yaml>\n```\n"
        % (a.name, len(jobs), a.mode, a.diffusion_samples, len(jobs),
           os.path.basename(sb), HERE, os.path.abspath(a.outdir), os.path.abspath(a.outdir)))
    print("SLURM array  -> %s   (sbatch --array=1-%d %s)" % (sb, len(jobs), sb))
    print("runbook      -> %s" % runbook)
    print("\nNext: review YAMLs, sync %s to the cluster, then `sbatch --array=1-%d %s`."
          % (a.outdir, len(jobs), os.path.basename(sb)))


if __name__ == "__main__":
    main()
