"""CLI wrapper for the Circos visualization entry point."""

import argparse
import shutil
import subprocess
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Generate Circos plot files for HMM hits on a single genome."
    )
    ap.add_argument(
        "--hmmscan",
        required=True,
        help="CSV (btex_hmm_summary.csv): sample,hmm,hits,total_genes,hit_headers from hmmsearch.py.",
    )
    ap.add_argument(
        "--contig-lengths",
        help=(
            "TSV: sample<TAB>contig<TAB>length. "
            "If omitted, contig_length.tsv is generated under the Circos output dir using --dna."
        ),
    )
    ap.add_argument("-o", "--outdir", required=True, help="Output dir")
    ap.add_argument("-s", "--sample", required=True, help="Sample/genome name to plot (must match sample column in --hmmscan)")
    ap.add_argument(
        "--dna",
        help=(
            "Path to the genome nucleotide FASTA (.fna/.fasta). "
            "Used to generate contig lengths when --contig-lengths is omitted."
        ),
    )
    return ap


def parse_args(argv=None):
    return build_parser().parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise SystemExit("Rscript not found in PATH.")

    repo_root = Path(__file__).resolve().parents[1]
    r_circos = repo_root / "visualization_scripts" / "circos.R"
    if not r_circos.exists():
        raise SystemExit(f"Circos R script not found: {r_circos}")

    cmd = [
        rscript,
        str(r_circos),
        "--hmmscan", args.hmmscan,
        "--outdir", args.outdir,
        "--sample", args.sample,
    ]
    if args.contig_lengths:
        cmd += ["--contig-lengths", args.contig_lengths]
    if args.dna:
        cmd += ["--dna", args.dna]

    subprocess.run(cmd, check=True)

    return 0


__all__ = ["build_parser", "parse_args", "main"]
