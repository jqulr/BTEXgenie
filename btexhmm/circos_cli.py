"""CLI wrapper for the Circos visualization entry point."""

import argparse
import shlex
import shutil
from pathlib import Path

try:
    from .logging_utils import command_logger, run_logged_command
except ImportError:
    from logging_utils import command_logger, run_logged_command


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Generate Circos plot files for HMM hits on a single genome."
    )
    ap.add_argument(
        "--hmmscan",
        required=True,
        help="CSV (btex_hmm_summary.csv): hit-level hmmscan output with sample, hmm, and hit_header columns.",
    )
    ap.add_argument(
        "-o", "--outdir", required=True, help="Output dir"
    )
    ap.add_argument("-s", "--sample", required=True, help="Sample/genome name to plot (must match sample column in --hmmscan)")
    ap.add_argument(
        "-g",
        "--genome",
        required=True,
        help=(
            "Path to the genome nucleotide FASTA (.fna/.fasta). "
            "Used to generate contig lengths for the Circos plot."
        ),
    )
    ap.add_argument(
        "--prodigal-gbk",
        help="Optional Prodigal GenBank file for mapping KOfam hits to genomic coordinates.",
    )
    ap.add_argument(
        "--kofam-output",
        dest="kofam_output",
        help="Optional KOfam hit table to pass through to the Circos R script.",
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
    if args.genome:
        cmd += ["--genome", args.genome]
    if args.prodigal_gbk:
        cmd += ["--prodigal-gbk", args.prodigal_gbk]
    if args.kofam_output:
        cmd += ["--kofam", args.kofam_output]

    log_path = Path(args.outdir).expanduser().resolve() / "btex_run_circos.log"
    with command_logger(log_path):
        print(f"[info] writing btex-run-circos log to {log_path}")
        print(f"[cmd] {shlex.join(['btex-run-circos', *[str(part) for part in cmd[2:]]])}")
        run_logged_command(cmd)

    return 0


__all__ = ["build_parser", "parse_args", "main"]
