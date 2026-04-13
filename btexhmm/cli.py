# btexhmm/cli.py
import argparse
from pathlib import Path

from .hmmscan import main as hmmscan_main

HERE = Path(__file__).resolve().parent
DEFAULT_HMM_LIB = HERE / "hmms" / "all_models.hmm"


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Annotate a directory of protein FASTA files with the BTEX HMM database. "
            "Runs hmmscan using the bundled BTEX HMM library."
        )
    )
    p.add_argument(
        "-p",
        "--proteins",
        dest="proteins",
        required=True,
        help=(
            "Directory containing protein FASTA files, or a single protein FASTA file "
            "(for example *.faa, *.fa, *.fasta)."
        ),
    )
    p.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory for results",
    )
    p.add_argument(
        "--cpus",
        type=int,
        default=8,
        help="Number of CPUs for hmmscan",
    )
    p.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip rerunning hmmscan if outputs already exist",
    )
    p.add_argument(
        "--evalue",
        default=None,
        help="Sequence E-value cutoff passed to hmmscan (e.g. 1e-5)",
    )
    p.add_argument(
        "--domevalue",
        default=None,
        help="Domain E-value cutoff passed to hmmscan (e.g. 1e-5)",
    )
    p.add_argument(
        "--min-bits",
        type=float,
        default=None,
        help="Optional minimum bit score when counting hits",
    )
    p.add_argument(
        "--max-ie",
        type=float,
        default=None,
        help="Optional maximum i-Evalue when counting hits",
    )
    p.add_argument(
        "--unique-per-seq",
        action="store_true",
        help="Count at most one hit per sequence per HMM",
    )
    p.add_argument(
        "--all-hits-per-protein",
        action="store_true",
        help="Count all passing HMM hits per protein (default is best hit per protein)",
    )
    p.add_argument(
        "--models",
        default=None,
        help="Optional comma/space-separated subset of HMM names to report",
    )
    return p.parse_args()


def main():
    args = parse_args()

    proteins = Path(args.proteins).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    out_csv = outdir / "btex_hmm_summary.csv"
    hmmscan_argv = [
        "--hmm-lib", str(DEFAULT_HMM_LIB),
        "--proteins", str(proteins),
        "--out", str(out_csv),
        "--cpus", str(args.cpus),
    ]
    if args.skip_existing:
        hmmscan_argv.append("--skip-existing")
    if args.evalue:
        hmmscan_argv.extend(["--evalue", str(args.evalue)])
    if args.domevalue:
        hmmscan_argv.extend(["--domevalue", str(args.domevalue)])
    if args.min_bits is not None:
        hmmscan_argv.extend(["--min-bits", str(args.min_bits)])
    if args.max_ie is not None:
        hmmscan_argv.extend(["--max-ie", str(args.max_ie)])
    if args.unique_per_seq:
        hmmscan_argv.append("--unique-per-seq")
    if args.all_hits_per_protein:
        hmmscan_argv.append("--all-hits-per-protein")
    if args.models:
        hmmscan_argv.extend(["--models", args.models])

    hmmscan_main(hmmscan_argv)
