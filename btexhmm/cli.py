# btexhmm/cli.py
import argparse
from pathlib import Path

from .hmmsearch import annotate_proteins

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
DEFAULT_HMM_DIR = DATA_DIR / "flat"               # your packaged HMMs
DEFAULT_CUTOFFS = DATA_DIR / "hmm_cutoffs.tsv"    # your packaged cutoffs


def parse_args():
    p = argparse.ArgumentParser(
        description="Annotate a directory of protein FASTA files with the BTEX HMM database"
    )
    p.add_argument(
        "--proteins-dir",
        required=True,
        help=(
            "Directory containing protein FASTA files "
            "(for example *.faa, *.fa, *.fasta). Each file is treated as one sample."
        ),
    )
    p.add_argument(
        "--outdir",
        required=True,
        help="Output directory for results",
    )
    p.add_argument(
        "--cpus",
        type=int,
        default=8,
        help="Number of CPUs for hmmsearch",
    )
    return p.parse_args()


def main():
    args = parse_args()

    proteins_dir = Path(args.proteins_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    out_csv = outdir / "btex_hmm_summary.csv"
    annotate_proteins(
        proteins_path=proteins_dir,
        out_csv=out_csv,
        hmm_dir=DEFAULT_HMM_DIR,
        cutoffs_path=DEFAULT_CUTOFFS,
        cpus=args.cpus,
        evalue=None,
        domevalue=None,
    )
    print(f"wrote {out_csv}")
