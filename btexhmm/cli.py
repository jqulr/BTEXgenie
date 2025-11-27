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
        description=(
            "Annotate a directory of protein FASTA files with the BTEX HMM database. "
            "If only nucleotide FASTAs are present, Prodigal can be run first to generate proteins."
        )
    )
    p.add_argument(
        "--proteins-dir",
        required=True,
        help=(
            "Directory containing protein or nucleotide FASTA files "
            "(for example *.faa, *.fa, *.fasta, *.fna). Each file is treated as one sample."
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
    p.add_argument(
        "--prot-glob",
        default="*proteins.faa,*.faa,*.fa,*.fasta,*.fna",
        help="Comma-separated glob(s) for input FASTA files (default includes protein and nucleotide FASTAs)",
    )
    p.add_argument(
        "--run-prodigal",
        action="store_true",
        help="Run Prodigal on nucleotide FASTAs to create proteins before hmmsearch (auto-enabled when only nucleotide FASTAs are found)",
    )
    p.add_argument(
        "--prodigal-out-dir",
        default=None,
        help="Directory for Prodigal protein FASTAs (default: <proteins-dir>/prodigal_proteins)",
    )
    p.add_argument(
        "--prodigal-mode",
        default="single",
        choices=["single", "meta"],
        help="Prodigal mode passed to -p (default: single)",
    )
    p.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip rerunning hmmsearch/prodigal if outputs already exist",
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
        prot_glob=tuple(p.strip() for p in args.prot_glob.split(",") if p.strip()),
        run_prodigal=args.run_prodigal,
        prodigal_out_dir=Path(args.prodigal_out_dir) if args.prodigal_out_dir else None,
        prodigal_mode=args.prodigal_mode,
        skip_existing=args.skip_existing,
    )
    print(f"wrote {out_csv}")
