import argparse
import shlex
from pathlib import Path

from .hmmscan import main as hmmscan_main
from .logging_utils import command_logger

HERE = Path(__file__).resolve().parent
DEFAULT_HMM_LIB = HERE / "hmms" / "all_models.hmm"


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Annotate a directory of genome DNA FASTA files or protein FASTA files "
            "with the BTEX HMM database. DNA inputs are first gene-called with Prodigal."
        )
    )
    p.add_argument(
        "-g",
        "--genome-dir",
        dest="genome_dir",
        metavar="GENOMES",
        required=True,
        help=(
            "Directory containing genome DNA FASTA files or protein FASTA files "
            "(for example *.fna, *.fa, *.fasta, *.faa)."
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
    prodigal_group = p.add_mutually_exclusive_group()
    prodigal_group.add_argument(
        "--meta",
        dest="prodigal_mode",
        action="store_const",
        const="meta",
        help="Use Prodigal meta mode for DNA genome inputs",
    )
    prodigal_group.add_argument(
        "--single",
        dest="prodigal_mode",
        action="store_const",
        const="single",
        help="Use Prodigal single mode for DNA genome inputs",
    )
    p.add_argument(
        "--evalue",
        default=None,
        help="Sequence E-value cutoff passed to hmmscan (e.g. 1e-5)",
    )
    p.add_argument(
        "--skip-kofam",
        action="store_true",
        help="Skip the KOfam HMM database download and search after hmmscan completes",
    )
    return p.parse_args()


def main():
    args = parse_args()

    genomes = Path(args.genome_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    out_csv = outdir / "btex_hmm_summary.csv"
    hmmscan_argv = [
        "--hmm-lib", str(DEFAULT_HMM_LIB),
        "--genomes-dir", str(genomes),
        "--out", str(out_csv),
        "--cpus", str(args.cpus),
    ]
    if args.prodigal_mode == "meta":
        hmmscan_argv.append("-meta")
    elif args.prodigal_mode == "single":
        hmmscan_argv.append("-single")
    if args.evalue:
        hmmscan_argv.extend(["--evalue", str(args.evalue)])
    if args.skip_kofam:
        hmmscan_argv.append("--skip-kofam")

    top_cmd = [
        "annotate-btex",
        "-g",
        str(genomes),
        "-o",
        str(outdir),
        "--cpus",
        str(args.cpus),
    ]
    if args.prodigal_mode == "meta":
        top_cmd.append("--meta")
    elif args.prodigal_mode == "single":
        top_cmd.append("--single")
    if args.evalue:
        top_cmd.extend(["--evalue", str(args.evalue)])
    if args.skip_kofam:
        top_cmd.append("--skip-kofam")

    log_path = outdir / "log_file_annotate-btex.txt"
    with command_logger(log_path):
        print(f"[info] writing annotate-btex log to {log_path}")
        print(f"[cmd] {shlex.join(top_cmd)}")
        hmmscan_main(hmmscan_argv)
