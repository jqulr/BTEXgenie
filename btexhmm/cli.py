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
            "Annotate a directory or single genome DNA FASTA files or protein FASTA files "
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
            "Directory or single file containing genome DNA FASTA or protein FASTA input "
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
        help="Number of CPUs for hmmscan (default: 8)",
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
        help="Use Prodigal single mode for DNA genome inputs (default mode)",
    )
    p.add_argument(
        "--evalue",
        default="1e-5",
        help="Sequence E-value cutoff passed to hmmscan (default: 1e-5)",
    )
    p.add_argument(
        "--kofam",
        action="store_true",
        help="Run the HMM search against the KOfam database in addition to the BTEX-HMM database.",
    )
    return p.parse_args()


def main():
    args = parse_args()
    prodigal_mode = args.prodigal_mode or "single"

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
    if prodigal_mode == "meta":
        hmmscan_argv.append("-meta")
    elif prodigal_mode == "single":
        hmmscan_argv.append("-single")
    if args.evalue:
        hmmscan_argv.extend(["--evalue", str(args.evalue)])
    if args.kofam:
        hmmscan_argv.append("--kofam")

    top_cmd = [
        "btex-annotate",
        "-g",
        str(genomes),
        "-o",
        str(outdir),
        "--cpus",
        str(args.cpus),
    ]
    if prodigal_mode == "meta":
        top_cmd.append("--meta")
    elif prodigal_mode == "single":
        top_cmd.append("--single")
    if args.evalue:
        top_cmd.extend(["--evalue", str(args.evalue)])
    if args.kofam:
        top_cmd.append("--kofam")

    log_path = outdir / "btex_annotate.log"
    with command_logger(log_path):
        print(f"[info] writing btex-annotate log to {log_path}")
        print(f"[cmd] {shlex.join(top_cmd)}")
        print("[info] running BTEXgenie on input genomes")
        hmmscan_main(hmmscan_argv)
        print("Done!")
