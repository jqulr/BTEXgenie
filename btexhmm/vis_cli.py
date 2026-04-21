import argparse
import shlex
import shutil
import subprocess
from pathlib import Path

try:
    from .logging_utils import command_logger, run_logged_command
except ImportError:
    from logging_utils import command_logger, run_logged_command


DEFAULT_PATHWAYS = "00642,00623,00622,00362"
DEFAULT_KO_MAP = Path(__file__).resolve().parent / "data" / "btex_to_KO_map.tsv"


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Generate KEGG pathway visualization links from BTEX hmmscan summary output."
        )
    )
    p.add_argument(
        "--hmmscan",
        default=None,
        help="Path to hmmscan output CSV (must contain at least sample and hmm columns).",
    )
    p.add_argument(
        "-g",
        "--genome-dir",
        dest="genome_dir",
        metavar="GENOMES",
        default=None,
        help="Optional genome directory with per-genome subdirectories containing kofam_abv_thres.tsv.",
    )
    p.add_argument(
        "-s",
        "--sample",
        default="all",
        help="Sample name from hmmscan CSV, or 'all' (default: all).",
    )
    p.add_argument(
        "--pathways",
        default=DEFAULT_PATHWAYS,
        help=(
        f"Comma-separated KEGG pathway IDs (default: {DEFAULT_PATHWAYS}). "
        "Note: Supplying IDs outside of BTEX degradation while using "
        "BTEX-HMMs scan may result in limited annotations."
    ),
    )
    p.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory (required).",
    )
    return p.parse_args()


def main():
    args = parse_args()

    if bool(args.hmmscan) == bool(args.genome_dir):
        raise SystemExit("Provide exactly one of --hmmscan or --genome-dir.")

    rscript = shutil.which("Rscript")
    if rscript is None:
        raise SystemExit("Rscript not found in PATH. Please install R and ensure Rscript is available.")

    repo_root = Path(__file__).resolve().parents[1]
    vis_script = repo_root / "visualization_scripts" / "visualize_pathways.R"
    if not vis_script.exists():
        raise SystemExit(f"Visualization script not found: {vis_script}")

    ko_map = DEFAULT_KO_MAP.resolve()
    outdir = Path(args.outdir).expanduser().resolve()

    if not ko_map.exists():
        raise SystemExit(f"Default KO map file not found: {ko_map}")

    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        rscript,
        str(vis_script),
        "--sample",
        str(args.sample),
        "--ko-map",
        str(ko_map),
        "--pathways",
        str(args.pathways),
        "--outdir",
        str(outdir),
    ]
    if args.hmmscan:
        hmmscan = Path(args.hmmscan).expanduser().resolve()
        if not hmmscan.exists():
            raise SystemExit(f"--hmmscan file not found: {hmmscan}")
        cmd[2:2] = ["--hmmscan", str(hmmscan)]
    if args.genome_dir:
        genome_dir = Path(args.genome_dir).expanduser().resolve()
        if not genome_dir.exists():
            raise SystemExit(f"--genome-dir not found: {genome_dir}")
        cmd[2:2] = ["--genome-dir", str(genome_dir)]
    log_path = outdir / "log_file_vis-btex.txt"
    try:
        with command_logger(log_path):
            print(f"[info] writing vis-btex log to {log_path}")
            print(f"[cmd] {shlex.join(['vis-btex', *[str(part) for part in cmd[2:]]])}")
            run_logged_command(cmd)
    except subprocess.CalledProcessError as exc:
        raise SystemExit(exc.returncode)
