import argparse
import shutil
import subprocess
from pathlib import Path


DEFAULT_PATHWAYS = "00642,00623,00622,00362"


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Generate KEGG pathway visualization links from BTEX hmmscan summary output."
        )
    )
    p.add_argument(
        "--hmmscan",
        required=True,
        help="Path to hmmscan summary CSV (must contain sample, hmm, and hits columns).",
    )
    p.add_argument(
        "-s",
        "--sample",
        default="all",
        help="Sample name from hmmscan CSV, or 'all' (default: all).",
    )
    p.add_argument(
        "--ko-map",
        dest="ko_map",
        required=True,
        help="Path to KO mapping TSV.",
    )
    p.add_argument(
        "--pathways",
        default=DEFAULT_PATHWAYS,
        help=f"Comma-separated KEGG pathway IDs (default: {DEFAULT_PATHWAYS}).",
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

    rscript = shutil.which("Rscript")
    if rscript is None:
        raise SystemExit("Rscript not found in PATH. Please install R and ensure Rscript is available.")

    repo_root = Path(__file__).resolve().parents[1]
    vis_script = repo_root / "visualization_scripts" / "visualize_pathways.R"
    if not vis_script.exists():
        raise SystemExit(f"Visualization script not found: {vis_script}")

    hmmscan = Path(args.hmmscan).expanduser().resolve()
    ko_map = Path(args.ko_map).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()

    if not hmmscan.exists():
        raise SystemExit(f"--hmmscan file not found: {hmmscan}")
    if not ko_map.exists():
        raise SystemExit(f"--ko-map file not found: {ko_map}")

    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        rscript,
        str(vis_script),
        "--hmmscan",
        str(hmmscan),
        "--sample",
        str(args.sample),
        "--ko-map",
        str(ko_map),
        "--pathways",
        str(args.pathways),
        "--outdir",
        str(outdir),
    ]
    subprocess.run(cmd, check=True)

