"""CLI wrapper for the Circos visualization entry point."""

import argparse

from visualization_scripts.circos import run


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
    ap.add_argument("--outdir", required=True, help="Output root dir")
    ap.add_argument("-s", "--sample", required=True, help="Sample/genome name to render")
    ap.add_argument(
        "--only-hit-contigs",
        action="store_true",
        help="Include only contigs that have >=1 hit",
    )
    ap.add_argument(
        "--operon",
        action="store_true",
        help="Switch to operon-centric plotting mode with layered tracks for completeness.",
    )
    ap.add_argument(
        "--dna",
        help=(
            "Path to the genome nucleotide FASTA (.fna/.fasta). "
            "Used to generate contig lengths when --contig-lengths is omitted and to populate GenBank sequences."
        ),
    )
    return ap


def parse_args(argv=None):
    return build_parser().parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    return run(args)


__all__ = ["build_parser", "parse_args", "main"]
