#!/usr/bin/env python3
"""
find_contig_length.py

Walk a directory of genome subdirectories, extract contig lengths from each
sample's FASTA, and write a combined TSV suitable for circos_from_hmm.py.

The TSV schema is: sample<TAB>contig<TAB>length

usage:
python /home/juneq/hmm/scripts/circos/find_contig_length.py \
--genomes-dir /home/juneq/hmm/validation_genomes/T12D_validated/all \
--out /home/juneq/hmm/validation_genomes/T12D_validated/all/contig_lengths.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

FIELDNAMES = ["sample", "contig", "length"]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Generate contig length table from genome FASTA files."
    )
    ap.add_argument(
        "--genomes-dir",
        required=True,
        help="Root directory containing one subdirectory per genome sample.",
    )
    ap.add_argument(
        "--fasta-glob",
        default="*.fna",
        help="Glob pattern (relative to each genome directory) for FASTA files. Default: %(default)s",
    )
    ap.add_argument(
        "--out",
        default="contig_lengths.tsv",
        help="Output TSV path (default: %(default)s)",
    )
    return ap.parse_args()


def read_fai(fai_path: Path) -> Iterable[Tuple[str, int]]:
    with fai_path.open("r") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            contig, length = parts[0], int(parts[1])
            yield contig, length


def read_fasta_lengths(fasta_path: Path) -> Iterable[Tuple[str, int]]:
    lengths: Dict[str, int] = {}
    cur = None
    with fasta_path.open("r") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    continue
                contig = header.split()[0]
                lengths.setdefault(contig, 0)
                cur = contig
            else:
                if cur is None:
                    continue
                lengths[cur] += len(line.strip())
    for contig, length in lengths.items():
        yield contig, length


def gather_contig_lengths(genomes_dir: Path, fasta_glob: str) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    top_level_fastas = [p for p in sorted(genomes_dir.glob(fasta_glob)) if p.is_file()]
    genome_dirs = [p for p in sorted(genomes_dir.iterdir()) if p.is_dir()]
    if not genome_dirs and not top_level_fastas:
        print(f"[err] no FASTA files or genome subdirectories found in {genomes_dir}", file=sys.stderr)
        return rows

    def _process_sample(sample: str, fasta: Path) -> None:
        fai = fasta.with_suffix(fasta.suffix + ".fai")

        if fai.exists():
            entries = list(read_fai(fai))
        else:
            print(f"[info] {fai.name} missing; parsing FASTA to compute lengths for {sample}", file=sys.stderr)
            entries = list(read_fasta_lengths(fasta))

        if not entries:
            print(f"[warn] no contigs found for sample {sample} (file {fasta})", file=sys.stderr)
            return

        for contig, length in entries:
            rows.append({"sample": sample, "contig": contig, "length": length})

    # Case 1: FASTAs placed directly under genomes_dir
    for fasta in top_level_fastas:
        sample = fasta.stem
        _process_sample(sample, fasta)

    # Case 2: One FASTA per subdirectory under genomes_dir
    for gdir in genome_dirs:
        fastas = sorted(gdir.glob(fasta_glob))
        if not fastas:
            print(f"[warn] no FASTA matching '{fasta_glob}' in {gdir}", file=sys.stderr)
            continue
        if len(fastas) > 1:
            print(f"[warn] multiple FASTA files in {gdir}; using first: {fastas[0].name}", file=sys.stderr)
        _process_sample(gdir.name, fastas[0])

    return rows


def main() -> None:
    args = parse_args()
    genomes_dir = Path(args.genomes_dir)
    if not genomes_dir.is_dir():
        print(f"[err] genomes directory not found: {genomes_dir}", file=sys.stderr)
        sys.exit(1)

    rows = gather_contig_lengths(genomes_dir, args.fasta_glob)
    if not rows:
        print("[err] no contig lengths gathered; no output written", file=sys.stderr)
        sys.exit(1)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[ok] wrote {len(rows)} contig entries to {out_path}")


if __name__ == "__main__":
    main()
