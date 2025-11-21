#!/usr/bin/env python3
"""
Run HMMER hmmsearch for a directory of HMMs against one or more protein FASTA files.

Features
  * Input HMMs: directory containing *.hmm
  * Target proteins: single FASTA file or a directory of FASTA files
  * For each protein FASTA, creates a subdirectory under --out-dir named after the FASTA stem
    and writes one domtblout per HMM:
        <out-dir>/<sample>/<hmm_name>.domtblout
  * Produces a summary CSV in --out-dir with columns:
        sample,hmm,hits,total_genes,hit_headers
    where:
        sample       = protein FASTA stem
        hmm          = HMM model stem
        hits         = number of matching sequences
        total_genes  = number of sequences in that FASTA
        hit_headers  = full FASTA header lines for matching sequences (joined by "|")

Example
  python /home/juneq/hmm/scripts/a-test/hmmsearch.py \
      --hmm-dir /home/juneq/hmm/archetypes/hmm_cutoffs/hmms/_packed/flat \
      --proteins /home/juneq/hmm/bioproject_PRJNA378566_Jiangsu/prodigal_out_unfiltered \
      --prot-glob "*.faa" \
      --out-dir /home/juneq/hmm/bioproject_PRJNA378566_Jiangsu/hmm_outputs_test_github
"""

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path
from collections import defaultdict


def check_bin(name: str) -> None:
    """Ensure an executable is available in PATH."""
    from shutil import which
    if which(name) is None:
        raise SystemExit(f"[err] '{name}' not found in PATH")


def sh(cmd: list[str]) -> None:
    """Run a command, raising on failure, silencing stdout/stderr."""
    print("[cmd]", " ".join(cmd), flush=True)
    subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def parse_domtbl(path: Path):
    """
    Parse an HMMER domtblout file.

    Returns
        list of (seqid, hmm_name, i_evalue, bit_score)
    """
    rows = []
    with open(path, "rt") as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line.strip())
            if len(parts) < 22:
                continue
            seqid = parts[0]        # target name
            hmm_name = parts[3]     # query name (model name)
            i_e = float(parts[12])  # i-Evalue
            bits = float(parts[13]) # bit score
            rows.append((seqid, hmm_name, i_e, bits))
    return rows


def load_headers_for_fasta(protein_faa: Path) -> dict[str, str]:
    """
    Build mapping from sequence id -> full header line (including '>').

    The sequence id is taken as the first whitespace-delimited token after '>'.
    """
    mp: dict[str, str] = {}
    with open(protein_faa, "rt", errors="ignore") as fh:
        for ln in fh:
            if not ln.startswith(">"):
                continue
            hdr = ln.rstrip("\n\r")
            seqid = hdr[1:].split()[0]
            mp.setdefault(seqid, hdr)
    return mp


def find_protein_fastas(proteins_path: Path, prot_glob: str) -> list[Path]:
    """
    Determine the set of protein FASTA files from the given path.

    If proteins_path is a file, returns [proteins_path].
    If proteins_path is a directory, returns all files matching prot_glob inside it.
    """
    if proteins_path.is_file():
        return [proteins_path]
    if proteins_path.is_dir():
        files = sorted(proteins_path.glob(prot_glob))
        return [p for p in files if p.is_file()]
    raise SystemExit(f"[err] --proteins path does not exist: {proteins_path}")


def run_hmmsearch_for_sample(
    sample: str,
    protein_faa: Path,
    hmm_files: list[Path],
    out_sample_dir: Path,
    cpus: int,
    evalue: str | None,
    domevalue: str | None,
) -> None:
    """
    Run hmmsearch for each HMM against a single protein FASTA.

    Writes one domtblout file per HMM in out_sample_dir.
    """
    out_sample_dir.mkdir(parents=True, exist_ok=True)

    for hmm in hmm_files:
        model = hmm.stem
        out_dom = out_sample_dir / f"{model}.domtblout"

        cmd = ["hmmsearch", "--cpu", str(cpus)]
        if evalue is not None:
            cmd += ["-E", str(evalue)]
        if domevalue is not None:
            cmd += ["--domE", str(domevalue)]
        cmd += ["--domtblout", str(out_dom), str(hmm), str(protein_faa)]

        try:
            sh(cmd)
        except subprocess.CalledProcessError as e:
            print(
                f"[err] hmmsearch failed for sample={sample} model={model} "
                f"with exit code {e.returncode}",
                file=sys.stderr,
            )


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Run HMMER hmmsearch for all *.hmm in a directory against one or more "
            "protein FASTA files, writing per-sample domtblout files and a summary CSV."
        )
    )
    ap.add_argument(
        "--hmm-dir",
        required=True,
        help="Directory containing *.hmm files.",
    )
    ap.add_argument(
        "--proteins",
        required=True,
        help="Protein FASTA file or directory of FASTA files.",
    )
    ap.add_argument(
        "--prot-glob",
        default="*.faa",
        help=(
            "Glob pattern for protein FASTA files if --proteins is a directory "
            "(default: '*.faa')."
        ),
    )
    ap.add_argument(
        "--out-dir",
        required=True,
        help="Output directory. Per-sample subdirectories and summary CSV will be written here.",
    )
    ap.add_argument(
        "--cpus",
        type=int,
        default=8,
        help="Number of CPU threads for hmmsearch (default: 8).",
    )
    ap.add_argument(
        "--evalue",
        type=str,
        default=None,
        help="Optional sequence E-value cutoff for reporting, for example 1e-5.",
    )
    ap.add_argument(
        "--domevalue",
        type=str,
        default=None,
        help="Optional domain E-value cutoff for reporting, for example 1e-5.",
    )
    ap.add_argument(
        "--min-bits",
        type=float,
        default=None,
        help="Optional minimum bit score filter when summarizing hits.",
    )
    ap.add_argument(
        "--max-ie",
        type=float,
        default=None,
        help="Optional maximum i-Evalue filter when summarizing hits.",
    )
    ap.add_argument(
        "--unique-per-seq",
        action="store_true",
        help="If set, count each sequence at most once per HMM (unique seqids).",
    )
    ap.add_argument(
        "--summary-name",
        default="hmmsearch_summary.csv",
        help="Name of the summary CSV written to --out-dir (default: hmmsearch_summary.csv).",
    )

    args = ap.parse_args()

    check_bin("hmmsearch")

    hmm_dir = Path(args.hmm_dir)
    proteins_path = Path(args.proteins)
    out_dir = Path(args.out_dir)

    out_dir.mkdir(parents=True, exist_ok=True)

    # Load HMM files
    hmm_files = sorted(hmm_dir.glob("*.hmm"))
    if not hmm_files:
        raise SystemExit(f"[err] no *.hmm files found in {hmm_dir}")

    print(f"[*] found {len(hmm_files)} HMMs in {hmm_dir}")

    # Find protein FASTA files
    protein_fastas = find_protein_fastas(proteins_path, args.prot_glob)
    if not protein_fastas:
        raise SystemExit(
            f"[err] no protein FASTAs found for path={proteins_path} with glob={args.prot_glob!r}"
        )

    print(f"[*] found {len(protein_fastas)} protein FASTA files to search")

    # Track per-sample information for summary
    header_maps: dict[str, dict[str, str]] = {}
    total_genes: dict[str, int] = {}
    samples: list[str] = []

    # Phase 1: run hmmsearch for each protein file
    for protein_faa in protein_fastas:
        sample = protein_faa.stem
        samples.append(sample)

        out_sample_dir = out_dir / sample
        print(f"[*] processing sample {sample} from {protein_faa}")

        run_hmmsearch_for_sample(
            sample=sample,
            protein_faa=protein_faa,
            hmm_files=hmm_files,
            out_sample_dir=out_sample_dir,
            cpus=args.cpus,
            evalue=args.evalue,
            domevalue=args.domevalue,
        )

        hdr_map = load_headers_for_fasta(protein_faa)
        header_maps[sample] = hdr_map
        total_genes[sample] = len(hdr_map)

    # Phase 2: parse domtbls and summarize counts
    counts = defaultdict(set) if args.unique_per_seq else defaultdict(int)
    headers_by_key = defaultdict(set)

    for protein_faa in protein_fastas:
        sample = protein_faa.stem
        sample_dir = out_dir / sample

        if not sample_dir.is_dir():
            print(f"[warn] no sample directory found for {sample} in {sample_dir}")
            continue

        seqid_to_header = header_maps.get(sample, {})

        for hmm in hmm_files:
            model = hmm.stem
            domtbl_path = sample_dir / f"{model}.domtblout"
            if not domtbl_path.exists():
                # no hits or run failed; skip
                continue

            for seqid, hmm_name, i_e, bits in parse_domtbl(domtbl_path):
                # basic filters
                if args.min_bits is not None and bits < args.min_bits:
                    continue
                if args.max_ie is not None and i_e > args.max_ie:
                    continue

                key = (sample, model)
                if args.unique_per_seq:
                    counts[key].add(seqid)
                else:
                    counts[key] += 1

                hdr = seqid_to_header.get(seqid, f">{seqid}")
                headers_by_key[key].add(hdr)

    # Phase 3: write summary CSV
    summary_path = out_dir / args.summary_name
    with open(summary_path, "w", newline="") as fw:
        w = csv.writer(fw)
        w.writerow(["sample", "hmm", "hits", "total_genes", "hit_headers"])

        for sample in samples:
            tg = total_genes.get(sample, 0)
            for hmm in sorted(hmm_files, key=lambda p: p.stem):
                model = hmm.stem
                key = (sample, model)
                if args.unique_per_seq:
                    val = len(counts.get(key, set()))
                else:
                    val = int(counts.get(key, 0))
                headers_list = sorted(headers_by_key.get(key, set()))
                headers_joined = "|".join(headers_list) if headers_list else ""
                w.writerow([sample, model, val, tg, headers_joined])

    print(f"[*] wrote summary: {summary_path}")


if __name__ == "__main__":
    main()
