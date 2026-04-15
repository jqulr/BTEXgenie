#!/usr/bin/env python3
import argparse
import csv
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


DEFAULT_CPUS = 8


def check_bin(name: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"[err] '{name}' not found in PATH")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Run KOfamScan exec_annotation on protein FASTA inputs and write only "
            "the above-threshold KO hits for each genome."
        )
    )
    parser.add_argument(
        "--genomes-dir",
        required=True,
        help=(
            "Directory of genome subfolders (each containing one .faa), or a "
            "directory of .faa files."
        ),
    )
    parser.add_argument(
        "--db-dir",
        default=os.environ.get("KOFAM_DB"),
        help=(
            "KOfam DB path. Defaults to the KOFAM_DB environment variable. "
            "The directory must contain profiles/ and ko_list or ko_list.gz."
        ),
    )
    parser.add_argument(
        "--cpus",
        type=int,
        default=DEFAULT_CPUS,
        help=f"Threads for exec_annotation (default: {DEFAULT_CPUS}).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun KOfamScan even if kofam_abv_thres.tsv already exists.",
    )
    return parser.parse_args(argv)


def resolve_ko_list(db_dir: Path) -> Path:
    ko_list = db_dir / "ko_list"
    ko_list_gz = db_dir / "ko_list.gz"
    if ko_list.is_file():
        return ko_list
    if ko_list_gz.is_file():
        return ko_list_gz
    raise SystemExit(
        f"[err] Neither '{ko_list}' nor '{ko_list_gz}' was found."
    )


def validate_kofam_db(db_dir: Path) -> Path:
    profiles_dir = db_dir / "profiles"
    if not profiles_dir.is_dir():
        raise SystemExit(
            f"[err] Invalid KOFAM_DB: expected '{profiles_dir}' to exist."
        )

    ko_list = resolve_ko_list(db_dir)
    print(f"[info] KOfam database detected: {db_dir}")
    print(f"[info] profiles dir: {profiles_dir}")
    print(f"[info] ko list     : {ko_list}")
    return ko_list


def iter_input_fastas(genomes_dir: Path):
    processed = False

    for genome_subdir in sorted(p for p in genomes_dir.iterdir() if p.is_dir()):
        faa_files = sorted(genome_subdir.glob("*.faa"))
        if not faa_files:
            continue
        processed = True
        yield genome_subdir.name, faa_files[0], genome_subdir

    for protein_file in sorted(genomes_dir.glob("*.faa")):
        processed = True
        out_dir = genomes_dir / protein_file.stem
        yield protein_file.stem, protein_file, out_dir

    if not processed:
        raise SystemExit(
            "[err] No input .faa files found. Expected either genome subdirectories "
            f"under '{genomes_dir}' or .faa files directly inside it."
        )


def run_exec_annotation(protein_file: Path, db_dir: Path, ko_list: Path, cpus: int, detail_out: Path):
    with tempfile.TemporaryDirectory(prefix=f"kofam_{protein_file.stem}_") as tmp_dir:
        cmd = [
            "exec_annotation",
            "-o",
            str(detail_out),
            "-f",
            "detail-tsv",
            "--cpu",
            str(cpus),
            "--tmp-dir",
            tmp_dir,
            "-p",
            str(db_dir / "profiles") + "/",
            "-k",
            str(ko_list),
            str(protein_file),
        ]
        print("[cmd]", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)


def write_above_threshold(detail_tsv: Path, out_tsv: Path) -> int:
    rows_written = 0
    with open(detail_tsv, "rt", newline="") as src, open(out_tsv, "wt", newline="") as dst:
        writer = csv.writer(dst, delimiter="\t")
        writer.writerow(["gene_name", "KO", "threshold", "score", "e_value"])
        for raw in src:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 6 and fields[0] == "*":
                writer.writerow(fields[1:6])
                rows_written += 1
    return rows_written


def main(argv=None):
    args = parse_args(argv)
    check_bin("exec_annotation")

    genomes_dir = Path(args.genomes_dir).resolve()
    if not genomes_dir.is_dir():
        raise SystemExit(f"[err] genomes dir not found: {genomes_dir}")

    if not args.db_dir:
        raise SystemExit(
            "[err] No KOfam database provided. Set KOFAM_DB or pass --db-dir."
        )

    db_dir = Path(args.db_dir).resolve()
    ko_list = validate_kofam_db(db_dir)

    print("Starting KOfamScan above-threshold annotation")
    print("-------------------------------------------")
    print(f"DB dir      : {db_dir}")
    print(f"Genomes dir : {genomes_dir}")
    print(f"CPUs        : {args.cpus}")
    print("-------------------------------------------")

    for genome_name, protein_file, out_dir in iter_input_fastas(genomes_dir):
        out_dir.mkdir(parents=True, exist_ok=True)
        final_out = out_dir / "kofam_abv_thres.tsv"
        if final_out.exists() and final_out.stat().st_size > 0 and not args.force:
            print(f"[{genome_name}] Reusing existing {final_out.name}")
            print("-------------------------------------------")
            continue

        print(f"[{genome_name}] Running exec_annotation")
        with tempfile.NamedTemporaryFile(
            mode="wt",
            suffix=".tsv",
            prefix=f"{genome_name}_kofam_detail_",
            delete=False,
        ) as handle:
            detail_tmp = Path(handle.name)

        try:
            run_exec_annotation(protein_file, db_dir, ko_list, args.cpus, detail_tmp)
            count = write_above_threshold(detail_tmp, final_out)
        finally:
            detail_tmp.unlink(missing_ok=True)

        print(f"[{genome_name}] Wrote {count} hits -> {final_out}")
        print("-------------------------------------------")

    print("All genomes processed.")


if __name__ == "__main__":
    main()
