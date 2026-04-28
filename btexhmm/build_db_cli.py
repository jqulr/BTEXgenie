#!/usr/bin/env python3
import argparse
import gzip
import os
import shutil
import tarfile
import urllib.error
import urllib.request
from pathlib import Path


KOFAM_KO_LIST_URL = "https://www.genome.jp/ftp/db/kofam/ko_list.gz"
KOFAM_PROFILES_URL = "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz"
APPROX_DOWNLOAD_SIZES = {
    KOFAM_KO_LIST_URL: "approximately 884K",
    KOFAM_PROFILES_URL: "approximately 1.4G",
}


def check_bin(name: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"[err] '{name}' not found in PATH")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Download and prepare the KOfam database for BTEX-HMM "
            "with a Conda activation hook that exports KOFAM_DB."
        )
    )
    parser.add_argument(
        "--db-dir",
        required=True,
        help="Directory for the KOfam database.",
    )
    parser.add_argument(
        "--conda-prefix",
        default=os.environ.get("CONDA_PREFIX"),
        help=(
            "Conda environment prefix where activate.d is written. "
            "Uses default CONDA_PREFIX when available."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Redownload database and overwrite extracted database files.",
    )
    return parser.parse_args(argv)


def download_file(url: str, dest: Path, force: bool) -> None:
    if dest.exists() and not force:
        print(f"[info] reusing existing download: {dest}")
        return

    dest.parent.mkdir(parents=True, exist_ok=True)
    approx_size = APPROX_DOWNLOAD_SIZES.get(url)
    if approx_size:
        print(f"[*] downloading {url} ({approx_size})")
    else:
        print(f"[*] downloading {url}")
    try:
        with urllib.request.urlopen(url) as response, open(dest, "wb") as fh:
            shutil.copyfileobj(response, fh)
    except urllib.error.URLError as exc:
        raise SystemExit(f"[err] failed to download {url}: {exc}") from exc
    print(f"[info] wrote {dest}")


def extract_ko_list(src_gz: Path, out_txt: Path, force: bool) -> None:
    if out_txt.exists() and not force:
        print(f"[info] reusing existing file: {out_txt}")
        return

    print(f"[info] decompressing {src_gz.name}")
    with gzip.open(src_gz, "rb") as src, open(out_txt, "wb") as dst:
        shutil.copyfileobj(src, dst)
    print(f"[info] wrote {out_txt}")


def extract_profiles(src_tar_gz: Path, db_dir: Path, force: bool) -> None:
    profiles_dir = db_dir / "profiles"
    if profiles_dir.is_dir() and any(profiles_dir.iterdir()) and not force:
        print(f"[info] reusing existing directory: {profiles_dir}")
        return

    if profiles_dir.exists() and force:
        shutil.rmtree(profiles_dir)

    print(f"[info] decompressing {src_tar_gz.name}")
    with tarfile.open(src_tar_gz, "r:gz") as tar:
        tar.extractall(path=db_dir)
    print(f"[info] extracted profiles to {profiles_dir}")


def validate_db(db_dir: Path) -> None:
    profiles_dir = db_dir / "profiles"
    ko_list = db_dir / "ko_list"
    if not profiles_dir.is_dir():
        raise SystemExit(f"[err] expected extracted profiles dir at {profiles_dir}")
    if not any(profiles_dir.iterdir()):
        raise SystemExit(f"[err] profiles dir is empty: {profiles_dir}")
    if not ko_list.is_file():
        raise SystemExit(f"[err] expected extracted ko list at {ko_list}")
    print(f"[info] validated database at {db_dir}")


def write_activate_hook(db_dir: Path, conda_prefix: Path) -> Path:
    hook_dir = conda_prefix / "etc" / "conda" / "activate.d"
    hook_dir.mkdir(parents=True, exist_ok=True)
    hook_path = hook_dir / "kofam.sh"
    hook_path.write_text(
        "# Added by btex-build-db\n"
        f"export KOFAM_DB={db_dir}\n",
        encoding="utf-8",
    )
    print(f"[info] wrote activate hook: {hook_path}")
    return hook_path


def main(argv=None):
    args = parse_args(argv)
    check_bin("exec_annotation")

    if not args.conda_prefix:
        raise SystemExit(
            "[err] no Conda prefix found. Activate the target Conda environment "
            "first or rerun with --conda-prefix."
        )

    db_dir = Path(args.db_dir).resolve()
    db_dir.mkdir(parents=True, exist_ok=True)

    ko_list_gz = db_dir / "ko_list.gz"
    profiles_tar_gz = db_dir / "profiles.tar.gz"
    ko_list_txt = db_dir / "ko_list"

    download_file(KOFAM_KO_LIST_URL, ko_list_gz, args.force)
    download_file(KOFAM_PROFILES_URL, profiles_tar_gz, args.force)
    extract_ko_list(ko_list_gz, ko_list_txt, args.force)
    extract_profiles(profiles_tar_gz, db_dir, args.force)
    validate_db(db_dir)
    write_activate_hook(db_dir, Path(args.conda_prefix).resolve())

    print("[info] database setup complete")
    print(f"[info] KOFAM_DB should be set to: {db_dir}")


if __name__ == "__main__":
    main()
