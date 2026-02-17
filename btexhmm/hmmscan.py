#!/usr/bin/env python3
import argparse, csv, re, subprocess, sys
from pathlib import Path
from collections import defaultdict

# Output columns:
#   sample,hmm,hits,scores,ievalues,cutoff_used,total_genes,hit_headers
#
# Example:
# python /home/juneq/Toluene-HMM/btexhmm/hmmscan.py \
#   --hmm-lib /home/juneq/Toluene-HMM/btexhmm/hmms/all_models.hmm \
#   --proteins-dir /home/juneq/Toluene-HMM/btexhmm/test_genomes \
#   --out /home/juneq/Toluene_test/test_output/hmmscan_summary.csv \
#   --cpus 8


def load_cutoffs(path: Path) -> dict[str, float]:
    """
    Reads a 2-column TSV with headers: hmm, cutoff
    Returns dict[str -> float]
    """
    if not path or not path.exists():
        return {}
    mp: dict[str, float] = {}
    with open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_hmm = header.index("hmm")
            i_cut = header.index("cutoff")
        except ValueError:
            raise SystemExit("[err] cutoff file must have headers: hmm<TAB>cutoff")
        for ln, line in enumerate(fh, 2):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_hmm, i_cut):
                continue
            name = parts[i_hmm].strip()
            if not name:
                continue
            try:
                val = float(parts[i_cut])
            except ValueError:
                raise SystemExit(f"[err] non-numeric cutoff on line {ln}: {parts[i_cut]!r}")
            mp[name] = val
    return mp


def sh(cmd, quiet=False):
    if not quiet:
        print("[cmd]", " ".join(cmd), flush=True)
    # Silence hmmscan output; rely on return codes for errors
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def check_bin(name):
    from shutil import which
    if which(name) is None:
        raise SystemExit(f"[err] '{name}' not found in PATH")

def read_names_from_hmm_lib(hmm_lib: Path) -> list[str]:
    """
    Parse concatenated HMM library to extract model names (lines starting with 'NAME').
    """
    names = []
    if not hmm_lib.exists():
        return names
    with open(hmm_lib, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("NAME"):
                parts = line.strip().split()
                if len(parts) >= 2:
                    names.append(parts[1])
    return names


def read_ga_thresholds_from_hmm_lib(hmm_lib: Path) -> dict[str, tuple[float, float]]:
    """
    Parse model-specific GA thresholds from a concatenated HMM library.
    Returns: model_name -> (seq_ga_bits, dom_ga_bits)
    """
    ga_by_model: dict[str, tuple[float, float]] = {}
    if not hmm_lib.exists():
        return ga_by_model

    cur_name: str | None = None
    with open(hmm_lib, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("NAME"):
                parts = line.strip().split()
                cur_name = parts[1] if len(parts) >= 2 else None
                continue
            if cur_name and line.startswith("GA"):
                nums = re.findall(r"[-+]?\d*\.?\d+", line)
                if not nums:
                    continue
                if len(nums) == 1:
                    seq_ga = dom_ga = float(nums[0])
                else:
                    seq_ga = float(nums[0])
                    dom_ga = float(nums[1])
                ga_by_model[cur_name] = (seq_ga, dom_ga)
    return ga_by_model


def filter_domtbl_to_ga(
    in_domtbl: Path,
    out_domtbl: Path,
    ga_by_model: dict[str, tuple[float, float]],
) -> None:
    """
    Create GA-filtered domtblout from an unfiltered domtblout using model GA
    thresholds parsed from the HMM library.
    """
    warned_missing_ga: set[str] = set()
    with open(in_domtbl, "rt") as fr, open(out_domtbl, "wt") as fw:
        for line in fr:
            if line.startswith("#") or not line.strip():
                fw.write(line)
                continue
            parts = re.split(r"\s+", line.strip())
            if len(parts) < 22:
                continue
            hmm_name = parts[0]
            try:
                full_bits = float(parts[7])   # full sequence score
                dom_bits = float(parts[13])   # this domain score
            except ValueError:
                continue
            ga = ga_by_model.get(hmm_name)
            if ga is None:
                if hmm_name not in warned_missing_ga:
                    print(f"[warn] no GA threshold found in HMM library for model '{hmm_name}'; keeping row")
                    warned_missing_ga.add(hmm_name)
                fw.write(line)
                continue
            seq_ga, dom_ga = ga
            if full_bits >= seq_ga and dom_bits >= dom_ga:
                fw.write(line)


def parse_domtbl(path: Path):
    """
    Returns list of (seqid, hmm_name, i_e, bits) tuples from an HMMER hmmscan domtblout.

    For hmmscan domain table:
      target name (col 1)   -> HMM (model) name
      query name  (col 4)   -> sequence id
      i-Evalue    (col 13)  -> parts[12]
      bit score   (col 14)  -> parts[13]
    """
    rows = []
    with open(path, "rt") as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line.strip())
            if len(parts) < 22:
                continue
            hmm_name = parts[0]     # target name = HMM
            seqid = parts[3]        # query name  = protein ID
            i_e = float(parts[12])  # i-Evalue
            bits = float(parts[13]) # bit score
            rows.append((seqid, hmm_name, i_e, bits))
    return rows


def find_input_files(proteins_path: Path) -> list[Path]:
    if proteins_path.is_file():
        return [proteins_path]
    if proteins_path.is_dir():
        return [p for p in sorted(proteins_path.iterdir()) if p.is_file()]
    return []


def validate_protein_fasta(path: Path) -> tuple[bool, str]:
    """
    Validate that a file is FASTA and contains only protein-like residue symbols.
    Returns (ok, reason_if_not_ok).
    """
    allowed = set("ABCDEFGHIKLMNPQRSTUVWXYZ*-.")
    nuc_chars = set("ACGTUN")
    saw_header = False
    saw_seq = False
    total_alpha = 0
    nuc_alpha = 0
    saw_protein_specific = False
    try:
        with open(path, "rt", errors="ignore") as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    saw_header = True
                    continue
                if not saw_header:
                    return False, "missing FASTA header line starting with '>'"
                saw_seq = True
                seq = line.replace(" ", "").upper()
                bad = [c for c in seq if c not in allowed]
                if bad:
                    return False, f"contains non-protein residue character(s): {''.join(sorted(set(bad)))}"
                for c in seq:
                    if c.isalpha():
                        total_alpha += 1
                        if c in nuc_chars:
                            nuc_alpha += 1
                        else:
                            saw_protein_specific = True
    except OSError as e:
        return False, f"cannot read file: {e}"

    if not saw_header:
        return False, "no FASTA headers found"
    if not saw_seq:
        return False, "no sequence lines found"
    # Reject likely nucleotide FASTA: sequence alphabet is almost entirely A/C/G/T/U/N
    # and we never observe protein-specific residue letters.
    if total_alpha > 0:
        nuc_frac = nuc_alpha / total_alpha
        if nuc_frac >= 0.95 and not saw_protein_specific:
            return False, "sequence composition looks nucleotide (A/C/G/T/U/N only)"
    return True, ""


def load_headers_for_fasta(protein_faa: Path) -> dict[str, str]:
    """
    Build mapping from sequence id -> full header line (including leading '>') for a single FASTA.
    The id is the first whitespace-delimited token after '>'.
    """
    mp: dict[str, str] = {}
    try:
        with open(protein_faa, "rt", errors="ignore") as fh:
            for ln in fh:
                if not ln.startswith(">"):
                    continue
                hdr = ln.rstrip("\n\r")
                seqid = hdr[1:].split()[0]
                mp.setdefault(seqid, hdr)  # keep first occurrence
    except FileNotFoundError:
        pass
    return mp


def run_hmmscan_for_sample(
    gdir: Path,
    protein_faa: Path,
    hmm_lib: Path,
    domtbl_dir: Path,
    cpus: int,
    evalue: str | None,
    domevalue: str | None,
    ga_by_model: dict[str, tuple[float, float]],
    skip_existing: bool,
):
    """
    Run hmmscan once per sample against the concatenated HMM library.
    Writes per-sample domtbl outputs:
      - <sample>.pre_ga.domtblout  (no --cut_ga)
      - <sample>.ga.domtblout      (filtered from pre-GA using model GA thresholds)
    """
    if not protein_faa or not protein_faa.exists():
        print(f"[warn] protein FASTA missing for {protein_faa} — skipping hmmscan")
        return

    if not hmm_lib.exists():
        raise SystemExit(f"[err] HMM library not found: {hmm_lib}")

    domtbl_dir.mkdir(parents=True, exist_ok=True)
    sample = protein_faa.stem
    out_dom_pre = domtbl_dir / f"{sample}.pre_ga.domtblout"
    out_dom_ga = domtbl_dir / f"{sample}.ga.domtblout"
    if skip_existing:
        if (
            out_dom_pre.exists() and out_dom_pre.stat().st_size > 0
            and out_dom_ga.exists() and out_dom_ga.stat().st_size > 0
        ):
            return

    # 1) Always write unfiltered domtbl (before GA).
    print(f"[info] running hmmscan on {protein_faa.name}")
    cmd_pre = ["hmmscan", "--cpu", str(cpus)]
    if evalue:
        cmd_pre += ["-E", str(evalue)]
    if domevalue:
        cmd_pre += ["--domE", str(domevalue)]
    cmd_pre += ["--domtblout", str(out_dom_pre), str(hmm_lib), str(protein_faa)]
    try:
        sh(cmd_pre, quiet=True)
    except subprocess.CalledProcessError as e:
        print(
            f"[err] hmmscan (pre-GA) failed for {protein_faa.name} with code {e.returncode}",
            file=sys.stderr,
        )
        return

    # 2) Derive GA-filtered domtbl from pre-GA output (no second hmmscan run).
    filter_domtbl_to_ga(out_dom_pre, out_dom_ga, ga_by_model)


def parse_models_arg(models_arg: str | None) -> set[str] | None:
    """
    Accepts comma/space separated names; returns a set or None if not provided.
    """
    if not models_arg:
        return None
    toks = re.split(r"[,\s]+", models_arg.strip())
    sel = {t for t in (x.strip() for x in toks) if t}
    return sel or None


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(
        description=(
            "Run hmmscan vs protein FASTAs using a concatenated HMM library, "
            "writing sample,hmm,hits,scores,ievalues,cutoff_used,total_genes,hit_headers. "
            "Sample name is taken from the FASTA basename."
        )
    )
    ap.add_argument("--hmm-dir", required=False, help="directory with individual *.hmm (for model name list; optional if --hmm-lib contains names)")
    ap.add_argument("--hmm-lib", required=True, help="concatenated HMM library for hmmscan")
    ap.add_argument(
        "--proteins-dir",
        "--genomes-dir",
        dest="proteins_dir",
        required=True,
        help=(
            "Directory containing protein FASTA files, or a single protein FASTA file. "
            "Supports legacy alias --genomes-dir."
        ),
    )
    ap.add_argument(
        "--domtbl-subdir",
        default="hmmscan_output",
        help="subfolder under the --out directory to store per-sample domtbl outputs (e.g., hmmscan_output/<sample>/<sample>.domtblout)",
    )
    ap.add_argument("--cpus", type=int, default=8, help="threads for hmmscan")
    ap.add_argument("--evalue", type=str, default=None, help="sequence E-value cutoff for reporting, e.g. 1e-5")
    ap.add_argument(
        "--domevalue",
        type=str,
        default=None,
        help="domain E-value cutoff for reporting and for filtering parsed hits, e.g. 1e-5",
    )
    ap.add_argument("--min-bits", type=float, default=None, help="optional min bit score when counting")
    ap.add_argument("--max-ie", type=float, default=None, help="optional max i-Evalue when counting")
    ap.add_argument("--unique-per-seq", action="store_true", help="count at most one hit per sequence per HMM")
    ap.add_argument(
        "--skip-existing",
        action="store_true",
        help="skip running hmmscan if sample domtbl already exists",
    )
    ap.add_argument(
        "--all-hits-per-protein",
        action="store_true",
        help=(
            "Count all HMM hits for a protein that pass filters. "
            "Default is to count only the single best-scoring HMM hit per protein."
        ),
    )
    ap.add_argument(
        "--out",
        required=True,
        help="CSV path for output counts (sample,hmm,hits,scores,ievalues,cutoff_used,total_genes,hit_headers)",
    )
    ap.add_argument("--cutoffs", default=None, help="TSV with columns: hmm<TAB>cutoff (may include many models)")
    ap.add_argument(
        "--models",
        default=None,
        help=(
            "Optional comma/space separated list of HMM names (stems) to count, e.g. 'TmoA,TmoB'. "
            "If omitted, reports all models found in --hmm-dir."
        ),
    )
    ap.add_argument(
        "--fallback-domtbl-root",
        default=None,
        help=(
            "Optional root containing per-sample domtbl folders: "
            "<root>/<sample>/domtbl/*.domtblout or *.domtbl (used if local domtbl missing)"
        ),
    )
    args = ap.parse_args(argv)

    check_bin("hmmscan")

    hmm_lib = Path(args.hmm_lib)
    proteins_path = Path(args.proteins_dir)
    ga_by_model = read_ga_thresholds_from_hmm_lib(hmm_lib)

    source_desc = ""
    out_path = Path(args.out)
    output_root = out_path.parent if out_path.parent != Path("") else Path(".")

    # Load model names either from hmm-dir (preferred) or from the concatenated hmm-lib
    if args.hmm_dir:
        hmm_dir = Path(args.hmm_dir)
        all_hmm_files = sorted(hmm_dir.glob("*.hmm"))
        if not all_hmm_files:
            raise SystemExit(f"[err] no *.hmm in {hmm_dir}")
        all_model_names = [h.stem for h in all_hmm_files]
        source_desc = f"hmm_dir={hmm_dir}"
    else:
        if not hmm_lib.exists():
            raise SystemExit(f"[err] HMM library not found: {hmm_lib}")
        all_model_names = read_names_from_hmm_lib(hmm_lib)
        if not all_model_names:
            raise SystemExit(f"[err] no model names found in HMM library: {hmm_lib}")
        source_desc = f"hmm_lib={hmm_lib}"
    all_model_set = set(all_model_names)

    # Resolve model selection (subset) if provided
    selected = parse_models_arg(args.models)
    if selected is not None:
        unknown = sorted(selected - all_model_set)
        if unknown:
            print(f"[warn] --models includes names not found in {source_desc}: {', '.join(unknown)}")
        model_set = selected & all_model_set
        if not model_set:
            raise SystemExit("[err] none of the requested --models were found in the provided HMM source")
    else:
        model_set = all_model_set

    print(f"[*] tracking {len(model_set)} HMM(s) from {source_desc}: {', '.join(sorted(model_set))}")

    protein_fastas = find_input_files(proteins_path)
    if not protein_fastas:
        raise SystemExit(
            f"[err] no input files found in {proteins_path}"
        )
    for protein_faa in protein_fastas:
        ok, reason = validate_protein_fasta(protein_faa)
        if not ok:
            print(f"[err] {protein_faa.name}")
            print(f"[err] {protein_faa.name} is not a proper protein sequence file: {reason}")
            raise SystemExit(1)

    # Load full cutoffs table (may include many models); use only for selected models
    cutmap_all = load_cutoffs(Path(args.cutoffs)) if args.cutoffs else {}
    cutmap = {m: cutmap_all[m] for m in model_set if m in cutmap_all}

    # Apply an i-Evalue threshold when counting: prefer explicit --max-ie, else reuse --domevalue if numeric
    max_ie = args.max_ie
    if max_ie is None and args.domevalue is not None:
        try:
            max_ie = float(args.domevalue)
        except ValueError:
            print(f"[warn] unable to parse --domevalue={args.domevalue!r} as float; ignoring for filtering")

    header_maps: dict[str, dict[str, str]] = {}
    total_genes: dict[str, int] = {}
    sample_parent: dict[str, Path] = {}
    samples: list[str] = []

    # Phase 1: ensure per-sample domtbl exists (one hmmscan per sample)
    for protein_faa in protein_fastas:
        sample = protein_faa.stem
        samples.append(sample)
        sample_parent[sample] = protein_faa.parent

        domtbl_dir = output_root / args.domtbl_subdir / sample
        run_hmmscan_for_sample(
            gdir=protein_faa.parent,
            protein_faa=protein_faa,
            hmm_lib=hmm_lib,
            domtbl_dir=domtbl_dir,
            cpus=args.cpus,
            evalue=args.evalue,
            domevalue=args.domevalue,
            ga_by_model=ga_by_model,
            skip_existing=args.skip_existing,
        )

        mp = load_headers_for_fasta(protein_faa)
        header_maps[sample] = mp
        total_genes[sample] = len(mp)

    # Phase 2: parse domtbls and build counts/scores
    counts = defaultdict(set) if args.unique_per_seq else defaultdict(int)
    headers_by_key = defaultdict(set)
    scores_by_key = defaultdict(list)
    evalues_by_key = defaultdict(list)
    thresh_by_key = {}
    warned_missing_cut = set()

    for sample in samples:
        gdir = sample_parent.get(sample)
        if gdir is None:
            continue

        domtbl_dir = output_root / args.domtbl_subdir / sample
        domtbl_path = domtbl_dir / f"{sample}.ga.domtblout"

        # Fallback to external root if local domtbl missing
        if not domtbl_path.exists() and args.fallback_domtbl_root:
            fbdir = Path(args.fallback_domtbl_root) / sample / "domtbl"
            if fbdir.exists():
                fb_candidates = (
                    sorted(fbdir.glob("*.ga.domtblout"))
                    or sorted(fbdir.glob("*.domtblout"))
                    or sorted(fbdir.glob("*.domtbl"))
                )
                if fb_candidates:
                    domtbl_path = fb_candidates[0]
                    print(f"[info] using fallback domtbl for {sample} from {domtbl_path}")
        if not domtbl_path.exists():
            print(f"[info] no domtbl for {sample}; skipping")
            continue

        seqid_to_header = header_maps.get(sample, {})

        if not args.all_hits_per_protein:
            hits_per_protein = defaultdict(list)
            for seqid, model_name, i_e, bits in parse_domtbl(domtbl_path):
                model = model_name
                if model not in model_set:
                    continue

                cutoff = cutmap_all.get(model)
                if cutoff is not None:
                    if bits < cutoff:
                        continue
                elif model not in warned_missing_cut and args.cutoffs:
                    print(
                        f"[warn] no cutoff found for HMM '{model}' in {args.cutoffs}; "
                        "using only other filters for this model"
                    )
                    warned_missing_cut.add(model)

                if args.min_bits is not None and bits < args.min_bits:
                    continue
                if max_ie is not None and i_e > max_ie:
                    continue

                diff_score = bits - (cutoff if cutoff is not None else 0.0)
                hits_per_protein[seqid].append(
                    {
                        "model": model,
                        "seqid": seqid,
                        "diff": diff_score,
                        "bits": bits,
                        "i_e": i_e,
                        "cutoff": cutoff,
                    }
                )

            for seqid, hit_list in hits_per_protein.items():
                if not hit_list:
                    continue
                best_hit = max(hit_list, key=lambda h: h["diff"])
                key = (sample, best_hit["model"])
                if args.unique_per_seq:
                    before = len(counts[key])
                    counts[key].add(seqid)
                    if len(counts[key]) > before:
                        scores_by_key[key].append(best_hit["bits"])
                        evalues_by_key[key].append(best_hit["i_e"])
                else:
                    counts[key] += 1
                    scores_by_key[key].append(best_hit["bits"])
                    evalues_by_key[key].append(best_hit["i_e"])
                hdr = seqid_to_header.get(seqid, f">{seqid}")
                headers_by_key[key].add(hdr)
                if best_hit["cutoff"] is not None and key not in thresh_by_key:
                    thresh_by_key[key] = best_hit["cutoff"]

        else:
            for seqid, model_name, i_e, bits in parse_domtbl(domtbl_path):
                model = model_name
                if model not in model_set:
                    continue

                cutoff = cutmap_all.get(model)
                if cutoff is not None:
                    if bits < cutoff:
                        continue
                elif model not in warned_missing_cut and args.cutoffs:
                    print(
                        f"[warn] no cutoff found for HMM '{model}' in {args.cutoffs}; "
                        "using only other filters for this model"
                    )
                    warned_missing_cut.add(model)

                if args.min_bits is not None and bits < args.min_bits:
                    continue
                if max_ie is not None and i_e > max_ie:
                    continue

                key = (sample, model)
                if args.unique_per_seq:
                    before = len(counts[key])
                    counts[key].add(seqid)
                    if len(counts[key]) > before:
                        scores_by_key[key].append(bits)
                        evalues_by_key[key].append(i_e)
                else:
                    counts[key] += 1
                    scores_by_key[key].append(bits)
                    evalues_by_key[key].append(i_e)
                hdr = seqid_to_header.get(seqid, f">{seqid}")
                headers_by_key[key].add(hdr)
                if cutoff is not None and key not in thresh_by_key:
                    thresh_by_key[key] = cutoff

    # Phase 3: write output
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fw:
        w = csv.writer(fw)
        w.writerow(
            [
                "sample",
                "hmm",
                "hits",
                "scores",
                "ievalues",
                "cutoff_used",
                "total_genes",
                "hit_headers",
            ]
        )
        for sample in samples:
            tg = total_genes.get(sample, 0)
            for hmm in sorted(model_set):
                key = (sample, hmm)
                if args.unique_per_seq:
                    val = len(counts.get(key, set()))
                else:
                    val = int(counts.get(key, 0))
                scores_list = scores_by_key.get(key, [])
                scores_joined = "|".join(str(s) for s in scores_list) if scores_list else ""
                evalues_list = evalues_by_key.get(key, [])
                evalues_joined = "|".join(str(ev) for ev in evalues_list) if evalues_list else ""
                thresh_val = thresh_by_key.get(key, "")
                headers_list = sorted(headers_by_key.get(key, set()))
                headers_joined = "|".join(headers_list) if headers_list else ""
                w.writerow([sample, hmm, val, scores_joined, evalues_joined, thresh_val, tg, headers_joined])

    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
