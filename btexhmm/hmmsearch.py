#!/usr/bin/env python3
import argparse, csv, re, subprocess, sys
from pathlib import Path
from collections import defaultdict

# Output columns:
#   sample,hmm,hits,score,thres_score,total_genes,hit_headers

def collect_domtbls_for_sample(gdir: Path, sample: str, models_set: set[str],
                               domtbl_subdir: str, fallback_root: str | None):
    """
    Returns domtbl Path objects to parse for this sample, preferring:
      <gdir>/<domtbl_subdir>/<sample>/*.domtblout (excluding pfam)
    If none found and fallback_root is provided, optionally fall back to:
      <fallback_root>/<sample>/domtbl/*.domtbl
    Only files whose stem is in models_set are returned.
    """
    primary = gdir / domtbl_subdir / sample
    primary_files = sorted(primary.glob("*.domtblout")) if primary.exists() else []
    primary_keep = [p for p in primary_files if p.stem in models_set and p.stem != "pfam"]
    if primary_keep:
        return primary_keep

    if fallback_root:
        fbdir = Path(fallback_root) / sample / "domtbl"
        if fbdir.exists():
            fb_files = sorted(fbdir.glob("*.domtbl"))
            fb_keep = [p for p in fb_files if p.stem in models_set]
            if fb_keep:
                print(f"[info] using fallback domtbls for {sample} from {fbdir}")
                return fb_keep
    return []

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
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def check_bin(name):
    from shutil import which
    if which(name) is None:
        raise SystemExit(f"[err] '{name}' not found in PATH")

def parse_domtbl(path: Path):
    """
    Returns list of (seqid, hmm_name, i_e, bits) tuples from an HMMER domtblout.
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
            hmm_name = parts[3]     # query name (model ID)
            i_e = float(parts[12])  # i-Evalue
            bits = float(parts[13]) # bit score
            rows.append((seqid, hmm_name if hmm_name != "-" else path.stem, i_e, bits))
    return rows

def find_protein_fastas(gdir: Path, prot_glob: str) -> list[Path]:
    files = sorted(gdir.glob(prot_glob))
    return [p for p in files if p.is_file()]

def run_prodigal_on_fastas(
    nuc_fastas: list[Path],
    out_dir: Path,
    mode: str = "single",
    skip_existing: bool = False,
) -> list[Path]:
    """
    Run Prodigal on nucleotide FASTA inputs, emitting protein FASTAs (.faa) into out_dir.
    Returns list of output paths.
    """
    if not nuc_fastas:
        return []
    check_bin("prodigal")
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for nf in nuc_fastas:
        out_prot = out_dir / f"{nf.stem}.faa"
        if skip_existing and out_prot.exists() and out_prot.stat().st_size > 0:
            print(f"[info] skipping prodigal for {nf.name}; found existing {out_prot.name}")
            outputs.append(out_prot)
            continue
        cmd = ["prodigal", "-i", str(nf), "-a", str(out_prot), "-p", mode, "-q"]
        try:
            sh(cmd, quiet=False)
            outputs.append(out_prot)
        except subprocess.CalledProcessError as e:
            print(f"[err] prodigal failed for {nf} with code {e.returncode}", file=sys.stderr)
    return outputs

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

def annotate_proteins(
    proteins_path: Path,
    out_csv: Path,
    hmm_dir: Path,
    cutoffs_path: Path | None = None,
    cpus: int = 8,
    evalue: str | None = None,
    domevalue: str | None = None,
    min_bits: float | None = None,
    max_ie: float | None = None,
    unique_per_seq: bool = False,
    all_hits_per_protein: bool = False,
    skip_existing: bool = False,
    domtbl_subdir: str = "domtbl",
    fallback_domtbl_root: str | None = None,
    prot_glob: tuple[str, ...] | str = ("*.faa", "*.fa", "*.fasta", "*.fna"),
    run_prodigal: bool = False,
    prodigal_out_dir: Path | None = None,
    prodigal_mode: str = "single",
    write_output: bool = True,
) -> list[list]:
    """
    Run hmmsearch for all protein FASTA files in a directory (or a single FASTA)
    against all HMMs in hmm_dir and write a combined summary CSV. If no protein
    FASTAs are found but nucleotide FASTAs are present (or run_prodigal is set),
    Prodigal is run first to generate .faa files.
    Returns list of output rows (including sample name).
    """
    check_bin("hmmsearch")

    proteins_path = Path(proteins_path)
    out_csv = Path(out_csv)
    hmm_dir = Path(hmm_dir)

    all_hmm_files = sorted(hmm_dir.glob("*.hmm"))
    if not all_hmm_files:
        raise SystemExit(f"[err] no *.hmm in {hmm_dir}")
    model_set = {h.stem for h in all_hmm_files}
    cutmap_all = load_cutoffs(Path(cutoffs_path)) if cutoffs_path else {}

    patterns = (prot_glob,) if isinstance(prot_glob, str) else tuple(prot_glob)
    protein_files: list[Path] = []
    nucleotide_files: list[Path] = []
    if proteins_path.is_dir():
        for pat in patterns:
            protein_files.extend(sorted(proteins_path.glob(pat)))
    elif proteins_path.is_file():
        protein_files = [proteins_path]
    else:
        raise SystemExit(f"[err] proteins_path not found: {proteins_path}")

    # Reclassify matches into protein vs nucleotide lists so we can optionally run Prodigal.
    classified_protein_files: list[Path] = []
    for p in protein_files:
        if not p.is_file():
            continue
        suffix = p.suffix.lower()
        if suffix in {".faa", ".aa", ".pep"}:
            classified_protein_files.append(p)
        elif suffix in {".fna", ".fa", ".fasta", ".ffn"}:
            nucleotide_files.append(p)
        else:
            classified_protein_files.append(p)
    protein_files = classified_protein_files

    auto_run_prodigal = not protein_files and nucleotide_files
    if run_prodigal or auto_run_prodigal:
        prod_out_dir = (
            Path(prodigal_out_dir)
            if prodigal_out_dir
            else (proteins_path.parent / "prodigal_proteins" if proteins_path.is_file() else proteins_path / "prodigal_proteins")
        )
        print(f"[info] running prodigal on {len(nucleotide_files)} nucleotide FASTA(s) -> {prod_out_dir}")
        prot_from_prodigal = run_prodigal_on_fastas(
            nuc_fastas=nucleotide_files,
            out_dir=prod_out_dir,
            mode=prodigal_mode,
            skip_existing=skip_existing,
        )
        protein_files.extend(prot_from_prodigal)

    if not protein_files:
        raise SystemExit(f"[err] no FASTA files matching {patterns} in {proteins_path}")

    def _annotate_single_fasta(protein_faa: Path) -> list[list]:
        sample = protein_faa.stem
        domtbl_dir = protein_faa.parent / domtbl_subdir / sample
        run_hmmsearch_for_sample(
            gdir=protein_faa.parent,
            protein_faa=protein_faa,
            hmm_files=all_hmm_files,
            domtbl_dir=domtbl_dir,
            cpus=cpus,
            evalue=evalue,
            domevalue=domevalue,
            skip_existing=skip_existing,
        )

        header_map = load_headers_for_fasta(protein_faa)
        total_genes = len(header_map)

        counts = defaultdict(set) if unique_per_seq else defaultdict(int)
        headers_by_key = defaultdict(set)
        scores_by_key = defaultdict(list)
        thresh_by_key = {}
        warned_missing_cut = set()

        domtbl_paths = collect_domtbls_for_sample(
            protein_faa.parent, sample, model_set, domtbl_subdir, fallback_domtbl_root
        )
        if not domtbl_paths:
            print(f"[info] no per-model domtbls for {sample}; skipping")
        else:
            if not all_hits_per_protein:
                # Count only the best-scoring HMM hit per protein
                hits_per_protein = defaultdict(list)
                for dom in domtbl_paths:
                    dom_model_stem = dom.stem
                    for seqid, model_name, i_e, bits in parse_domtbl(dom):
                        model = model_name if model_name != "-" else dom_model_stem
                        if model not in model_set:
                            continue
                        cutoff = cutmap_all.get(model)
                        if cutoff is not None:
                            if bits < cutoff:
                                continue
                        elif model not in warned_missing_cut:
                            print(f"[warn] no cutoff found for HMM '{model}' in {cutoffs_path}; using only other filters for this model")
                            warned_missing_cut.add(model)
                        if min_bits is not None and bits < min_bits:
                            continue
                        if max_ie is not None and i_e > max_ie:
                            continue
                        diff_score = bits - (cutoff if cutoff is not None else 0.0)
                        hits_per_protein[seqid].append({
                            "model": model,
                            "seqid": seqid,
                            "diff": diff_score,
                            "bits": bits,
                            "cutoff": cutoff
                        })
                for seqid, hit_list in hits_per_protein.items():
                    if not hit_list:
                        continue
                    best_hit = max(hit_list, key=lambda h: h['diff'])
                    key = (sample, best_hit['model'])
                    if unique_per_seq:
                        before = len(counts[key])
                        counts[key].add(seqid)
                        if len(counts[key]) > before:
                            scores_by_key[key].append(best_hit["bits"])
                    else:
                        counts[key] += 1
                        scores_by_key[key].append(best_hit["bits"])
                    hdr = header_map.get(seqid, f">{seqid}")
                    headers_by_key[key].add(hdr)
                    if best_hit["cutoff"] is not None and key not in thresh_by_key:
                        thresh_by_key[key] = best_hit["cutoff"]
            else:
                # Count all HMM hits per protein that pass filters
                for dom in domtbl_paths:
                    dom_model_stem = dom.stem
                    for seqid, model_name, i_e, bits in parse_domtbl(dom):
                        model = model_name if model_name != "-" else dom_model_stem
                        if model not in model_set:
                            continue
                        cutoff = cutmap_all.get(model)
                        if cutoff is not None:
                            if bits < cutoff:
                                continue
                        elif model not in warned_missing_cut:
                            print(f"[warn] no cutoff found for HMM '{model}' in {cutoffs_path}; using only other filters for this model")
                            warned_missing_cut.add(model)
                        if min_bits is not None and bits < min_bits:
                            continue
                        if max_ie is not None and i_e > max_ie:
                            continue
                        key = (sample, model)
                        if unique_per_seq:
                            before = len(counts[key])
                            counts[key].add(seqid)
                            if len(counts[key]) > before:
                                scores_by_key[key].append(bits)
                        else:
                            counts[key] += 1
                            scores_by_key[key].append(bits)
                        hdr = header_map.get(seqid, f">{seqid}")
                        headers_by_key[key].add(hdr)
                        if cutoff is not None and key not in thresh_by_key:
                            thresh_by_key[key] = cutoff

        rows: list[list] = []
        for hmm in sorted(model_set):
            key = (sample, hmm)
            if unique_per_seq:
                val = len(counts.get(key, set()))
            else:
                val = int(counts.get(key, 0))
            scores_list = scores_by_key.get(key, [])
            scores_joined = "|".join(str(s) for s in scores_list) if scores_list else ""
            thresh_val = thresh_by_key.get(key, "")
            headers_list = sorted(headers_by_key.get(key, set()))
            headers_joined = "|".join(headers_list) if headers_list else ""
            row = [sample, hmm, val, scores_joined, thresh_val, total_genes, headers_joined]
            rows.append(row)
        return rows

    all_rows: list[list] = []
    for pf in protein_files:
        all_rows.extend(_annotate_single_fasta(pf))

    if all_rows and write_output:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        with open(out_csv, "w", newline="") as fw:
            w = csv.writer(fw)
            w.writerow(["sample","hmm","hits","score","thres_score","total_genes","hit_headers"])
            w.writerows(all_rows)

    return all_rows

def run_hmmsearch_for_sample(gdir: Path, protein_faa: Path, hmm_files: list[Path],
                             domtbl_dir: Path, cpus: int,
                             evalue: str | None, domevalue: str | None,
                             skip_existing: bool):
    """
    For each selected HMM, run hmmsearch against this single protein FASTA.
    Writes <domtbl_dir>/<model>.domtblout
    """
    if not protein_faa or not protein_faa.exists():
        print(f"[warn] protein FASTA missing for {protein_faa} — skipping hmmsearch")
        return

    domtbl_dir.mkdir(parents=True, exist_ok=True)

    for hmm in hmm_files:
        model = hmm.stem
        out_dom = domtbl_dir / f"{model}.domtblout"
        if skip_existing and out_dom.exists() and out_dom.stat().st_size > 0:
            continue

        cmd = ["hmmsearch", "--cpu", str(cpus)]
        if evalue:
            cmd += ["-E", str(evalue)]
        if domevalue:
            cmd += ["--domE", str(domevalue)]
        cmd += ["--domtblout", str(out_dom), str(hmm), str(protein_faa)]

        try:
            sh(cmd, quiet=False)
        except subprocess.CalledProcessError as e:
            print(f"[err] hmmsearch failed for {protein_faa.name} × {model} with code {e.returncode}", file=sys.stderr)

def main():
    ap = argparse.ArgumentParser(
        description="Run hmmsearch vs protein FASTAs, writing sample,hmm,hits,total_genes,hit_headers. "
                    "Sample name is taken from the FASTA basename."
    )
    ap.add_argument("--hmm-dir", required=True, help="directory with *.hmm")
    ap.add_argument("--genomes-dir", required=True, help="directory with protein FASTA files")
    ap.add_argument("--prot-glob", default="*proteins.faa", help="glob for protein FASTA filenames (comma-separated for multiple patterns)")
    ap.add_argument("--domtbl-subdir", default="domtbl", help="subfolder to store domtbl outputs alongside each FASTA")
    ap.add_argument("--cpus", type=int, default=8, help="threads for hmmsearch")
    ap.add_argument("--evalue", type=str, default=None, help="sequence E-value cutoff for reporting, e.g. 1e-5")
    ap.add_argument("--domevalue", type=str, default=None, help="domain E-value cutoff for reporting, e.g. 1e-5")
    ap.add_argument("--min-bits", type=float, default=None, help="optional min bit score when counting")
    ap.add_argument("--max-ie", type=float, default=None, help="optional max i-Evalue when counting")
    ap.add_argument("--unique-per-seq", action="store_true", help="count at most one hit per sequence per HMM")
    ap.add_argument("--skip-existing", action="store_true", help="skip running hmmsearch if domtbl already exists")
    ap.add_argument("--all-hits-per-protein", action="store_true",
                    help="Count all HMM hits for a protein that pass filters. "
                         "Default is to count only the single best-scoring HMM hit per protein.")
    ap.add_argument("--run-prodigal", action="store_true",
                    help="Run Prodigal on nucleotide FASTA inputs to create proteins before hmmsearch. "
                         "Auto-enabled when no protein FASTAs are found but FASTA/FNA files exist.")
    ap.add_argument("--prodigal-out-dir", default=None,
                    help="Directory to write Prodigal protein FASTAs (default: <genomes-dir>/prodigal_proteins).")
    ap.add_argument("--prodigal-mode", default="single", choices=["single", "meta"],
                    help="Prodigal -p mode (default: single).")
    ap.add_argument("--out", required=True, help="CSV path for output counts (sample,hmm,hits,total_genes,hit_headers)")
    ap.add_argument("--cutoffs", default=None, help="TSV with columns: hmm<TAB>cutoff (may include many models)")
    ap.add_argument("--fallback-domtbl-root", default=None,
                    help="Optional root containing per-sample domtbl folders: <root>/<sample>/domtbl/*.domtbl")
    args = ap.parse_args()

    check_bin("hmmsearch")

    hmm_dir = Path(args.hmm_dir)
    proteins_dir = Path(args.genomes_dir)
    prot_globs = tuple(p.strip() for p in args.prot_glob.split(",")) if args.prot_glob else ("*proteins.faa",)

    rows = annotate_proteins(
        proteins_path=proteins_dir,
        out_csv=Path(args.out),
        hmm_dir=hmm_dir,
        cutoffs_path=Path(args.cutoffs) if args.cutoffs else None,
        cpus=args.cpus,
        evalue=args.evalue,
        domevalue=args.domevalue,
        min_bits=args.min_bits,
        max_ie=args.max_ie,
        unique_per_seq=args.unique_per_seq,
        all_hits_per_protein=args.all_hits_per_protein,
        skip_existing=args.skip_existing,
        domtbl_subdir=args.domtbl_subdir,
        fallback_domtbl_root=args.fallback_domtbl_root,
        prot_glob=prot_globs,
        run_prodigal=args.run_prodigal,
        prodigal_out_dir=Path(args.prodigal_out_dir) if args.prodigal_out_dir else None,
        prodigal_mode=args.prodigal_mode,
    )

    if rows:
        print(f"wrote {args.out}")
