#!/usr/bin/env python3
"""
make_circos_for_one_genome.py

Build one Circos project for a single genome (sample) so that:
- all contigs from that genome appear on the SAME plot (unless --only-hit-contigs),
- each HMM subunit is drawn in ONE consistent (lowercase) color,
- only HMMs with >=1 hit for that genome are included,
- colors are taken from the bundled pathway map TSV.

MODIFIED VERSION:
- In --operon mode, calculates operon completeness for each potential operon instance.
- Creates multiple, concentric tracks (rings) in the plot.
- The most complete operons are drawn on the outermost track.
- Less complete operons are drawn on progressively smaller, inner tracks.
- Colors can be grouped by pathway using the bundled default pathway map when present.

Usage example
-------------
usage:
python /home/juneq/BTEX-HMMs/visualization_scripts/circos.py \
  --hmmscan /home/juneq/BTEX_test/test_output_2/btex_hmm_summary.csv \
  --outdir /home/juneq/BTEX_test/test_output_circos \
  --dna /home/juneq/BTEX-HMMs/btexhmm/test_genomes/Aromatoleum_bremense_PbN1T.fna \
  --sample "Aromatoleum_bremense_PbN1T"

If --contig-lengths is omitted, supply --genomes-dir (and optionally --fasta-glob)
to auto-write contig_length.tsv under the Circos output directory using find_contig_length logic.

The Circos etc directory is auto-detected using the circos executable on PATH,
following: ETCDIR=$(dirname "$(dirname "$(which circos)")")/etc

For stable pathway coloring across runs, the script auto-uses
`btexhmm/data/pathway_map.tsv` when present.

The script renders the plot via `circos -conf circos.conf` in the output directory.
"""
import csv
import os
import re
import sys
import shutil
import subprocess
from pathlib import Path
from collections import OrderedDict, defaultdict

try:
    from visualization_scripts.find_contig_length import (
        FIELDNAMES as CONTIG_TSV_FIELDS,
        read_fai,
        read_fasta_lengths,
    )
except ImportError:
    # Fallback for running as a standalone script from the source tree.
    from find_contig_length import (
        FIELDNAMES as CONTIG_TSV_FIELDS,
        read_fai,
        read_fasta_lengths,
    )

try:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
except ImportError:
    print("[ERROR] Biopython is required. Please install it: pip install biopython", file=sys.stderr)
    sys.exit(1)

def find_circos_etc_dir() -> Path:
    circos_exe = shutil.which("circos")
    if not circos_exe:
        raise SystemExit("[ERROR] circos executable not found in PATH; cannot locate its etc directory.")

    etc_dir = Path(circos_exe).resolve().parent.parent / "etc"
    if not etc_dir.exists():
        raise SystemExit(
            f"[ERROR] Expected Circos etc directory not found at {etc_dir}. "
            "Ensure circos is installed and available on PATH."
        )
    return etc_dir

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def hex_to_rgb_str(hexcode):
    h = hexcode.strip().lstrip("#")
    if len(h) != 6: raise ValueError(f"Bad hex color: {hexcode}")
    r = int(h[0:2], 16); g = int(h[2:4], 16); b = int(h[4:6], 16)
    return f"{r},{g},{b}"

def sanitize_color_name(token):
    return re.sub(r"[^a-z0-9]+", "_", token.lower())

def normalize_color_token(token):
    t = token.strip()
    if "," in t:
        parts = [p.strip() for p in t.split(",")]
        if len(parts) != 3: raise ValueError
        r,g,b = (int(parts[0]), int(parts[1]), int(parts[2])); [v for v in (r,g,b) if v < 0 or v > 255]
        return f"{r},{g},{b}"
    if t.startswith("#"): return hex_to_rgb_str(t)
    raise ValueError(f"Unsupported color token: {token}")

def build_color_map(hmms_with_hits, pathway_map, pathway_colors):
    """
    Builds an HMM -> color map using only the bundled pathway map color column.
    """
    pathway_map_lower = {k.lower(): v for k, v in pathway_map.items()}
    pathway_colors_lower = {k.lower(): normalize_color_token(v) for k, v in pathway_colors.items()}
    final_cmap = {}
    for hmm in hmms_with_hits:
        hmm_lower = hmm.lower()
        pathway = pathway_map_lower.get(hmm_lower)
        if not pathway:
            raise SystemExit(f"[ERROR] No pathway mapping found for HMM: {hmm}")
        rgb = pathway_colors_lower.get(pathway.lower())
        if not rgb:
            raise SystemExit(
                f"[ERROR] No color defined in pathway_map.tsv for pathway '{pathway}' (HMM: {hmm})"
            )
        final_cmap[hmm] = (sanitize_color_name(pathway), rgb)
    return final_cmap

def read_contig_lengths(path, genome):
    out = OrderedDict()
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            if row["sample"] == genome: out[row["contig"]] = int(row["length"])
    return out

def contig_lengths_from_dna(dna_path: Path):
    fai = dna_path.with_suffix(dna_path.suffix + ".fai")
    if fai.exists():
        return list(read_fai(fai))
    print(f"[info] {fai.name} missing; parsing FASTA to compute lengths for {dna_path.name}", file=sys.stderr)
    return list(read_fasta_lengths(dna_path))

def generate_contig_lengths_table_from_dna(dna_path: Path, genome: str, out_path: Path) -> Path:
    if not dna_path.exists():
        raise SystemExit(f"[ERROR] DNA FASTA not found: {dna_path}")

    if dna_path.suffix.lower() not in {".fna", ".fasta", ".fa"}:
        print(f"[warn] DNA file has unexpected extension ({dna_path.suffix}); continuing anyway.", file=sys.stderr)

    entries = contig_lengths_from_dna(dna_path)
    if not entries:
        raise SystemExit(f"[ERROR] No contig lengths gathered from {dna_path}")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CONTIG_TSV_FIELDS, delimiter="\t")
        writer.writeheader()
        for contig, length in entries:
            writer.writerow({"sample": genome, "contig": contig, "length": length})
    print(f"[info] wrote {len(entries)} contig entries to {out_path}", file=sys.stderr)
    return out_path

def ensure_contig_lengths(args, output_dir: Path) -> Path:
    """
    Return a Path to a contig-length TSV under output_dir, generating it from --dna when needed.
    If --contig-lengths is provided elsewhere, it is copied into output_dir for use.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    default_out = output_dir / "contig_length.tsv"

    if args.contig_lengths:
        candidate = Path(args.contig_lengths)
        if candidate.exists():
            if candidate.resolve() != default_out.resolve():
                shutil.copy(candidate, default_out)
                print(f"[info] Copied contig lengths to {default_out}", file=sys.stderr)
            else:
                print(f"[info] Using contig lengths at {default_out}", file=sys.stderr)
            return default_out
        if not args.dna:
            raise SystemExit(
                f"[ERROR] contig lengths file not found: {candidate}. "
                "Provide --dna so it can be generated automatically."
            )
        print(f"[info] contig lengths not found at {candidate}; generating -> {default_out}", file=sys.stderr)
        return generate_contig_lengths_table_from_dna(Path(args.dna), args.sample, default_out)

    if not args.dna:
        raise SystemExit("[ERROR] Provide --dna or an existing --contig-lengths file.")

    print(f"[info] Auto-generating contig lengths -> {default_out}", file=sys.stderr)
    return generate_contig_lengths_table_from_dna(Path(args.dna), args.sample, default_out)

def normalize_contig_id(token):
    if not token: return ""
    name = token.strip().lstrip(">").strip()
    if " #" in name: name = name.split(" #", 1)[0].strip()
    if "|" in name: name = name.split("|")[-1]
    return name.strip()

def parse_hit_headers(sample, hmm, headers):
    out = []
    if not headers: return out
    for raw in headers.split("|"):
        if not raw.strip(): continue
        left = raw.split(" #", 1)[0] if " #" in raw else raw
        contig = normalize_contig_id(left)
        # Prodigal-style headers append an ORF index as a trailing "_<number>"
        # to the contig name (e.g., NZ_CP059467.1_1976). Remove only that suffix.
        if " #" in raw:
            contig = re.sub(r"_\d+$", "", contig)
        if not contig: continue
        m = re.findall(r"#\s*(\d+)", raw)
        if len(m) < 2: m = re.findall(r"\b(\d+)\b", raw)
        if len(m) < 2: continue
        a, b = int(m[0]), int(m[1])
        out.append({"sample": sample, "hmm": hmm, "contig": contig, "start": min(a,b), "end": max(a,b)})
    return out

def find_operon_instances(features, operon_defs_path):
    operon_defs = {}
    with open(operon_defs_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            op_id = row['operon_id']
            members = {m.strip() for m in row['members'].split(',')}
            max_span = int(row['max_span_bp'])
            operon_defs[op_id] = {'members': members, 'max_span': max_span}

    hits_by_contig = defaultdict(list)
    for f in features:
        hits_by_contig[f["contig"]].append(f)

    found_operons = []
    for operon_id, op_def in operon_defs.items():
        required_members = op_def['members']
        max_span = op_def['max_span']

        for contig, contig_hits in hits_by_contig.items():
            member_hmms_lower = {m.lower() for m in required_members}
            member_hits = [h for h in contig_hits if h['hmm'].lower() in member_hmms_lower]
            
            if not member_hits:
                continue

            member_hits.sort(key=lambda x: x['start'])
            best_cluster_on_contig = None

            for i in range(len(member_hits)):
                start_hit = member_hits[i]
                current_cluster_hits = [start_hit]
                
                for j in range(i + 1, len(member_hits)):
                    next_hit = member_hits[j]
                    if (next_hit['end'] - start_hit['start']) <= max_span:
                        current_cluster_hits.append(next_hit)
                    else:
                        break

                cluster_hmms_found = {h['hmm'].lower() for h in current_cluster_hits}
                completeness = len(cluster_hmms_found) / len(required_members)
                span = current_cluster_hits[-1]['end'] - current_cluster_hits[0]['start']

                is_better = False
                if best_cluster_on_contig is None:
                    is_better = True
                else:
                    if completeness > best_cluster_on_contig['completeness']:
                        is_better = True
                    elif completeness == best_cluster_on_contig['completeness'] and span < best_cluster_on_contig['span']:
                        is_better = True
                
                if is_better:
                    best_cluster_on_contig = {
                        'operon_id': operon_id,
                        'contig': contig,
                        'start': current_cluster_hits[0]['start'],
                        'end': current_cluster_hits[-1]['end'],
                        'span': span,
                        'completeness': completeness,
                        'hits': current_cluster_hits
                    }

            if best_cluster_on_contig:
                completeness_pct = best_cluster_on_contig['completeness'] * 100
                print(f"[info] Found instance of '{operon_id}' on {contig} "
                      f"(Completeness: {completeness_pct:.1f}%, Span: {best_cluster_on_contig['span']} bp)", file=sys.stderr)
                found_operons.append(best_cluster_on_contig)

    return found_operons


def write_colors_conf(path, cmap):
    with open(path, "w") as fh:
        fh.write("<colors>\n")
        # MODIFIED: Use cmap.values() to handle cases where multiple HMMs share a color name
        # This prevents duplicate color name definitions in the conf file.
        unique_colors = {val[0]: val[1] for val in cmap.values()}
        for cname, rgb in unique_colors.items():
            fh.write(f"  {cname} = {rgb}\n")
        fh.write("</colors>\n")

def write_hmm_colors_tsv(path, cmap):
    with open(path, "w") as fh:
        fh.write("hmm\tcolor_name\trgb\n")
        for hmm, (cname, rgb) in cmap.items(): fh.write(f"{hmm}\t{cname}\t{rgb}\n")

def write_karyotype(path, contig_lengths):
    with open(path, "w") as fh:
        for contig, length in contig_lengths.items():
            fh.write(f"chr - {contig} {contig} 0 {length} grey\n")

def write_genes(path, features, cmap):
    with open(path, "w") as fh:
        for f in features:
            if f["hmm"] in cmap:
                cname = cmap[f["hmm"]][0]
                fh.write(f"{f['contig']} {f['start']} {f['end']} color={cname}\n")

def write_gene_labels(path, features):
    seen_labels = set()
    with open(path, "w") as fh:
        for f in features:
            label_key = (f['contig'], f['start'], f['end'], f['hmm'])
            if label_key not in seen_labels:
                fh.write(f"{f['contig']} {f['start']} {f['end']} {f['hmm']}\n")
                seen_labels.add(label_key)

def write_links_placeholder(path):
    open(path, "w").close()

def write_genbank_file(path, features, contig_lengths, fna_path=None):
    seq_records = {}
    if fna_path and fna_path.exists():
        try:
            for record in SeqIO.parse(fna_path, "fasta"):
                seq_records[record.id] = record
            print(f"[info] Loaded {len(seq_records)} sequences from {fna_path}", file=sys.stderr)
        except Exception as e:
            print(f"[warn] Could not read FASTA file {fna_path}: {e}. Falling back to 'N' sequence.", file=sys.stderr)
            seq_records = {}

    records = []
    for contig, length in contig_lengths.items():
        if contig in seq_records:
            record = seq_records[contig]
        else:
            seq = Seq("N" * length)
            record = SeqRecord(seq, id=contig, name=contig, description="")
        
        record.annotations["molecule_type"] = "DNA"
        
        contig_features = []
        for f in features:
            if f["contig"] == contig:
                feature = SeqFeature(
                    FeatureLocation(f["start"] - 1, f["end"]),
                    type="CDS",
                    qualifiers={"locus_tag": f["hmm"], "product": f["hmm"]}
                )
                contig_features.append(feature)
        
        record.features.extend(sorted(contig_features, key=lambda f: f.location.start))
        records.append(record)
    
    with open(path, "w") as fh:
        SeqIO.write(records, fh, "genbank")


def write_circos_conf(path, etcdir, tracks_info=None, operon_regions=None, operon_flag=False):
    plots_block = ""
    zooms_block = ""
    ideogram_radius = 0.85

    if operon_flag and tracks_info:
        outermost_track_r1 = ideogram_radius - 0.05
        track_thickness = 0.05
        track_gap = 0.07

        for i, track in enumerate(tracks_info):
            r1_highlight = outermost_track_r1 - i * (track_thickness + track_gap)
            r0_highlight = r1_highlight - track_thickness

            label_padding = 0.005
            r0_label = r1_highlight + label_padding
            r1_label = r1_highlight + track_gap - label_padding

            plots_block += f"""
  # Track {i+1} Highlights
  <plot>
    type   = highlight
    file   = {track['highlight_file']}
    r0     = {r0_highlight:.3f}r
    r1     = {r1_highlight:.3f}r
    stroke_thickness = 0p
    z = {10 - i}
    <rules>
      <rule>
        condition  = 1
        fill_color = eval(var(color))
      </rule>
    </rules>
  </plot>
"""
            plots_block += f"""
  # Track {i+1} Labels
  <plot>
    type         = text
    file         = {track['label_file']}
    r0           = {r0_label:.3f}r
    r1           = {r1_label:.3f}r
    label_size   = 13p
    label_font   = light
    padding      = 0p
    rpadding     = 0p
    show_links   = no
    label_snuggle = yes
    max_snuggle_distance = 5r
    <rules>
      <rule>
        importance = {100 - i*10}
        condition = 1
        color = black
      </rule>
    </rules>
  </plot>
"""
        if operon_regions:
            zooms_block = "<zooms>\n"
            for r in operon_regions:
                zooms_block += f"""  <zoom>
    chr   = {r['contig']}
    start = {r['start']}b
    end   = {r['end']}b
    scale = 100
  </zoom>
"""
            zooms_block += "</zooms>\n"

    else:
        plots_block = """
  <plot>
    type   = highlight
    file   = genes.txt
    r0     = 0.70r
    r1     = 0.82r
    stroke_thickness = 0p
    <rules>
      <rule>
        condition  = 1
        fill_color = eval(var(color))
      </rule>
    </rules>
  </plot>
  <plot>
    type         = text
    file         = gene_labels.txt
    r0           = 0.84r
    r1           = 0.90r
    label_size   = 13p
    label_font   = default
    padding      = 0p
    rpadding     = 0p
    show_links   = no
    label_snuggle = yes
    max_snuggle_distance = 5r
    <rules>
      <rule>
        importance = 100
        condition = 1
        color = black
      </rule>
    </rules>
  </plot>
"""

    with open(path, "w") as fh:
        fh.write(f"""
<<include {etcdir}/colors_fonts_patterns.conf>>
<<include {etcdir}/housekeeping.conf>>
<image>
<<include {etcdir}/image.conf>>
</image>
<<include colors.conf>>
karyotype = karyotype.txt
chromosomes_units = 1000000
<ideogram>
  <spacing>
    default = 0.005r
  </spacing>
  radius    = {ideogram_radius}r
  thickness = 40p
  fill      = yes
  show_label       = yes
  label_font       = default
  label_radius     = dims(ideogram,radius) + 0.075r
  label_size       = 30
</ideogram>
show_ticks          = yes
show_tick_labels    = yes
<ticks>
  radius           = dims(ideogram,radius_outer)
  color            = black
  thickness        = 2p
  <tick>
    spacing        = 0.5u
    size           = 15p
    show_label     = yes
    label_size     = 20p
    label_offset   = 10p
    label_format   = eval(sprintf("%d Mb", var(value)))
  </tick>
  <tick>
    spacing        = 0.1u
    size           = 8p
    show_label     = no
  </tick>
</ticks>
<plots>
{plots_block}
</plots>
{zooms_block}
""")

def run(args):
    genome = args.sample
    etc_dir = find_circos_etc_dir()
    repo_root = Path(__file__).resolve().parents[1]
    default_pathway_map = repo_root / "btexhmm" / "data" / "pathway_map.tsv"
    default_operon_defs = repo_root / "btexhmm" / "data" / "operons.tsv"
    if default_pathway_map.exists():
        print(f"[info] Using default pathway map: {default_pathway_map}", file=sys.stderr)
    else:
        raise SystemExit(f"[ERROR] Default pathway map not found: {default_pathway_map}")

    hmmscan_path = Path(args.hmmscan)
    if not hmmscan_path.exists():
        raise SystemExit(f"[ERROR] hmmscan file not found: {hmmscan_path}")

    sample_dir = Path(args.outdir) / f"{genome}_circos_plot"
    ensure_dir(sample_dir)

    operon_defs_path = None
    if args.operon:
        operon_defs_path = default_operon_defs
        if not operon_defs_path.exists():
            raise SystemExit(
                f"[ERROR] Operon definitions file not found: {operon_defs_path}. "
                "Add the bundled btexhmm/data/operons.tsv file."
            )

    contig_lengths_path = ensure_contig_lengths(args, sample_dir)
    contig_lengths = read_contig_lengths(contig_lengths_path, genome)
    if not contig_lengths:
        raise SystemExit(f"[ERROR] No contig lengths found for genome: {genome}")

    features, hmms_with_hits, contigs_with_hits = [], set(), set()
    with open(args.hmmscan, newline="") as fh:
        rdr = csv.DictReader(fh)
        for row in rdr:
            if row.get("sample") != genome or int(row.get("hits", "0")) < 1:
                continue
            hmm, headers = row["hmm"], row.get("hit_headers", "")
            parsed = parse_hit_headers(genome, hmm, headers)
            if parsed:
                features.extend(parsed)
                hmms_with_hits.add(hmm)
                for rec in parsed: contigs_with_hits.add(rec["contig"])

    if args.only_hit_contigs:
        if not contigs_with_hits:
            raise SystemExit(f"[ERROR] --only-hit-contigs set but no hits found for {genome}.")
        filtered = OrderedDict((c, l) for c, l in contig_lengths.items() if c in contigs_with_hits or normalize_contig_id(c) in contigs_with_hits)
        if filtered: contig_lengths = filtered
        else: print(f"[warn] No contig names matched hit headers for {genome}; keeping all contigs.", file=sys.stderr)

    pathway_map = {}
    pathway_colors = {}
    if default_pathway_map.exists():
        try:
            with open(default_pathway_map, newline="") as fh:
                rdr = csv.DictReader(fh, delimiter="\t")
                required = {"hmm", "pathway", "color"}
                if not required.issubset(set(rdr.fieldnames or [])):
                    missing = required - set(rdr.fieldnames or [])
                    raise SystemExit(f"[ERROR] Pathway map missing required columns: {', '.join(sorted(missing))}")
                for row in rdr:
                    hmm_val = row.get('hmm', '').strip()
                    pathway_val = row.get('pathway', '').strip()
                    color_val = row.get('color', '').strip()
                    if not hmm_val or not pathway_val or not color_val:
                        raise SystemExit(
                            f"[ERROR] Pathway map row must include hmm, pathway, and color values: {row}"
                        )
                    pathway_map[hmm_val] = pathway_val
                    pathway_colors[pathway_val] = color_val
        except FileNotFoundError:
            raise SystemExit(f"[ERROR] Pathway map file not found: {default_pathway_map}")

    cmap = build_color_map(hmms_with_hits, pathway_map, pathway_colors) if hmms_with_hits else {}

    if args.operon:
        operon_instances = find_operon_instances(features, operon_defs_path)
        if not operon_instances:
            print(f"[warn] No operons found for {genome}; Circos plot may be empty or minimal.", file=sys.stderr)
            write_circos_conf(sample_dir / "circos.conf", etc_dir)
        else:
            completeness_levels = sorted({inst['completeness'] for inst in operon_instances}, reverse=True)
            level_to_track_idx = {level: i for i, level in enumerate(completeness_levels)}

            tracks_data = defaultdict(list)
            for inst in operon_instances:
                track_idx = level_to_track_idx[inst['completeness']]
                tracks_data[track_idx].extend(inst['hits'])

            tracks_info_for_conf = []
            for track_idx in sorted(tracks_data.keys()):
                seen_hits = set()
                unique_hits = []
                for hit in tracks_data[track_idx]:
                    hit_tuple = (hit['contig'], hit['start'], hit['end'], hit['hmm'])
                    if hit_tuple not in seen_hits:
                        unique_hits.append(hit)
                        seen_hits.add(hit_tuple)

                highlight_filename = f"genes_track_{track_idx}.txt"
                label_filename = f"gene_labels_track_{track_idx}.txt"

                write_genes(sample_dir / highlight_filename, unique_hits, cmap)
                write_gene_labels(sample_dir / label_filename, unique_hits)

                tracks_info_for_conf.append({
                    "highlight_file": highlight_filename,
                    "label_file": label_filename
                })

            write_circos_conf(
                sample_dir / "circos.conf",
                etc_dir,
                tracks_info=tracks_info_for_conf,
                operon_regions=operon_instances,
                operon_flag=True
            )

    else:
        write_genes(sample_dir / "genes.txt", features, cmap)
        write_gene_labels(sample_dir / "gene_labels.txt", features)
        write_circos_conf(sample_dir / "circos.conf", etc_dir)

        gbk_path = sample_dir / "gene_hits.gbk"
        dna_file_path = Path(args.dna) if args.dna else None
        write_genbank_file(gbk_path, features, contig_lengths, fna_path=dna_file_path)
        print(f"[OK] GenBank hits file for {genome} -> {gbk_path}")

    write_colors_conf(sample_dir / "colors.conf", cmap)
    write_hmm_colors_tsv(sample_dir / "hmm_colors.tsv", cmap)
    write_karyotype(sample_dir / "karyotype.txt", contig_lengths)
    write_links_placeholder(sample_dir / "links.txt")

    print(f"[OK] Circos project for {genome} -> {sample_dir}")

    try:
        result = subprocess.run(
            ["circos", "-conf", "circos.conf"],
            cwd=sample_dir,
            capture_output=True,
            text=True,
            check=True,
        )
        sys.stdout.write(result.stdout)
        sys.stderr.write(result.stderr)
        print(f"[OK] Circos rendering complete in {sample_dir}")
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stdout or "")
        sys.stderr.write(e.stderr or "")
        raise SystemExit(f"[ERROR] circos rendering failed with exit code {e.returncode}. See logs above.")


def main(args):
    return run(args)

if __name__ == "__main__":
    from btexhmm.circos_cli import main as cli_main

    cli_main()
