#!/usr/bin/env python3
"""
make_circos_for_one_genome.py

Build one Circos project for a single genome (sample) so that:
- all contigs from that genome appear on the SAME plot (unless --only-hit-contigs),
- each HMM subunit is drawn in ONE consistent (lowercase) color,
- only HMMs with >=1 hit for that genome are included,
- colors may be user-specified via --hmm-colors TSV or auto-assigned.

MODIFIED VERSION:
- In --operon mode, calculates operon completeness for each potential operon instance.
- Creates multiple, concentric tracks (rings) in the plot.
- The most complete operons are drawn on the outermost track.
- Less complete operons are drawn on progressively smaller, inner tracks.
- NEW: Can group HMM colors by pathway using --pathway-map.

Usage example
-------------
usage:
python /home/juneq/BTEX-HMMs/btexhmm/circos.py \
  --hits /home/juneq/BTEX_test/test_output_2/btex_hmm_summary.csv \
  --outdir /home/juneq/BTEX_test/test_output_circos \
  --dna /home/juneq/BTEX-HMMs/btexhmm/test_genomes/Aromatoleum_bremense_PbN1T.fna \
  --genome "Aromatoleum_bremense_PbN1T" 

# TODO: check operon functionality 
  --operon-defs /home/juneq/hmm/archetypes/hmm_cutoffs/operons_dist.tsv \

If --contig-lengths is omitted, supply --genomes-dir (and optionally --fasta-glob)
to auto-write contig_length.tsv under the Circos output directory using find_contig_length logic.

The Circos etc directory is auto-detected using the circos executable on PATH,
following: ETCDIR=$(dirname "$(dirname "$(which circos)")")/etc

For stable pathway coloring across runs, supply a pathway map TSV with columns:
  hmm<TAB>pathway<TAB>color
The color column is optional; when present, all HMMs in that pathway share the
specified color (hex #RRGGBB or r,g,b).

If --pathway-map is not provided, circos.py will auto-use btexhmm/data/pathway_map.tsv
when present (relative to this script).

By default, the script renders the plot via `circos -conf circos.conf` in the output
directory; use --no-render to skip rendering.
"""
import argparse
import csv
import os
import re
import hashlib
import sys
import shutil
import subprocess
import itertools
from pathlib import Path
from collections import OrderedDict, defaultdict

try:
    from btexhmm.find_contig_length import (
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


# Auto palette (rgb triplets as strings)
AUTO_RGBS = [
    "0,114,178","213,94,0","0,158,115","204,121,167","86,180,233","230,159,0",
    "240,228,66","173,73,74","0,0,0","102,0,204","0,153,136","255,127,14",
    "31,119,180","44,160,44","214,39,40","148,103,189","140,86,75","227,119,194",
    "127,127,127","188,189,34","23,190,207","255,152,150","197,176,213","196,156,148"
]

def parse_args():
    ap = argparse.ArgumentParser(description="Generate Circos plot files for HMM hits on a single genome.")
    ap.add_argument(
        "--hits",
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
    ap.add_argument("--genome", required=True, help="Sample/genome name to render")
    ap.add_argument("--hmm-colors", help="Optional TSV: hmm<TAB>color (rgb 'r,g,b' or hex '#RRGGBB')")
    ap.add_argument("--only-hit-contigs", action="store_true", help="Include only contigs that have >=1 hit")
    ap.add_argument("--operon", action="store_true", help="Switch to operon-centric plotting mode with layered tracks for completeness.")
    ap.add_argument("--operon-defs", help="Required if --operon is used. TSV with header: operon_id,members,max_span_bp")
    ap.add_argument("--no-render", action="store_true", help="Skip running circos -conf circos.conf after generating the project.")
    ap.add_argument(
        "--dna",
        help=(
            "Path to the genome nucleotide FASTA (.fna/.fasta). "
            "Used to generate contig lengths when --contig-lengths is omitted and to populate GenBank sequences."
        ),
    )
    # NEW: Argument for the pathway map file.
    ap.add_argument(
        "--pathway-map",
        help=(
            "Optional TSV: hmm<TAB>pathway[<TAB>color]. "
            "Groups HMMs by pathway for coloring; a pathway color column can be supplied to keep colors stable across runs."
        ),
    )
    return ap.parse_args()

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

def hashed_rgb(name):
    h = hashlib.md5(name.encode("utf-8")).hexdigest()
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    def bump(x): return min(255, int(0.6*255 + 0.4*x))
    return f"{bump(r)},{bump(g)},{bump(b)}"

# MODIFIED: This function is rewritten to handle pathway-based and case-insensitive coloring.
def build_color_map(hmms_with_hits, pathway_map, pathway_colors, user_colors_path=None):
    """
    Builds a color map for HMMs with a new priority order and case-insensitivity:
    1. User-specified colors from --hmm-colors (highest priority).
    2. Pathway-based colors from --pathway-map.
    3. Unique auto-assigned colors for any remaining HMMs (fallback).
    """
    # 1. Load user-specified colors with lowercase HMM keys for case-insensitive matching.
    user_rgb_lower = {}
    if user_colors_path:
        with open(user_colors_path, newline="") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            for row in rdr:
                user_rgb_lower[row['hmm'].strip().lower()] = normalize_color_token(row['color'])

    # 2. Create a lowercase version of the pathway map for consistent matching.
    pathway_map_lower = {k.lower(): v for k, v in pathway_map.items()}
    pathway_colors_lower = {k.lower(): normalize_color_token(v) for k, v in pathway_colors.items()}

    # 3. Determine the set of pathways and un-mapped HMMs to assign colors to.
    entities_to_color = set()
    for hmm in hmms_with_hits:
        hmm_lower = hmm.lower() # Use a lowercase version for lookups
        if hmm_lower in user_rgb_lower:
            continue  # Skip, this color is already defined by the user.
        
        pathway = pathway_map_lower.get(hmm_lower)
        if pathway:
            entities_to_color.add(pathway)
        else:
            entities_to_color.add(hmm) # Fallback: the entity is the HMM itself.

    # 4. Assign a color to each entity (pathway name or individual HMM name).
    auto_iter = iter(AUTO_RGBS)
    entity_color_map = {
        p_lower: (sanitize_color_name(p_lower), rgb) for p_lower, rgb in pathway_colors_lower.items()
    }
    for entity in sorted(list(entities_to_color), key=lambda x: x.lower()):
        key = entity.lower()
        if key in entity_color_map:
            continue
        try:
            rgb = next(auto_iter)
        except StopIteration:
            rgb = hashed_rgb(entity)
        
        # Color names are also made from lowercase entity names for consistency
        color_name = sanitize_color_name(entity)
        entity_color_map[key] = (color_name, rgb)

    # 5. Build the final HMM -> color map using lowercase lookups.
    final_cmap = {}
    for hmm in hmms_with_hits:
        hmm_lower = hmm.lower()
        if hmm_lower in user_rgb_lower:
            # Priority 1: User-defined color.
            cname = sanitize_color_name(hmm_lower) # Create a unique name
            rgb = user_rgb_lower[hmm_lower]
            final_cmap[hmm] = (cname, rgb)
        else:
            pathway = pathway_map_lower.get(hmm_lower)
            if pathway:
                # Priority 2: Pathway color.
                final_cmap[hmm] = entity_color_map[pathway.lower()]
            else:
                # Priority 3: Fallback unique color for the HMM.
                final_cmap[hmm] = entity_color_map[hmm_lower]
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
        return generate_contig_lengths_table_from_dna(Path(args.dna), args.genome, default_out)

    if not args.dna:
        raise SystemExit("[ERROR] Provide --dna or an existing --contig-lengths file.")

    print(f"[info] Auto-generating contig lengths -> {default_out}", file=sys.stderr)
    return generate_contig_lengths_table_from_dna(Path(args.dna), args.genome, default_out)

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

def main():
    args = parse_args()
    genome = args.genome
    etc_dir = find_circos_etc_dir()
    default_pathway_map = Path(__file__).resolve().parent / "data" / "pathway_map.tsv"
    if not args.pathway_map:
        if default_pathway_map.exists():
            args.pathway_map = str(default_pathway_map)
            print(f"[info] Using default pathway map: {args.pathway_map}", file=sys.stderr)
        else:
            print(f"[warn] No --pathway-map provided and default not found at {default_pathway_map}; colors will be auto-assigned.", file=sys.stderr)

    hits_path = Path(args.hits)
    if not hits_path.exists():
        raise SystemExit(f"[ERROR] hits file not found: {hits_path}")

    sample_dir = Path(args.outdir) / f"{genome}_circos_plot"
    ensure_dir(sample_dir)

    if args.operon and not args.operon_defs:
        raise SystemExit("[ERROR] --operon-defs is required when using the --operon flag.")

    contig_lengths_path = ensure_contig_lengths(args, sample_dir)
    contig_lengths = read_contig_lengths(contig_lengths_path, genome)
    if not contig_lengths:
        raise SystemExit(f"[ERROR] No contig lengths found for genome: {genome}")

    features, hmms_with_hits, contigs_with_hits = [], set(), set()
    with open(args.hits, newline="") as fh:
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

    # NEW: Read the pathway map file if provided.
    pathway_map = {}
    pathway_colors = {}
    if args.pathway_map:
        try:
            with open(args.pathway_map, newline="") as fh:
                rdr = csv.DictReader(fh, delimiter="\t")
                required = {"hmm", "pathway"}
                if not required.issubset(set(rdr.fieldnames or [])):
                    missing = required - set(rdr.fieldnames or [])
                    raise SystemExit(f"[ERROR] Pathway map missing required columns: {', '.join(sorted(missing))}")
                if "color" not in (rdr.fieldnames or []):
                    print(f"[warn] Pathway map {args.pathway_map} has no 'color' column; using auto-colors per pathway.", file=sys.stderr)
                for row in rdr:
                    hmm_val = row.get('hmm', '').strip()
                    pathway_val = row.get('pathway', '').strip()
                    if not hmm_val or not pathway_val:
                        print(f"[warn] Skipping pathway map row with empty hmm/pathway: {row}", file=sys.stderr)
                        continue
                    pathway_map[hmm_val] = pathway_val
                    color_val = row.get('color', '').strip()
                    if color_val:
                        pathway_colors[pathway_val] = color_val
        except FileNotFoundError:
            raise SystemExit(f"[ERROR] Pathway map file not found: {args.pathway_map}")

    # MODIFIED: Pass the pathway map to the color building function.
    cmap = build_color_map(hmms_with_hits, pathway_map, pathway_colors, args.hmm_colors) if hmms_with_hits else {}

    if args.operon:
        operon_instances = find_operon_instances(features, args.operon_defs)
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
    if args.no_render:
        print(f"[info] Rendering skipped (--no-render set). To render manually: cd '{sample_dir}' && circos -conf circos.conf")
        return

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

if __name__ == "__main__":
    main()
