#!/usr/bin/env python3
import argparse, re
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Expected input columns in both hmm CSVs:
#   sample, hmm, hits, total_genes
#
# BTEX map TSV columns:
#   subunit, ids, rule
#
# BTEX pathway map TSV columns:
#   subunit, pathway
#
# PFAM map TSV columns:
#   subunit, pfam_ids, rule[, pfam_names]
#
# KO map TSV columns:
#   pathway, ko_ids, display_label

# usage
# python /home/juneq/hmm/scripts/pfam/benchmark_7.py \
#   --btex-hmm /home/juneq/hmm/validation_genomes/atested/abtested_results/hmmsearch_summary_all_hits.csv \
#   --pfam-hmm /home/juneq/hmm/validation_genomes/atested/abtested_results/pfam_hmmsearch_summary.csv \
#   --btex-pathway-map /home/juneq/hmm/archetypes/hmm_cutoffs/pathway_map_aryl_benz_deleted.tsv \
#   --pfam-map /home/juneq/hmm/archetypes/hmm_cutoffs/master_pfam_map.tsv \
#   --kofam-hits-dir /home/juneq/hmm/validation_genomes/atested/genomes \
#   --ko-map /home/juneq/hmm/archetypes/hmm_cutoffs/ko_map.tsv \
#   --out-svg /home/juneq/hmm/validation_genomes/atested/figures/ko_output_percent_1119.png \
#   --out-tsv /home/juneq/hmm/validation_genomes/atested/figures/ko_output_percent_1119.matrix.tsv \
#   --value-mode percent \
#   --include-pfam-ko \
#   --out-order /home/juneq/hmm/validation_genomes/atested/figures/ko_output_1119.samples_order.txt

def _split_tokens(s: str):
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return []
    return [t.strip() for t in re.split(r"[,\s]+", str(s)) if t.strip()]

def read_hmm_csv(path: Path, is_pfam: bool = False) -> pd.DataFrame:
    df = pd.read_csv(path)
    need = {"sample","hmm","hits","total_genes"}
    if not need.issubset(df.columns):
        raise SystemExit(f"{path} is missing required columns: {need}")
    df["sample"] = df["sample"].astype(str)
    df["hmm"] = df["hmm"].astype(str)
    if is_pfam:
        df["hmm"] = df["hmm"].str.split(".").str[0]
    df["hits"] = pd.to_numeric(df["hits"], errors="coerce").fillna(0.0)
    df["total_genes"] = pd.to_numeric(df["total_genes"], errors="coerce").fillna(0.0)
    return df

def load_btex_pathway_map(path: Path):
    try:
        df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1])
        headerless = True
    except Exception:
        df = pd.read_csv(path, sep="\t")
        if not {"subunit", "pathway"}.issubset(df.columns):
            raise SystemExit("BTEX pathway map must have two columns: subunit, pathway")
        headerless = False

    sub2path = {}
    path_order = []
    if headerless:
        for t in df.itertuples(index=False, name=None):
            sub = str(t[0]).strip()
            pth = str(t[1]).strip()
            if not sub or not pth:
                continue
            sub2path[sub] = pth
            if pth not in path_order:
                path_order.append(pth)
    else:
        for r in df.itertuples(index=False):
            sub = str(getattr(r, "subunit")).strip()
            pth = str(getattr(r, "pathway")).strip()
            if not sub or not pth:
                continue
            sub2path[sub] = pth
            if pth not in path_order:
                path_order.append(pth)
    return sub2path, path_order

def load_pfam_map(path: Path, name_joiner: str = "+"):
    df = pd.read_csv(path, sep="\t")
    need = {"subunit","pfam_ids","rule"}
    if not need.issubset(df.columns):
        raise SystemExit("PFAM map must have columns: subunit, pfam_ids, rule")
    has_names = "pfam_names" in df.columns

    order, sub2pfams, sub2label = [], {}, {}
    for r in df.itertuples(index=False):
        sub = str(getattr(r, "subunit")).strip()
        if not sub:
            continue
        pfids = [re.split(r"\.", p.strip())[0] for p in _split_tokens(getattr(r, "pfam_ids"))]
        if sub not in sub2pfams:
            order.append(sub)
            sub2pfams[sub] = []
        sub2pfams[sub].extend(pfids)

        if has_names:
            raw_names = _split_tokens(getattr(r, "pfam_names"))
            if raw_names and len(raw_names) == len(pfids):
                label = name_joiner.join(raw_names)
            else:
                label = name_joiner.join(pfids)
        else:
            label = name_joiner.join(pfids)
        label = re.sub(r"\s+", "_", label).strip()
        sub2label[sub] = label

    for k in sub2pfams:
        sub2pfams[k] = sorted(set(sub2pfams[k]))
    return order, sub2pfams, sub2label

def load_ko_map(path: Path):
    df = pd.read_csv(path, sep="\t")
    need = {"pathway", "ko_ids", "display_label"}
    if not need.issubset(df.columns):
        raise SystemExit(f"KO map must have columns: {need}")
    order, pathway2kos, pathway2label = [], {}, {}
    for r in df.itertuples(index=False):
        pathway = str(getattr(r, "pathway")).strip()
        if not pathway:
            continue
        raw_kos = _split_tokens(getattr(r, "ko_ids"))
        kos = []
        for k in raw_kos:
            k2 = str(k).strip().upper()
            if re.fullmatch(r"\d{5}", k2):
                k2 = "K" + k2
            kos.append(k2)
        label = str(getattr(r, "display_label")).strip()
        if pathway not in pathway2kos:
            order.append(pathway)
            pathway2kos[pathway] = []
        pathway2kos[pathway].extend(kos)
        pathway2label[pathway] = label
    return order, pathway2kos, pathway2label

def _first_faa_basename(genome_dir: Path) -> str | None:
    for pat in ("*.faa", "*.fa", "*.fasta", "*protein*.fa*", "*/*.faa"):
        hits = list(genome_dir.glob(pat))
        if hits:
            return hits[0].stem
    return None

def read_kofam_hits(kofam_dir: Path) -> pd.DataFrame:
    all_dfs = []
    hit_files = list(kofam_dir.glob("*/*kofam_btex_hits.tsv"))
    for f in hit_files:
        genome_dir = f.parent
        sample_name = _first_faa_basename(genome_dir) or genome_dir.name
        try:
            df = pd.read_csv(f, sep="\t")
        except (pd.errors.EmptyDataError, FileNotFoundError):
            continue

        cols = {c.lower(): c for c in df.columns}
        ko_col = cols.get("ko") or cols.get("ko_id") or cols.get("k")
        prot_col = cols.get("protein_id") or cols.get("protein") or cols.get("gene") or cols.get("query")

        if ko_col is None or prot_col is None:
            continue

        df = df.rename(columns={ko_col: "KO", prot_col: "Protein_ID"})[["KO", "Protein_ID"]].copy()
        df["KO"] = df["KO"].astype(str).str.strip().str.upper()
        df.loc[~df["KO"].str.startswith("K") & df["KO"].str.match(r"^\d{5}$"), "KO"] = "K" + df["KO"]
        df["sample"] = sample_name
        all_dfs.append(df)

    if not all_dfs:
        return pd.DataFrame(columns=["KO", "Protein_ID", "sample"])

    return pd.concat(all_dfs, ignore_index=True)

def _normalize_samples_to_faa_basename(df: pd.DataFrame, genomes_root: Path) -> pd.DataFrame:
    dir_to_faa = {}
    for d in genomes_root.iterdir():
        if d.is_dir():
            dir_to_faa[d.name] = _first_faa_basename(d) or d.name
    df = df.copy()
    df["sample"] = df["sample"].map(lambda s: dir_to_faa.get(str(s), str(s)))
    return df

def matrix_from_ko_map(kofam_df: pd.DataFrame, totals: pd.Series, pathway_map: dict, label_map: dict, order: list, genomes: list, value_mode: str = "percent") -> pd.DataFrame:
    rows = [label_map.get(p, p) for p in order]
    out = pd.DataFrame(index=rows, columns=genomes, data=0.0)

    for pathway in order:
        label = label_map.get(pathway, pathway)
        kos = pathway_map.get(pathway, [])
        if not kos:
            continue

        dfp = kofam_df[kofam_df["KO"].isin(kos)]
        if dfp.empty:
            continue

        gsum = dfp.groupby("sample")["Protein_ID"].nunique().astype(float)
        if value_mode == "percent":
            gN = totals.reindex(gsum.index).replace(0, np.nan)
            values = (gsum / gN * 100.0).replace([np.inf, -np.inf], np.nan).fillna(0.0)
        else:
            values = gsum

        for g, v in values.items():
            if g in out.columns:
                out.at[label, g] = float(v)
    return out.fillna(0.0)

def matrix_from_pfam(df: pd.DataFrame, mapping: dict, row_namer, genomes: list[str], value_mode: str = "percent") -> pd.DataFrame:
    agg = df.groupby(["sample", "hmm"], as_index=False).agg({"hits": "sum", "total_genes": "max"})
    rows = [row_namer(s) for s in mapping.keys()]
    out = pd.DataFrame(index=rows, columns=genomes, data=0.0)

    for sub, models in mapping.items():
        sub_df = agg[agg["hmm"].isin(models)]
        if sub_df.empty:
            continue
        gsum = sub_df.groupby("sample")["hits"].sum().astype(float)
        if value_mode == "percent":
            gN = sub_df.groupby("sample")["total_genes"].max().replace(0, np.nan)
            values = (gsum / gN * 100.0).replace([np.inf, -np.inf], np.nan).fillna(0.0)
        else:
            values = gsum
        row = row_namer(sub)
        for g, v in values.items():
            if g in out.columns:
                out.at[row, g] = float(v)
    return out.fillna(0.0)

def matrix_from_btex_no_map(df: pd.DataFrame, genomes: list[str], value_mode: str = "raw") -> pd.DataFrame:
    agg = df.groupby(["sample", "hmm"], as_index=False).agg({"hits": "sum", "total_genes": "max"})
    rows = list(dict.fromkeys(agg["hmm"].tolist()))
    out = pd.DataFrame(index=rows, columns=genomes, data=0.0)

    for sub in rows:
        sub_df = agg[agg["hmm"] == sub]
        if sub_df.empty:
            continue
        gsum = sub_df.groupby("sample")["hits"].sum().astype(float)
        if value_mode == "percent":
            gN = sub_df.groupby("sample")["total_genes"].max().replace(0, np.nan)
            values = (gsum / gN * 100.0).replace([np.inf, -np.inf], np.nan).fillna(0.0)
        else:
            values = gsum
        for g, v in values.items():
            if g in out.columns:
                out.at[sub, g] = float(v)
    return out.fillna(0.0)

def _figure_and_scale(n_rows, n_cols, cell_w_in=0.40, cell_h_in=0.36, max_fill=0.78, legend_width_in=3.2):
    fig_w = 1.6 + n_cols * cell_w_in + legend_width_in + 0.8
    fig_h = 1.2 + n_rows * cell_h_in + 1.2
    cell_w_pts = cell_w_in * 72.0
    cell_h_pts = cell_h_in * 72.0
    max_diam_pts = max_fill * min(cell_w_pts, cell_h_pts)
    max_radius_pts = max_diam_pts / 2.0
    return (fig_w, fig_h), max_radius_pts

def _areas_vmax_scaled(values, vmax, max_radius_pts, bubble_scale: float):
    v = np.asarray(values, dtype=float).ravel()
    v = np.clip(v, 0, None)
    if vmax is None or vmax <= 0 or v.size == 0:
        return np.zeros_like(v, dtype=float)
    radii = max_radius_pts * np.sqrt(v / float(vmax))
    areas = (radii ** 2) * np.pi * float(bubble_scale)
    areas[v == 0.0] = 0.0
    return areas

def _areas_linear(values, size_per_unit: float, base_area: float):
    v = np.asarray(values, dtype=float).ravel()
    v = np.clip(v, 0, None)
    areas = float(base_area) + float(size_per_unit) * v
    areas[v == 0.0] = 0.0
    return areas

def bubble_plot(matrix_vals: pd.DataFrame, genomes: list[str], features: list[str], out_svg: Path, value_mode: str, bubble_scale: float, sep_index: int | None = None, row_spacing: float = 2.6, ytick_pad: float = 10.0, title: str = "", legend_steps: list[float] | None = None, cell_w_in: float = 0.40, cell_h_in: float = 0.36, legend_width_in: float = 3.2, max_fill: float = 0.78, bubble_vmax: float | None = None, vmax_percentile: float | None = None, print_vmax: bool = False, bubble_size_per_hit: float = 0.0, bubble_base_area: float = 0.0):
    if matrix_vals.empty or len(features) == 0:
        print("Warning: Matrix is empty, skipping plot generation.")
        return

    M = matrix_vals.loc[features, genomes].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    ny, nx = len(features), len(genomes)
    vals = M.to_numpy(dtype=float)

    if bubble_size_per_hit and bubble_size_per_hit > 0:
        vmax = None
    else:
        if bubble_vmax is not None and bubble_vmax > 0:
            vmax = float(bubble_vmax)
        else:
            flat = vals.ravel()
            if vmax_percentile is not None and 0 < vmax_percentile <= 100 and flat.size > 0:
                vmax = float(np.percentile(flat, vmax_percentile))
                if vmax <= 0:
                    vmax = float(np.nanmax(flat)) if flat.size else 1.0
            else:
                vmax = float(np.nanmax(vals)) if vals.size else 1.0
            if not np.isfinite(vmax) or vmax <= 0:
                vmax = 1.0

    if print_vmax and vmax is not None:
        print(f"bubble_vmax_used={vmax}")

    (fig_w, fig_h), max_r = _figure_and_scale(n_rows=ny, n_cols=nx, cell_w_in=cell_w_in, cell_h_in=cell_h_in, max_fill=max_fill, legend_width_in=legend_width_in)

    X, Y, Svals = [], [], []
    y_pos = np.arange(ny, dtype=float) * float(row_spacing)
    for j, g in enumerate(genomes):
        col = M[g].to_numpy(dtype=float)
        for i, v in enumerate(col):
            X.append(j)
            Y.append(y_pos[i])
            Svals.append(v)

    Svals = np.asarray(Svals, dtype=float)

    if bubble_size_per_hit and bubble_size_per_hit > 0:
        A = _areas_linear(Svals, size_per_unit=float(bubble_size_per_hit), base_area=float(bubble_base_area))
    else:
        A = _areas_vmax_scaled(Svals, vmax, max_r, bubble_scale)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=160)
    ax.scatter(X, Y, s=A)
    ax.set_xticks(np.arange(nx))
    ax.set_xticklabels(genomes, rotation=90)
    ax.set_xlabel("Genomes")
    ymin = -0.5 * float(row_spacing)
    ymax = y_pos[-1] + 0.5 * float(row_spacing) if ny > 0 else 0.5
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(features)
    ax.tick_params(axis="y", pad=ytick_pad)
    ax.invert_yaxis()
    ax.set_ylabel("Features")
    ax.set_title(title, fontsize=12)
    ax.grid(True, axis="both", linestyle=":", linewidth=0.5)
    if sep_index is not None and 0 < sep_index < nx:
        ax.axvline(x=sep_index - 0.5, linewidth=2)

    if legend_steps:
        legend_vals = [float(v) for v in legend_steps if float(v) > 0]
    else:
        # automatic legend based on data
        # use vmax (for vmax-scaling mode) or matrix max (for linear mode)
        if bubble_size_per_hit and bubble_size_per_hit > 0:
            raw_max = np.nanmax(vals) if vals.size else 1.0
        else:
            raw_max = vmax if "vmax" in locals() and vmax is not None else (
                np.nanmax(vals) if vals.size else 1.0
            )

        if not np.isfinite(raw_max) or raw_max <= 0:
            raw_max = 1.0

        # make the top legend value slightly larger than the max
        # and round to a nice scale
        import math
        order = math.floor(math.log10(raw_max)) if raw_max > 0 else 0.0
        scale = 10.0 ** order
        top = math.ceil((raw_max / scale) * 1.05) * scale

        # three legend values: 1/4, 1/2, and full top
        legend_vals = [top / 4.0, top / 2.0, top]

    handles, labels = [], []
    for v in legend_vals:
        if bubble_size_per_hit and bubble_size_per_hit > 0:
            a = _areas_linear([v], size_per_unit=float(bubble_size_per_hit), base_area=float(bubble_base_area))[0]
        else:
            a = _areas_vmax_scaled([v], vmax, max_r, bubble_scale)[0]
        h = ax.scatter([], [], s=a, edgecolors="black", facecolors="none")
        labels.append(f"{v:g}" if value_mode != "percent" else f"{v:.1f}%")
        handles.append(h)

    if handles:
        legend_title = "Bubble = hit count" if value_mode != "percent" else "Bubble = % of total genes"
        ax.legend(handles, labels, title=legend_title, loc="center left",
                  bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, frameon=True)

    pad_cells = 0.02
    ax.set_xlim(-0.5 - pad_cells, nx - 0.5 + pad_cells)
    fig.savefig(out_svg, dpi=300, bbox_inches="tight")
    plt.close(fig)

def aggregate_rows_by_groups(mat: pd.DataFrame, row_to_group: dict, group_order: list[str]) -> pd.DataFrame:
    if mat.empty or not row_to_group:
        return mat

    group2rows = {}
    for r in mat.index:
        g = row_to_group.get(r)
        if g is None:
            continue
        group2rows.setdefault(g, []).append(r)

    frames = []
    for g in group_order:
        rows = group2rows.get(g)
        if rows:
            frames.append(pd.DataFrame([mat.loc[rows].sum()], index=[g]))

    if not frames:
        return pd.DataFrame(index=[], columns=mat.columns).fillna(0.0)
    out = pd.concat(frames, axis=0)
    return out.reindex(columns=mat.columns).fillna(0.0)

def _read_order_file(p: Path) -> list[str]:
    lines = []
    with p.open() as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            lines.append(s)
    return lines

def main():
    ap = argparse.ArgumentParser(description="Bubble plot of BTEX PFAM KO hits without BTEX custom map, with optional BTEX pathway aggregation; zero values draw no bubble")
    ap.add_argument("--btex-hmm", required=True, help="BTEX hmmsearch summary CSV (sample,hmm,hits,total_genes)")
    ap.add_argument("--pfam-hmm", required=True, help="PFAM hmmsearch summary CSV (sample,hmm,hits,total_genes)")
    ap.add_argument("--kofam-hits-dir", required=True, help="Directory with per-genome subdirs containing kofam_btex_hits.tsv")
    ap.add_argument("--btex-pathway-map", default=None, help="BTEX pathway map TSV (subunit, pathway) to aggregate BTEX subunits to pathway rows")
    ap.add_argument("--include-pfam-ko", action="store_true", help="Include PFAM and KO rows in the plot even when BTEX pathway map is provided")
    ap.add_argument("--pfam-map", required=True, help="PFAM map TSV (subunit, pfam_ids, rule[, pfam_names])")
    ap.add_argument("--ko-map", required=True, help="KO map TSV (pathway, ko_ids, display_label)")
    ap.add_argument("--out-svg", required=True, help="Output image path")
    ap.add_argument("--out-tsv", default=None, help="Optional matrix TSV used for plotting")
    ap.add_argument("--out-order", default=None, help="Optional path to write plotted sample order; defaults to <out-svg>.samples_order.txt")
    ap.add_argument("--pathway-filter", default=None, help="Optional KO pathway filter by pathway name in KO map")
    ap.add_argument("--name-joiner", default="+", help="Joiner for pfam_names if multiple per subunit")
    ap.add_argument("--row-spacing", type=float, default=2.6)
    ap.add_argument("--ytick-pad", type=float, default=10.0)
    ap.add_argument("--legend-steps", type=str, default="")
    ap.add_argument("--cell-width", type=float, default=0.40)
    ap.add_argument("--cell-height", type=float, default=0.36)
    ap.add_argument("--legend-width", type=float, default=3.2)
    ap.add_argument("--max-fill", type=float, default=0.78)
    ap.add_argument("--bubble-scale", type=float, default=1.4, help="Multiplier for bubble areas in vmax scaling mode")
    ap.add_argument("--value-mode", choices=("percent","raw"), default="raw", help="Percent of total genes or raw hit counts")
    ap.add_argument("--bubble-vmax", type=float, default=None, help="Reference value for full size bubble in vmax scaling mode")
    ap.add_argument("--vmax-percentile", type=float, default=None, help="Percentile to compute vmax when bubble-vmax not given")
    ap.add_argument("--print-vmax", action="store_true", help="Print the vmax used in vmax scaling mode")
    ap.add_argument("--bubble-size-per-hit", type=float, default=0.0, help="Area points^2 added per unit value for linear bubble sizing")
    ap.add_argument("--bubble-base-area", type=float, default=0.0, help="Base area points^2 when value is zero in linear mode; zero ensures no bubble for zero")
    args = ap.parse_args()

    btex = read_hmm_csv(Path(args.btex_hmm), is_pfam=False)
    pfam = read_hmm_csv(Path(args.pfam_hmm), is_pfam=True)

    genomes_root = Path(args.kofam_hits_dir)
    btex = _normalize_samples_to_faa_basename(btex, genomes_root)
    pfam = _normalize_samples_to_faa_basename(pfam, genomes_root)

    totals_b = (
        btex[["sample","total_genes"]]
        .drop_duplicates("sample").set_index("sample")["total_genes"]
        .astype(float)
    )
    totals_p = (
        pfam[["sample","total_genes"]]
        .drop_duplicates("sample").set_index("sample")["total_genes"]
        .astype(float)
    )
    totals = totals_b.combine_first(totals_p)

    pfam_order, sub2pfams, sub2label = load_pfam_map(Path(args.pfam_map), name_joiner=args.name_joiner)
    ko_order, pathway2kos, pathway2label = load_ko_map(Path(args.ko_map))

    sub2path = {}
    path_order = []
    if args.btex_pathway_map:
        sub2path, path_order = load_btex_pathway_map(Path(args.btex_pathway_map))

    ko_hits_df = read_kofam_hits(Path(args.kofam_hits_dir))

    btex_samples = set(btex["sample"].unique())
    pfam_samples = set(pfam["sample"].unique())
    ko_samples = set(ko_hits_df["sample"].unique())

    btex_pfam_union = btex_samples | pfam_samples

    # genomes in BTEX or PFAM but with no KO hits at all
    missing_ko = sorted(btex_pfam_union - ko_samples)

    # genomes that have KO hits but do not appear in BTEX or PFAM
    extra_ko = sorted(ko_samples - btex_pfam_union)

    samples_ref = sorted(btex_pfam_union | ko_samples)
    ko_hits_df = ko_hits_df[ko_hits_df["sample"].isin(samples_ref)]

    print("BTEX samples:", len(btex_samples))
    print("PFAM samples:", len(pfam_samples))
    print("KO unique samples before filter:", len(ko_samples))
    print("Example KO samples:", sorted(ko_samples)[:5])
    print("Totals index size:", totals.index.size)
    print("Genomes with BTEX/PFAM but no KO hits:", missing_ko)
    print("Genomes with KO hits but no BTEX/PFAM:", extra_ko)

    if args.value_mode == "percent":
        missing_totals = [s for s in samples_ref if s not in totals.index]
        if missing_totals:
            print("Warning: samples missing totals for percent mode:", missing_totals[:5], "...")

    genomes = samples_ref

    btex_vals = matrix_from_btex_no_map(btex, genomes=genomes, value_mode=args.value_mode)
    pfam_vals = matrix_from_pfam(
        pfam, sub2pfams,
        row_namer=lambda s: sub2label.get(s, s),
        genomes=genomes, value_mode=args.value_mode
    )
    ko_vals = matrix_from_ko_map(
        ko_hits_df, totals, pathway2kos, pathway2label,
        ko_order, genomes, value_mode=args.value_mode
    )

    plot_title = "BTEX PFAM KO Hits as percent of total genes" if args.value_mode == "percent" else "BTEX PFAM KO Hits raw counts"

    if args.pathway_filter:
        pkeep = args.pathway_filter
        if pkeep not in pathway2label:
            raise SystemExit(f"Error: Pathway '{pkeep}' not found in KO map first column")
        lkeep = pathway2label[pkeep]
        plot_title = f"BTEX PFAM and {pkeep} Hits"
        if args.value_mode == "percent":
            plot_title += " as percent of total genes"
        ko_vals = ko_vals.loc[[lkeep]] if lkeep in ko_vals.index else pd.DataFrame()

    btex_vals = btex_vals.groupby(level=0).sum().T.groupby(level=0).sum().T
    pfam_vals = pfam_vals.groupby(level=0).sum().T.groupby(level=0).sum().T

    if sub2path:
        mapped = set(sub2path.keys())
        keep_rows = [row for row in btex_vals.index if row in mapped]
        btex_vals = btex_vals.loc[keep_rows] if keep_rows else pd.DataFrame(index=[], columns=btex_vals.columns).fillna(0.0)
        row_to_group = {row: sub2path[row] for row in keep_rows}
        btex_vals = aggregate_rows_by_groups(btex_vals, row_to_group=row_to_group, group_order=path_order)
        btex_rows = [p for p in path_order if p in btex_vals.index]
        if args.include_pfam_ko:
            pfam_rows = [sub2label.get(s, s) for s in pfam_order if sub2label.get(s, s) in pfam_vals.index]
            ko_rows = [pathway2label.get(p, p) for p in ko_order if pathway2label.get(p, p) in ko_vals.index]
            matrix_vals = pd.concat([btex_vals, pfam_vals, ko_vals], axis=0).fillna(0.0)
            ordered_pref = btex_rows + pfam_rows + ko_rows
        else:
            matrix_vals = btex_vals.copy()
            ordered_pref = btex_rows
    else:
        btex_rows = list(btex_vals.index)
        pfam_rows = [sub2label.get(s, s) for s in pfam_order if sub2label.get(s, s) in pfam_vals.index]
        ko_rows = [pathway2label.get(p, p) for p in ko_order if pathway2label.get(p, p) in ko_vals.index]
        matrix_vals = pd.concat([btex_vals, pfam_vals, ko_vals], axis=0).fillna(0.0)
        ordered_pref = btex_rows + pfam_rows + ko_rows

    extras = [r for r in matrix_vals.index if r not in ordered_pref]
    features = list(dict.fromkeys(ordered_pref + extras))

    matrix_vals = matrix_vals.groupby(level=0).sum().reindex(index=[f for f in features if f in matrix_vals.index])
    matrix_vals = matrix_vals.T.groupby(level=0).sum().T
    # determine default genomes from columns present
    default_genomes = list(matrix_vals.columns)
    matrix_vals = matrix_vals.reindex(columns=default_genomes, fill_value=0.0)
    matrix_vals = matrix_vals.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # sample order file logic
    order_path = Path(args.out_order) if args.out_order else Path(args.out_svg).with_suffix(".samples_order.txt")
    order_path_parent = order_path.parent
    order_path_parent.mkdir(parents=True, exist_ok=True)

    if order_path.exists():
        desired = _read_order_file(order_path)
        # keep only those present, then append any missing at end in their current order
        genomes = [g for g in desired if g in matrix_vals.columns] + [g for g in matrix_vals.columns if g not in desired]
        print(f"using existing sample order: {order_path}")
        # do not overwrite the file
    else:
        genomes = default_genomes
        with open(order_path, "w") as f:
            for s in genomes:
                f.write(s + "\n")
        print(f"wrote {order_path}")

    # reorder columns to match chosen genomes order
    matrix_vals = matrix_vals.reindex(columns=genomes)

    legend_steps = [float(x) for x in args.legend_steps.split(",")] if args.legend_steps.strip() else None

    bubble_plot(
        matrix_vals=matrix_vals,
        genomes=genomes,
        features=list(matrix_vals.index),
        out_svg=Path(args.out_svg),
        value_mode=args.value_mode,
        bubble_scale=max(0.1, float(args.bubble_scale)),
        sep_index=None,
        row_spacing=float(args.row_spacing),
        ytick_pad=float(args.ytick_pad),
        title=plot_title,
        legend_steps=legend_steps,
        cell_w_in=float(args.cell_width),
        cell_h_in=float(args.cell_height),
        legend_width_in=float(args.legend_width),
        max_fill=float(args.max_fill),
        bubble_vmax=args.bubble_vmax,
        vmax_percentile=args.vmax_percentile,
        print_vmax=args.print_vmax,
        bubble_size_per_hit=float(args.bubble_size_per_hit),
        bubble_base_area=float(args.bubble_base_area),
    )

    if args.out_tsv:
        matrix_vals.to_csv(args.out_tsv, sep="\t")

if __name__ == "__main__":
    main()
