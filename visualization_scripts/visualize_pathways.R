#!/usr/bin/env Rscript

# plot_kegg_pathways_links_only.R
#
# Purpose
# Generate KEGG show_pathway URLs that color KOs from hmmscan hits.
# This version DOES NOT download PNGs or do any browser/screenshot capture.
#
# Usage
# Rscript /home/juneq/hmm/scripts/kegg/plot_kegg_pathways_2.R \
#   --hmmscan /home/juneq/BTEX_test/test_output_2/btex_hmm_summary.csv \
#   --sample all \
#   --ko-map /home/juneq/hmm/archetypes/hmm_cutoffs/btex_ko_list.tsv \
#   --pathways 00642,00623,00622,00362 \
#   --outdir /home/juneq/BTEX_test/test_output_2/pathway_figures

# Rscript /home/juneq/hmm/scripts/kegg/plot_kegg_pathways_2.R \
#   --hmmscan /home/juneq/hmm/validation_genomes/atested/abtested_results/hmmscan_summary_all_hits_260216_BTEX.csv \
#   --sample burkholderia_vietnamiensis_g4 \
#   --ko-map /home/juneq/hmm/archetypes/hmm_cutoffs/btex_ko_list.tsv \
#   --pathways 00623,00642,00622,00362 \
#   --outdir /home/juneq/hmm/validation_genomes/atested/genomes/burkholderia_vietnamiensis_g4__GCA_000016205.1/pathview_outputs


# Outputs in --outdir
# KEGG_MAP_LINKS.txt
# sample_color_legend.tsv
# split_color_inputs.<pid>.tsv

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--hmmscan"), type = "character",
              help = "hmmscan CSV with at least sample,hmm columns"),
  make_option(c("--sample"), type = "character",
              help = "Sample name to plot. Use 'all' for split plots across all samples."),
  make_option(c("--ko-map"), type = "character", dest = "ko_map",
              default = "ko_map.tsv",
              help = "KO mapping table (gene/hmm/subunit + ko column)"),
  make_option(c("--pathways"), type = "character",
              help = "Comma separated KEGG map IDs, e.g. 00623,00642"),
  make_option(c("--outdir"), type = "character", default = "kegg_urls",
              help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

set_kegg_ua <- function() {
  options(url.method = "libcurl")
  options(download.file.method = "libcurl")
  options("HTTPUserAgent" = "R-KEGG-Client/links-only")
}

read_ko_map <- function(path) {
  tab <- read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  nms_l <- tolower(names(tab))

  hmm_idx <- which(nms_l %in% c("gene", "hmm", "subunit"))[1]
  ko_idx  <- which(nms_l %in% c("ko", "ko_id"))[1]
  if (is.na(hmm_idx) || is.na(ko_idx)) {
    stop("KO map must have a column named gene/hmm/subunit and a column named ko/ko_id")
  }

  hmm_col <- names(tab)[hmm_idx]
  ko_col  <- names(tab)[ko_idx]

  out <- list()
  for (i in seq_len(nrow(tab))) {
    key <- tolower(trimws(tab[i, hmm_col]))
    kos <- trimws(unlist(strsplit(as.character(tab[i, ko_col]), ",")))
    kos <- kos[grepl("^K[0-9]{5}$", kos)]
    if (length(kos) > 0) out[[key]] <- unique(c(out[[key]], kos))
  }
  out
}

get_pathway_kos_from_kgml <- function(pid) {
  kgml_url <- paste0("https://rest.kegg.jp/get/ko", pid, "/kgml")
  kgml_lines <- tryCatch(readLines(kgml_url, warn = FALSE), error = function(e) character(0))
  kgml_txt <- paste(kgml_lines, collapse = "\n")
  if (!nzchar(kgml_txt)) return(character(0))

  entries <- regmatches(kgml_txt, gregexpr("<entry\\b[^>]*>", kgml_txt, perl = TRUE))[[1]]
  if (length(entries) == 0) return(character(0))
  ortholog_entries <- entries[grepl('type="ortholog"', ortholog_entries <- entries, fixed = TRUE)]
  if (length(ortholog_entries) == 0) return(character(0))

  name_matches <- regmatches(ortholog_entries, regexec('name="([^"]+)"', ortholog_entries, perl = TRUE))
  name_attrs <- vapply(name_matches, function(m) if (length(m) >= 2) m[2] else "", character(1))
  name_attrs <- name_attrs[nzchar(name_attrs)]
  if (length(name_attrs) == 0) return(character(0))

  toks <- unlist(strsplit(paste(name_attrs, collapse = " "), "\\s+"), use.names = FALSE)
  kos <- sub("^ko:", "", toks)
  kos <- kos[grepl("^K[0-9]{5}$", kos)]
  unique(kos)
}

get_pathway_ko_groups_from_kgml <- function(pid) {
  kgml_url <- paste0("https://rest.kegg.jp/get/ko", pid, "/kgml")
  kgml_lines <- tryCatch(readLines(kgml_url, warn = FALSE), error = function(e) character(0))
  kgml_txt <- paste(kgml_lines, collapse = "\n")
  if (!nzchar(kgml_txt)) return(list())

  entries <- regmatches(kgml_txt, gregexpr("<entry\\b[^>]*>", kgml_txt, perl = TRUE))[[1]]
  if (length(entries) == 0) return(list())
  ortholog_entries <- entries[grepl('type="ortholog"', entries, fixed = TRUE)]
  if (length(ortholog_entries) == 0) return(list())

  name_matches <- regmatches(ortholog_entries, regexec('name="([^"]+)"', ortholog_entries, perl = TRUE))
  name_attrs <- vapply(name_matches, function(m) if (length(m) >= 2) m[2] else "", character(1))
  name_attrs <- name_attrs[nzchar(name_attrs)]
  if (length(name_attrs) == 0) return(list())

  groups <- lapply(name_attrs, function(x) {
    toks <- unlist(strsplit(x, "\\s+"), use.names = FALSE)
    kos <- sub("^ko:", "", toks)
    unique(kos[grepl("^K[0-9]{5}$", kos)])
  })
  groups <- groups[vapply(groups, length, integer(1)) > 0]
  groups
}

get_white_kos_for_non_btex_nodes <- function(pathway_ko_groups, ko_universe) {
  if (length(pathway_ko_groups) == 0) return(character(0))
  ko_universe <- unique(ko_universe[grepl("^K[0-9]{5}$", ko_universe)])

  white_groups <- pathway_ko_groups[vapply(pathway_ko_groups, function(g) {
    length(intersect(g, ko_universe)) == 0
  }, logical(1))]

  if (length(white_groups) == 0) return(character(0))
  unique(unlist(white_groups, use.names = FALSE))
}

get_sample_colors <- function(samples) {
  palette <- c(
    "#E41A1C","#377EB8","#FF7F00","#4DAF4A","#984EA3",
    "#00BCD4","#F781BF","#A65628","#FFD92F","#FB9A99",
    "#999999","#66C2A5","#1B9E77","#D95F02","#7570B3",
    "#E7298A","#66A61E","#E6AB02","#A6761D","#1F78B4"
  )
  if (length(samples) > length(palette)) {
    stop("Too many samples for the hardcoded palette. Extend palette in get_sample_colors().")
  }
  setNames(palette[seq_along(samples)], samples)
}

urlenc_spaces_newlines_only <- function(x) {
  # Keep KEGG token syntax (ko:Kxxxxx) intact; only encode characters needed in query payload.
  x <- gsub("%", "%25", x, fixed = TRUE)
  x <- gsub("#", "%23", x, fixed = TRUE)
  x <- gsub("\r\n", "\n", x, fixed = TRUE)
  x <- gsub("\n", "%0A", x, fixed = TRUE)
  x <- gsub(" ", "%20", x, fixed = TRUE)
  x
}

build_multi_query_btex_outline_split <- function(pathway_kos, ko_universe, kos_by_sample, sample_colors, white_kos = character(0)) {
  samples <- names(sample_colors)
  pathway_kos <- unique(pathway_kos[grepl("^K[0-9]{5}$", pathway_kos)])
  ko_universe <- unique(ko_universe[grepl("^K[0-9]{5}$", ko_universe)])
  target_kos <- sort(intersect(pathway_kos, ko_universe))
  white_kos <- sort(unique(white_kos[grepl("^K[0-9]{5}$", white_kos)]))
  if (length(target_kos) == 0 && length(white_kos) == 0) return("")

  lines <- character(0)
  if (length(white_kos) > 0) {
    lines <- c(lines, vapply(white_kos, function(ko) {
      paste0("ko:", ko, " #FFFFFF,#000000")
    }, character(1)))
  }

  for (ko in target_kos) {
    present_samples <- samples[vapply(samples, function(s) ko %in% kos_by_sample[[s]], logical(1))]
    if (length(present_samples) > 0) {
      cols <- paste0(unname(sample_colors[present_samples]), ",#FF0000")
      lines <- c(lines, paste(c(paste0("ko:", ko), cols), collapse = " "))
    } else {
      lines <- c(lines, paste0("ko:", ko, " #FFFFFF,#FF0000"))
    }
  }

  urlenc_spaces_newlines_only(paste(lines, collapse = "\n"))
}

build_multi_query_btex_outline_single <- function(pathway_kos, ko_universe, sample_kos, sample_color, white_kos = character(0)) {
  pathway_kos <- unique(pathway_kos[grepl("^K[0-9]{5}$", pathway_kos)])
  ko_universe <- unique(ko_universe[grepl("^K[0-9]{5}$", ko_universe)])
  target_kos <- sort(intersect(pathway_kos, ko_universe))
  white_kos <- sort(unique(white_kos[grepl("^K[0-9]{5}$", white_kos)]))
  if (length(target_kos) == 0 && length(white_kos) == 0) return("")

  sample_kos <- unique(sample_kos[grepl("^K[0-9]{5}$", sample_kos)])
  sample_color_red <- paste0(sample_color, ",#FF0000")
  lines <- character(0)

  if (length(white_kos) > 0) {
    lines <- c(lines, vapply(white_kos, function(ko) {
      paste0("ko:", ko, " #FFFFFF,#000000")
    }, character(1)))
  }

  for (ko in target_kos) {
    if (ko %in% sample_kos) {
      lines <- c(lines, paste0("ko:", ko, " ", sample_color_red))
    } else {
      lines <- c(lines, paste0("ko:", ko, " #FFFFFF,#FF0000"))
    }
  }

  urlenc_spaces_newlines_only(paste(lines, collapse = "\n"))
}

write_sample_color_legend <- function(path, sample_colors) {
  df <- data.frame(sample = names(sample_colors),
                   color = unname(sample_colors),
                   stringsAsFactors = FALSE)
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

write_split_color_inputs <- function(path, pathway_kos, kos_by_sample, sample_colors) {
  samples <- names(sample_colors)
  mat <- sapply(samples, function(s) pathway_kos %in% kos_by_sample[[s]])
  storage.mode(mat) <- "integer"

  present_samples <- apply(mat, 1, function(row) {
    ss <- samples[as.logical(row)]
    if (length(ss) == 0) "" else paste(ss, collapse = ",")
  })

  present_colors <- apply(mat, 1, function(row) {
    ss <- samples[as.logical(row)]
    if (length(ss) == 0) "" else paste(unname(sample_colors[ss]), collapse = ",")
  })

  out <- data.frame(
    ko = pathway_kos,
    present_samples = present_samples,
    present_colors = present_colors,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(samples)) out[[samples[i]]] <- mat[, i]

  write.table(out, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

set_kegg_ua()

hm <- read.csv(opt$hmmscan, stringsAsFactors = FALSE, check.names = FALSE)
ko_map <- read_ko_map(opt$ko_map)
ko_universe <- sort(unique(unlist(ko_map, use.names = FALSE)))
ko_universe <- ko_universe[grepl("^K[0-9]{5}$", ko_universe)]

sample_col_idx <- grep("sample", names(hm), ignore.case = TRUE)[1]
hmm_col_idx    <- grep("hmm", names(hm), ignore.case = TRUE)[1]
hits_col_idx   <- grep("^hits$|hit_count|num_hits|n_hits", names(hm), ignore.case = TRUE)[1]
if (is.na(sample_col_idx) || is.na(hmm_col_idx) || is.na(hits_col_idx)) {
  stop("hmmscan CSV must have columns that match /sample/i, /hmm/i, and /hits|hit_count|num_hits|n_hits/i")
}
sample_col <- names(hm)[sample_col_idx]
hmm_col    <- names(hm)[hmm_col_idx]
hits_col   <- names(hm)[hits_col_idx]

kos_by_sample <- lapply(split(hm, hm[[sample_col]]), function(df) {
  # Only map KOs from HMMs with observed hits in this sample.
  hit_vals <- suppressWarnings(as.numeric(df[[hits_col]]))
  df_hits <- df[!is.na(hit_vals) & hit_vals > 0, , drop = FALSE]

  kos <- unlist(lapply(df_hits[[hmm_col]], function(h) {
    key <- tolower(trimws(h))
    ko_map[[key]]
  }), use.names = FALSE)
  kos <- unique(kos)
  kos[grepl("^K[0-9]{5}$", kos)]
})

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(opt$outdir, "KEGG_MAP_LINKS.txt")

pids <- trimws(strsplit(opt$pathways, ",")[[1]])
cat("Generated KEGG show_pathway Links\n=================================\n", file = log_file)

if (!is.null(opt$sample) && tolower(opt$sample) != "all") {
  samp <- opt$sample
  if (!samp %in% names(kos_by_sample)) stop("Requested --sample not found in hmmscan CSV: ", samp)
  run_samples <- c(samp)
} else {
  all_samples <- sort(names(kos_by_sample))
  n_detected <- length(all_samples)
  message("Detected ", n_detected, " samples in hmmscan input.")

  if (n_detected > 20) {
    run_samples <- all_samples[seq_len(20)]
    message("Using top 20 samples (alphabetical) for downstream processing due to palette limit.")
  } else {
    run_samples <- all_samples
  }
}
message("Generating pathway visualizations on KEGG", paste(pids, collapse = ","), " maps.")

sample_colors <- get_sample_colors(run_samples)
if (length(run_samples) == 1) {
  sample_colors <- setNames("#377EB8", run_samples)
}
legend_path <- file.path(opt$outdir, "sample_color_legend.tsv")
write_sample_color_legend(legend_path, sample_colors)

for (pid in pids) {
  kos_by_sample_run <- kos_by_sample[run_samples]
  pathway_kos <- get_pathway_kos_from_kgml(pid)
  target_kos <- if (length(pathway_kos) > 0) pathway_kos else ko_universe

  if (!is.null(opt$sample) && tolower(opt$sample) != "all") {
    mq <- build_multi_query_btex_outline_single(
      pathway_kos = target_kos,
      ko_universe = ko_universe,
      sample_kos = kos_by_sample[[run_samples[1]]],
      sample_color = sample_colors[[run_samples[1]]]
    )

    if (!nzchar(mq)) {
      url <- paste0("https://www.kegg.jp/kegg-bin/show_pathway?map=map", pid)
      cat(paste0("\nPathway ", pid, " (single sample: ", run_samples[1], ", no BTEX KOs on pathway):\n", url, "\n"),
          file = log_file, append = TRUE)
      next
    }

    show_url <- paste0("https://www.kegg.jp/kegg-bin/show_pathway?map=map", pid, "&multi_query=", mq)

    cat(paste0("\nPathway ", pid, " (single sample: ", run_samples[1], ", BTEX outlined + split hits):\n", show_url, "\n"),
        file = log_file, append = TRUE)

  } else {
    mq <- build_multi_query_btex_outline_split(
      pathway_kos = target_kos,
      ko_universe = ko_universe,
      kos_by_sample = kos_by_sample_run,
      sample_colors = sample_colors
    )

    if (!nzchar(mq)) {
      url <- paste0("https://www.kegg.jp/kegg-bin/show_pathway?map=map", pid)
      cat(paste0("\nPathway ", pid, " (split across samples, no BTEX KOs on pathway):\n", url, "\n"),
          file = log_file, append = TRUE)
      next
    }

    show_url <- paste0("https://www.kegg.jp/kegg-bin/show_pathway?map=map", pid, "&multi_query=", mq)

    cat(paste0("\nPathway ", pid, " (split across samples, BTEX outlined + split hits):\n", show_url, "\n"),
        file = log_file, append = TRUE)
  }

  if (length(pathway_kos) > 0) {
    split_inputs_path <- file.path(opt$outdir, paste0("split_color_inputs.", pid, ".tsv"))
    write_split_color_inputs(split_inputs_path, pathway_kos, kos_by_sample_run, sample_colors)
  }
}

message("\nDone!")
message("Links: ", log_file)
message("Legend: ", legend_path)
message("Presence tables: split_color_inputs.<pid>.tsv in ", opt$outdir)
