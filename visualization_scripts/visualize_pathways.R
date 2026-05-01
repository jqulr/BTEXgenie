#!/usr/bin/env Rscript

#
# This script writes HTML launchers and tsv files for inspection, it does not download images.
#
# Output directory contains:
# {pathway_name}_{pathway_ID}.html
# sample_color_legend.tsv

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--hmmscan"), type = "character",
              help = "hmmscan CSV with at least sample,hmm columns"),
  make_option(c("--genome-dir"), type = "character", dest = "genome_dir",
              help = "Optional genome directory with each genome subdirectory containing kofam_abv_thres.tsv"),
  make_option(c("--sample"), type = "character",
              help = "Sample name to plot. Use 'all' for split plots across all samples."),
  make_option(c("--ko-map"), type = "character", dest = "ko_map",
              default = "ko_map.tsv",
              help = "KO mapping table (gene/hmm/subunit + ko column)"),
  make_option(c("--pathways"), type = "character",
              help = "Comma separated KEGG map IDs, e.g. 00623,00642"),
  make_option(c("--pathways-supplied"), action = "store_true", dest = "pathways_supplied", default = FALSE,
              help = "Internal flag indicating --pathways was explicitly supplied"),
  make_option(c("--outdir"), type = "character", default = "kegg_urls",
              help = "Output directory"),
  make_option(c("--fg-color"), type = "character", dest = "fg_color",
              default = "#FF0000",
              help = "Outline color for all KOs on the pathway (default: #FF0000 red). Use 'none' to disable outlines.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if ((is.null(opt$hmmscan) || !nzchar(opt$hmmscan)) && (is.null(opt$genome_dir) || !nzchar(opt$genome_dir))) {
  stop("Provide either --hmmscan or --genome-dir")
}
if (!is.null(opt$hmmscan) && nzchar(opt$hmmscan) && !is.null(opt$genome_dir) && nzchar(opt$genome_dir)) {
  stop("Use only one input mode: --hmmscan or --genome-dir")
}

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

sanitize_filename_part <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

normalize_pathway_id <- function(pid) {
  pid <- trimws(pid)
  pid <- sub("^(map|ko)", "", pid, ignore.case = TRUE)
  pid
}

pathway_metadata_cache <- new.env(parent = emptyenv())

get_pathway_metadata <- function(pid) {
  pid <- normalize_pathway_id(pid)
  if (exists(pid, envir = pathway_metadata_cache, inherits = FALSE)) {
    return(get(pid, envir = pathway_metadata_cache, inherits = FALSE))
  }

  entry_id <- paste0("map", pid)
  entry_url <- paste0("https://rest.kegg.jp/get/", entry_id)
  entry_lines <- tryCatch(readLines(entry_url, warn = FALSE), error = function(e) character(0))

  title <- ""
  if (length(entry_lines) > 0) {
    name_idx <- grep("^NAME\\s+", entry_lines)
    if (length(name_idx) > 0) {
      title <- trimws(sub("^NAME\\s+", "", entry_lines[name_idx[1]]))
      title <- sub("\\s*-\\s*reference pathway$", "", title, ignore.case = TRUE)
    }
  }

  if (!nzchar(title)) {
    title <- paste("pathway", pid)
  }

  filename_stub <- sanitize_filename_part(title)
  if (!nzchar(filename_stub)) {
    filename_stub <- paste0("pathway_", pid)
  }

  out <- list(pid = pid, title = title, filename_stub = filename_stub)
  assign(pid, out, envir = pathway_metadata_cache)
  out
}

read_kos_from_kofam_detail <- function(path) {
  tab <- read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  ko_idx <- which(tolower(names(tab)) %in% c("ko", "ko_id"))[1]
  if (is.na(ko_idx)) {
    stop("kofam_detail_abv_thres.tsv must contain a KO or ko_id column")
  }
  kos <- trimws(as.character(tab[[ko_idx]]))
  unique(kos[grepl("^K[0-9]{5}$", kos)])
}

read_kos_by_genome_from_dir <- function(genome_dir) {
  if (!dir.exists(genome_dir)) {
    stop("genome directory not found: ", genome_dir)
  }

  subdirs <- list.dirs(genome_dir, recursive = FALSE, full.names = TRUE)
  if (length(subdirs) == 0) {
    stop("No genome subdirectories found in genome directory: ", genome_dir)
  }

  out <- list()
  for (subdir in sort(subdirs)) {
    detail_path <- file.path(subdir, "kofam_abv_thres.tsv")
    if (!file.exists(detail_path)) next
    genome_name <- basename(subdir)
    out[[genome_name]] <- read_kos_from_kofam_detail(detail_path)
  }

  if (length(out) == 0) {
    stop("No kofam_abv_thres.tsv files found under genome directory: ", genome_dir)
  }
  if (length(out) > 7) {
    warning(
      "Detected ", length(out), " genomes under --genome-dir. More than 7 genomes may not display properly on the KEGG pathway split-color view.",
      call. = FALSE
    )
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

get_sample_colors <- function(samples) {
  palette <- c(
    "#F08A8B","#377EB8","#FF7F00","#4DAF4A","#984EA3",
    "#00BCD4","#F781BF","#A65628","#FFD92F","#a4df0f",
    "#999999","#66C2A5","#1B9E77","#D95F02","#7570B3",
    "#E7298A","#66A61E","#E6AB02","#A6761D","#1F78B4"
  )
  if (length(samples) > length(palette)) {
    stop("Too many samples for the hardcoded palette. Extend palette in get_sample_colors().")
  }
  setNames(palette[seq_along(samples)], samples)
}

urlenc_spaces_newlines_only <- function(x) {
  x <- gsub("%", "%25", x, fixed = TRUE)
  x <- gsub("#", "%23", x, fixed = TRUE)
  x <- gsub("\r\n", "\n", x, fixed = TRUE)
  x <- gsub("\n", "%0A", x, fixed = TRUE)
  x <- gsub(" ", "%20", x, fixed = TRUE)
  x
}

build_multi_query_btex_outline_single <- function(pathway_kos, ko_universe, sample_kos, sample_color, pathway_ko_groups = list(), fg_color = "#FF0000") {
  pathway_kos <- unique(pathway_kos[grepl("^K[0-9]{5}$", pathway_kos)])
  ko_universe <- unique(ko_universe[grepl("^K[0-9]{5}$", ko_universe)])
  target_kos <- sort(intersect(pathway_kos, ko_universe))
  if (length(target_kos) == 0) return("")

  sample_kos <- unique(sample_kos[grepl("^K[0-9]{5}$", sample_kos)])
  use_outline <- !is.null(fg_color) && nzchar(fg_color) && tolower(fg_color) != "none"

  lines <- character(0)
  group_keys <- character(0)
  target_groups <- list()

  if (length(pathway_ko_groups) > 0) {
    for (g in pathway_ko_groups) {
      gg <- sort(unique(intersect(g, target_kos)))
      if (length(gg) == 0) next
      k <- paste(gg, collapse = "|")
      if (k %in% group_keys) next
      group_keys <- c(group_keys, k)
      target_groups <- c(target_groups, list(gg))
    }
  }

  covered <- if (length(target_groups) > 0) unique(unlist(target_groups, use.names = FALSE)) else character(0)
  for (ko in setdiff(target_kos, covered)) {
    target_groups <- c(target_groups, list(ko))
  }

  for (grp in target_groups) {
    if (!any(grp %in% sample_kos)) {
      if (use_outline) {
        lines <- c(lines, vapply(grp, function(ko)
          paste0("ko:", ko, " #FFFFFF,", fg_color), character(1)))
      }
    } else {
      col <- if (use_outline) paste0(sample_color, ",", fg_color) else sample_color
      lines <- c(lines, vapply(grp, function(ko) paste0("ko:", ko, " ", col), character(1)))
    }
  }

  if (length(lines) == 0) return("")
  urlenc_spaces_newlines_only(paste(lines, collapse = "\n"))
}

build_multi_query_btex_outline_split <- function(pathway_kos, ko_universe, kos_by_sample, sample_colors, pathway_ko_groups = list(), fg_color = "#FF0000") {
  samples <- names(sample_colors)
  pathway_kos <- unique(pathway_kos[grepl("^K[0-9]{5}$", pathway_kos)])
  ko_universe <- unique(ko_universe[grepl("^K[0-9]{5}$", ko_universe)])
  target_kos <- sort(intersect(pathway_kos, ko_universe))
  if (length(target_kos) == 0) return("")

  use_outline <- !is.null(fg_color) && nzchar(fg_color) && tolower(fg_color) != "none"

  lines <- character(0)
  group_keys <- character(0)
  target_groups <- list()

  if (length(pathway_ko_groups) > 0) {
    for (g in pathway_ko_groups) {
      gg <- sort(unique(intersect(g, target_kos)))
      if (length(gg) == 0) next
      k <- paste(gg, collapse = "|")
      if (k %in% group_keys) next
      group_keys <- c(group_keys, k)
      target_groups <- c(target_groups, list(gg))
    }
  }

  covered <- if (length(target_groups) > 0) unique(unlist(target_groups, use.names = FALSE)) else character(0)
  for (ko in setdiff(target_kos, covered)) {
    target_groups <- c(target_groups, list(ko))
  }

  for (grp in target_groups) {
    present_samples <- samples[vapply(samples, function(s) any(grp %in% kos_by_sample[[s]]), logical(1))]
    if (length(present_samples) == 0) {
      if (use_outline) {
        # on pathway, absent from all samples: white fill + outline
        lines <- c(lines, vapply(grp, function(ko)
          paste0("ko:", ko, " #FFFFFF,", fg_color), character(1)))
      }
      # if no outline, skip absent KOs entirely (original behaviour)
    } else {
      if (use_outline) {
        cols <- paste(vapply(unname(sample_colors[present_samples]),
                             function(col) paste0(col, ",", fg_color), character(1)),
                      collapse = " ")
      } else {
        cols <- paste(unname(sample_colors[present_samples]), collapse = " ")
      }
      lines <- c(lines, vapply(grp, function(ko) paste0("ko:", ko, " ", cols), character(1)))
    }
  }

  if (length(lines) == 0) return("")
  urlenc_spaces_newlines_only(paste(lines, collapse = "\n"))
}

if (!is.null(opt$sample) && tolower(opt$sample) != "all") {
  samp <- opt$sample
  if (!samp %in% names(kos_by_sample)) stop("Requested --sample not found in input: ", samp)
  run_samples <- c(samp)
} else {
  all_samples <- sort(names(kos_by_sample))
  n_detected <- length(all_samples)
  message("Detected ", n_detected, " samples in input.")

  if (n_detected > 20) {
    run_samples <- all_samples[seq_len(20)]
    message("Using top 20 samples (alphabetical) for downstream processing due to palette limit.")
  } else {
    run_samples <- all_samples
  }
}

sample_colors <- get_sample_colors(run_samples)
if (length(run_samples) == 1) {
  sample_colors <- setNames("#377EB8", run_samples)
}
legend_path <- file.path(opt$outdir, "sample_color_legend.tsv")
write_sample_color_legend(legend_path, sample_colors)

write_kegg_html_launcher <- function(path, pid, multi_query_raw) {

  # Decode in reverse order of encoding (% must be decoded last)
  body_val <- gsub("%0A", "\n", multi_query_raw, fixed = TRUE)
  body_val <- gsub("%20", " ",  body_val,        fixed = TRUE)
  body_val <- gsub("%23", "#",  body_val,        fixed = TRUE)
  body_val <- gsub("%25", "%",  body_val,        fixed = TRUE)  # always last

  # Escape the three characters that are special inside JS template literals
  js_safe <- gsub("\\", "\\\\", body_val, fixed = TRUE)  # must be first
  js_safe <- gsub("`",  "\\`",  js_safe,  fixed = TRUE)
  js_safe <- gsub("$", "\\$", js_safe, fixed = TRUE)
  
  html <- paste0(
      '<!DOCTYPE html>\n',
      '<html><head><meta charset="utf-8"><title>KEGG map ', pid, '</title></head>\n',
      '<body>\n',
      '<p>Opening KEGG pathway map ', pid, '&#8230;</p>\n',
      '<p><button id="btn" onclick="go()">Click here if it did not open automatically</button></p>\n',
      '<form id="f" method="POST" action="https://www.kegg.jp/kegg-bin/show_pathway" target="_blank">\n',
      '  <input type="hidden" name="map"         value="map', pid, '">\n',
      '  <input type="hidden" name="multi_query" id="mq"  value="">\n',
      '</form>\n',
      '<script>\n',
      'function go(){\n',
      '  document.getElementById("mq").value = `', js_safe, '`;\n',
      '  document.getElementById("f").submit();\n',
      '}\n',
      'try { go(); } catch(e) {}\n',
      '</script>\n',
      '</body></html>\n'
    )

  writeLines(html, con = path)
}

for (pid in pids) {
  pathway_meta <- get_pathway_metadata(pid)
  if (isTRUE(opt$pathways_supplied)) {
    input_label <- if (identical(input_mode, "genome_dir")) "genome" else "input"
    message("[info] Scanning ", input_label, " for ", pathway_meta$filename_stub, ", ", pid)
  }
  fg_color_effective <- if (identical(input_mode, "genome_dir")) "none" else opt$fg_color
  kos_by_sample_run <- kos_by_sample[run_samples]
  pathway_kos <- get_pathway_kos_from_kgml(pid)
  pathway_ko_groups <- get_pathway_ko_groups_from_kgml(pid)

  if (length(pathway_kos) > 0) {
    target_kos <- sort(intersect(pathway_kos, ko_universe))
  } else {
    target_kos <- ko_universe
  }

  # ── build multi_query (same as before) ──────────────────────────────────────
  if (!is.null(opt$sample) && tolower(opt$sample) != "all") {
    mq <- build_multi_query_btex_outline_single(
      pathway_kos       = target_kos,
      ko_universe       = ko_universe,
      sample_kos        = kos_by_sample[[run_samples[1]]],
      sample_color      = sample_colors[[run_samples[1]]],
      pathway_ko_groups = pathway_ko_groups,
      fg_color          = fg_color_effective
    )
    label <- paste0("single sample: ", run_samples[1])
  } else {
        mq <- build_multi_query_btex_outline_split(
      pathway_kos       = target_kos,
      ko_universe       = ko_universe,
      kos_by_sample     = kos_by_sample_run,
      sample_colors     = sample_colors,
      pathway_ko_groups = pathway_ko_groups,
      fg_color          = fg_color_effective
    )
    label <- "split across samples"
  }

  # ── empty result ─────────────────────────────────────────────────────────────
  if (!nzchar(mq)) {
    message(
      "[warning] No hits detected for ",
      pathway_meta$filename_stub, "_", pid,
      " (", label, "); no HTML output generated."
    )
    next
  }

  # ── always write the HTML POST launcher ──────────────────────────────────────
  html_path <- file.path(opt$outdir, paste0(pathway_meta$filename_stub, "_", pid, ".html"))
  write_kegg_html_launcher(html_path, pid, mq)

  # ── presence table (unchanged) ───────────────────────────────────────────────
}

message("Done!")
message("Legend: ", legend_path)
