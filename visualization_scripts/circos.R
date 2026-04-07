#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(circlize)
  library(GenomicRanges)
})

ensure_dir <- function(p) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

hex_to_rgb_str <- function(hexcode) {
  h <- gsub("^#", "", trimws(hexcode))
  if (nchar(h) != 6) stop(sprintf("Bad hex color: %s", hexcode), call. = FALSE)
  vals <- c(
    strtoi(substr(h, 1, 2), 16L),
    strtoi(substr(h, 3, 4), 16L),
    strtoi(substr(h, 5, 6), 16L)
  )
  paste(vals, collapse = ",")
}

sanitize_color_name <- function(token) {
  token |>
    tolower() |>
    gsub("[^a-z0-9]+", "_", x = _)
}

normalize_color_token <- function(token) {
  t <- trimws(token)
  if (grepl(",", t, fixed = TRUE)) {
    parts <- trimws(strsplit(t, ",", fixed = TRUE)[[1]])
    if (length(parts) != 3) stop(sprintf("Unsupported color token: %s", token), call. = FALSE)
    vals <- suppressWarnings(as.integer(parts))
    if (any(is.na(vals)) || any(vals < 0 | vals > 255)) {
      stop(sprintf("Unsupported color token: %s", token), call. = FALSE)
    }
    return(paste(vals, collapse = ","))
  }
  if (startsWith(t, "#")) return(hex_to_rgb_str(t))
  stop(sprintf("Unsupported color token: %s", token), call. = FALSE)
}

rgb_str_to_hex <- function(rgb_str) {
  parts <- trimws(strsplit(rgb_str, ",", fixed = TRUE)[[1]])
  vals <- as.integer(parts)
  if (length(vals) != 3 || any(is.na(vals))) stop(sprintf("Bad rgb string: %s", rgb_str), call. = FALSE)
  grDevices::rgb(vals[1], vals[2], vals[3], maxColorValue = 255)
}

infer_default_pathway_map <- function() {
  candidate1 <- file.path("/home/juneq/BTEX-HMMs/btexhmm/data/pathway_map.tsv")
  if (file.exists(candidate1)) return(candidate1)
  stop(
    sprintf("[ERROR] Default pathway map not found. Tried:\n  %s\nUse --pathway-map to provide it explicitly.", candidate1),
    call. = FALSE
  )
}

read_fasta_lengths <- function(dna_path) {
  lines <- readLines(dna_path, warn = FALSE)
  headers <- which(startsWith(lines, ">"))
  if (length(headers) == 0) stop(sprintf("[ERROR] No FASTA headers found in %s", dna_path), call. = FALSE)

  contigs <- vector("list", length(headers))
  for (i in seq_along(headers)) {
    start <- headers[i]
    end <- if (i < length(headers)) headers[i + 1] - 1 else length(lines)
    header <- sub("^>", "", lines[start])
    contig_name <- strsplit(trimws(header), "\\s+")[[1]][1]
    seq_lines <- if (end > start) lines[(start + 1):end] else character(0)
    seq_len <- sum(nchar(gsub("\\s+", "", seq_lines)))
    contigs[[i]] <- data.frame(contig = contig_name, length = seq_len, stringsAsFactors = FALSE)
  }
  bind_rows(contigs)
}

read_fasta_sequences <- function(dna_path) {
  lines <- readLines(dna_path, warn = FALSE)
  headers <- which(startsWith(lines, ">"))
  if (length(headers) == 0) stop(sprintf("[ERROR] No FASTA headers found in %s", dna_path), call. = FALSE)

  seqs <- vector("list", length(headers))
  for (i in seq_along(headers)) {
    start <- headers[i]
    end <- if (i < length(headers)) headers[i + 1] - 1 else length(lines)
    header <- sub("^>", "", lines[start])
    contig_name <- strsplit(trimws(header), "\\s+")[[1]][1]
    seq_lines <- if (end > start) lines[(start + 1):end] else character(0)
    seq <- paste0(gsub("\\s+", "", seq_lines), collapse = "")
    seqs[[i]] <- tibble(contig = contig_name, sequence = seq)
  }
  bind_rows(seqs)
}

generate_contig_lengths_table_from_dna <- function(dna_path, genome, out_path) {
  if (!file.exists(dna_path)) stop(sprintf("[ERROR] DNA FASTA not found: %s", dna_path), call. = FALSE)
  entries <- read_fasta_lengths(dna_path)
  if (nrow(entries) == 0) stop(sprintf("[ERROR] No contig lengths gathered from %s", dna_path), call. = FALSE)
  out_df <- entries |>
    mutate(sample = genome) |>
    select(sample, contig, length)
  ensure_dir(dirname(out_path))
  write_tsv(out_df, out_path)
  message(sprintf("[info] wrote %d contig entries to %s", nrow(out_df), out_path))
  out_path
}

ensure_contig_lengths <- function(opt, output_dir) {
  default_out <- file.path(output_dir, "contig_length.tsv")
  ensure_dir(output_dir)

  if (!is.null(opt$contig_lengths) && nzchar(opt$contig_lengths)) {
    candidate <- opt$contig_lengths
    if (file.exists(candidate)) {
      if (normalizePath(candidate, winslash = "/") != normalizePath(default_out, winslash = "/", mustWork = FALSE)) {
        file.copy(candidate, default_out, overwrite = TRUE)
        message(sprintf("[info] Copied contig lengths to %s", default_out))
      } else {
        message(sprintf("[info] Using contig lengths at %s", default_out))
      }
      return(default_out)
    }
    if (is.null(opt$dna) || !nzchar(opt$dna)) {
      stop(sprintf("[ERROR] contig lengths file not found: %s. Provide --dna so it can be generated automatically.", candidate), call. = FALSE)
    }
    message(sprintf("[info] contig lengths not found at %s; generating -> %s", candidate, default_out))
    return(generate_contig_lengths_table_from_dna(opt$dna, opt$sample, default_out))
  }

  if (is.null(opt$dna) || !nzchar(opt$dna)) {
    stop("[ERROR] Provide --dna or an existing --contig-lengths file.", call. = FALSE)
  }

  message(sprintf("[info] Auto-generating contig lengths -> %s", default_out))
  generate_contig_lengths_table_from_dna(opt$dna, opt$sample, default_out)
}

read_contig_lengths <- function(path, genome) {
  read_tsv(path, show_col_types = FALSE) |>
    filter(sample == genome) |>
    mutate(length = as.integer(length)) |>
    select(contig, length)
}

normalize_contig_id <- function(token) {
  if (is.null(token) || !nzchar(token)) return("")
  name <- token |>
    trimws() |>
    sub("^>", "", x = _) |>
    trimws()
  if (grepl(" #", name, fixed = TRUE)) name <- strsplit(name, " #", fixed = TRUE)[[1]][1]
  if (grepl("|", name, fixed = TRUE)) {
    parts <- strsplit(name, "|", fixed = TRUE)[[1]]
    name <- parts[length(parts)]
  }
  trimws(name)
}

parse_hit_headers <- function(sample, hmm, headers) {
  if (is.na(headers) || !nzchar(headers)) return(tibble())
  raw_parts <- strsplit(headers, "|", fixed = TRUE)[[1]]
  out <- vector("list", length(raw_parts))
  keep <- 0

  for (raw in raw_parts) {
    if (!nzchar(trimws(raw))) next
    hash_fields <- trimws(strsplit(raw, "#", fixed = TRUE)[[1]])
    left <- hash_fields[1]
    contig <- normalize_contig_id(left)
    if (grepl(" #", raw, fixed = TRUE)) {
      contig <- sub("_\\d+$", "", contig)
    }
    if (!nzchar(contig)) next

    if (length(hash_fields) >= 4) {
      nums <- suppressWarnings(as.integer(hash_fields[2:3]))
      strand_num <- suppressWarnings(as.integer(hash_fields[4]))
    } else {
      nums <- as.integer(stringr::str_extract_all(raw, "\\b(\\d+)\\b")[[1]])
      strand_num <- if (length(nums) >= 2 && nums[1] > nums[2]) -1L else 1L
    }
    if (length(nums) < 2 || any(is.na(nums[1:2]))) next
    if (is.na(strand_num) || !strand_num %in% c(-1L, 1L)) {
      strand_num <- if (nums[1] > nums[2]) -1L else 1L
    }

    a <- nums[1]
    b <- nums[2]
    keep <- keep + 1
    out[[keep]] <- tibble(
      sample = sample,
      hmm = hmm,
      contig = contig,
      start = min(a, b),
      end = max(a, b),
      strand = strand_num
    )
  }

  if (keep == 0) return(tibble())
  bind_rows(out[seq_len(keep)])
}

load_pathway_map <- function(pathway_map_path) {
  pm <- read_tsv(pathway_map_path, show_col_types = FALSE)
  required <- c("hmm", "pathway", "color")
  missing <- setdiff(required, colnames(pm))
  if (length(missing) > 0) {
    stop(sprintf("[ERROR] Pathway map missing required columns: %s", paste(sort(missing), collapse = ", ")), call. = FALSE)
  }

  pm <- pm |>
    mutate(
      hmm = trimws(hmm),
      pathway = trimws(pathway),
      color = trimws(color)
    )

  bad <- pm |>
    filter(!nzchar(hmm) | !nzchar(pathway) | !nzchar(color))
  if (nrow(bad) > 0) {
    stop("[ERROR] Pathway map row must include hmm, pathway, and color values.", call. = FALSE)
  }
  pm
}

build_color_map <- function(hmms_with_hits, pathway_df) {
  pathway_map <- pathway_df |>
    select(hmm, pathway) |>
    distinct()
  pathway_colors <- pathway_df |>
    select(pathway, color) |>
    distinct() |>
    mutate(rgb = map_chr(color, normalize_color_token))

  final <- lapply(hmms_with_hits, function(hmm) {
    pathway_row <- pathway_map |> filter(tolower(hmm) == tolower(!!hmm))
    if (nrow(pathway_row) == 0) {
      stop(sprintf("[ERROR] No pathway mapping found for HMM: %s", hmm), call. = FALSE)
    }
    pathway <- pathway_row$pathway[1]
    color_row <- pathway_colors |> filter(tolower(pathway) == tolower(!!pathway))
    if (nrow(color_row) == 0) {
      stop(sprintf("[ERROR] No color defined in pathway_map.tsv for pathway '%s' (HMM: %s)", pathway, hmm), call. = FALSE)
    }
    tibble(
      hmm = hmm,
      color_name = sanitize_color_name(pathway),
      rgb = color_row$rgb[1],
      color_hex = rgb_str_to_hex(color_row$rgb[1]),
      pathway = pathway
    )
  })
  bind_rows(final)
}

compute_gc_skew_windows <- function(dna_path, contigs_to_keep = NULL, window_size = 5000L, step_size = 5000L) {
  seq_df <- read_fasta_sequences(dna_path)
  if (!is.null(contigs_to_keep)) {
    seq_df <- seq_df |>
      filter(contig %in% contigs_to_keep)
  }

  out <- list()
  idx <- 1L

  for (i in seq_len(nrow(seq_df))) {
    contig <- seq_df$contig[i]
    seq <- toupper(seq_df$sequence[i])
    n <- nchar(seq)
    if (n == 0) next

    starts <- seq.int(1L, n, by = step_size)
    for (s in starts) {
      e <- min(s + window_size - 1L, n)
      subseq <- substr(seq, s, e)
      chars <- strsplit(subseq, "", fixed = TRUE)[[1]]
      g_count <- sum(chars == "G")
      c_count <- sum(chars == "C")
      denom <- g_count + c_count
      skew <- if (denom == 0L) 0 else (g_count - c_count) / denom

      out[[idx]] <- tibble(
        contig = contig,
        start = s,
        end = e,
        gc_skew = skew
      )
      idx <- idx + 1L
    }
  }

  if (length(out) == 0) {
    return(tibble(contig = character(), start = integer(), end = integer(), gc_skew = numeric()))
  }
  bind_rows(out)
}

write_hmm_colors_tsv <- function(path, cmap) {
  if (nrow(cmap) == 0) {
    write_tsv(tibble(hmm = character(), color_name = character(), rgb = character()), path)
    return(invisible(NULL))
  }
  write_tsv(cmap |> select(hmm, color_name, rgb), path)
}

write_karyotype_tsv <- function(path, contig_lengths) {
  write_tsv(contig_lengths, path)
}

write_gene_hits_tsv <- function(path, features, cmap) {
  out <- features |>
    left_join(cmap |> select(hmm, color_name, rgb, pathway), by = "hmm")
  write_tsv(out, path)
}

write_gene_labels_tsv <- function(path, features) {
  out <- features |>
    distinct(contig, start, end, hmm, strand)
  write_tsv(out, path)
}

write_genbank_like_hits_tsv <- function(path, features, contig_lengths) {
  out <- features |>
    left_join(contig_lengths, by = c("contig" = "contig")) |>
    arrange(contig, start, end, hmm)
  write_tsv(out, path)
}

write_gc_skew_tsv <- function(path, gc_skew_df) {
  write_tsv(gc_skew_df, path)
}

prepare_sector_data <- function(contig_lengths) {
  contig_lengths |>
    mutate(contig = factor(contig, levels = contig))
}

plot_ideogram_track <- function(sector_df, tick_major = 5e5, tick_minor = 1e5) {
  circos.trackPlotRegion(
    ylim = c(0, 1),
    bg.border = NA,
    track.height = 0.08,
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")

      circos.rect(xlim[1], 0, xlim[2], 1, col = "lightblue", border = "black", lwd = 1)

      max_x <- xlim[2]
      major_ticks <- seq(0, max_x, by = tick_major)
      minor_ticks <- seq(0, max_x, by = tick_minor)

      if (length(minor_ticks) > 0) {
        circos.segments(minor_ticks, 1, minor_ticks, 1.08, col = "black", lwd = 0.5)
      }
      if (length(major_ticks) > 0) {
        circos.segments(major_ticks, 1, major_ticks, 1.15, col = "black", lwd = 1)
        circos.text(
          x = major_ticks,
          y = rep(1.22, length(major_ticks)),
          labels = paste0(round(major_ticks / 1e3, 1), " kb"),
          cex = 0.45,
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0, 0.5)
        )
      }

      circos.text(
        x = mean(xlim),
        y = 2.2,
        labels = sector.index,
        cex = 0.8,
        facing = "bending.outside",
        niceFacing = TRUE,
        font = 1
      )
    }
  )
}

plot_gc_skew_track <- function(gc_skew_df, track_height = 0.10) {
  if (nrow(gc_skew_df) == 0) return(invisible(NULL))

  max_abs_skew <- max(abs(gc_skew_df$gc_skew), na.rm = TRUE)
  plot_limit <- max(0.05, max_abs_skew)

  circos.genomicTrackPlotRegion(
    gc_skew_df[, c("contig", "start", "end", "gc_skew")],
    ylim = c(-plot_limit, plot_limit),
    bg.border = "black",
    track.height = track_height,
    panel.fun = function(region, value, ...) {
      if (nrow(region) == 0) return()
      skew_vals <- value[[1]]
      colors <- ifelse(skew_vals >= 0, "#E67E22", "#2E86C1")
      circos.genomicRect(
        region,
        value,
        ybottom = pmin(skew_vals, 0),
        ytop = pmax(skew_vals, 0),
        col = colors,
        border = NA,
        ...
      )
      xlim <- get.cell.meta.data("xlim")
      circos.segments(xlim[1], 0, xlim[2], 0, lty = 2, col = "grey35")
    }
  )
}

plot_feature_track <- function(track_features, cmap, track_height = 0.08, label_cex = 0.45, show_links = FALSE) {
  if (nrow(track_features) == 0) {
    circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = track_height, panel.fun = function(x, y) {})
    return(invisible(NULL))
  }

  feat <- track_features |>
    left_join(cmap |> select(hmm, color_hex), by = "hmm")

  circos.genomicTrackPlotRegion(
    feat[, c("contig", "start", "end", "hmm", "color_hex")],
    ylim = c(0, 1),
    bg.border = NA,
    track.height = track_height,
    panel.fun = function(region, value, ...) {
      if (nrow(region) == 0) return()
      circos.genomicRect(region, value, ytop = 1, ybottom = 0, col = value[["color_hex"]], border = NA, ...)
    }
  )

  lab_df <- feat |>
    distinct(contig, start, end, hmm) %>%
    mutate(chr = contig) %>%
    select(chr, start, end, hmm)

  if (nrow(lab_df) > 0) {
    circos.genomicLabels(
      lab_df,
      labels.column = 4,
      labels.side = "inside",
      cex = label_cex,
      col = "black",
      line_col = "black",
      connection_height = mm_h(5)
    )
  }
}

make_circlize_plot <- function(sample_dir, genome, contig_lengths, features, cmap, gc_skew_df = NULL, pdf_width = 10, pdf_height = 10) {
  pdf_path <- file.path(sample_dir, "circos_plot.pdf")
  png_path <- file.path(sample_dir, "circos_plot.png")

  if (!is.numeric(pdf_width) || length(pdf_width) != 1 || is.na(pdf_width) || !is.finite(pdf_width) || pdf_width <= 0) {
    stop(sprintf("[ERROR] --pdf-width must be a positive finite number; got: %s", as.character(pdf_width)), call. = FALSE)
  }
  if (!is.numeric(pdf_height) || length(pdf_height) != 1 || is.na(pdf_height) || !is.finite(pdf_height) || pdf_height <= 0) {
    stop(sprintf("[ERROR] --pdf-height must be a positive finite number; got: %s", as.character(pdf_height)), call. = FALSE)
  }

  sector_df <- prepare_sector_data(contig_lengths)
  xlim_mat <- cbind(rep(0, nrow(sector_df)), sector_df$length)
  rownames(xlim_mat) <- sector_df$contig

  draw_once <- function(device_fun) {
    device_fun()
    on.exit(dev.off(), add = TRUE)

    circos.clear()
    circos.par(
      start.degree = 90,
      gap.degree = rep(1, nrow(sector_df)),
      cell.padding = c(0, 0, 0, 0),
      track.margin = c(0.01, 0.01),
      points.overflow.warning = FALSE
    )

    circos.initialize(factors = as.character(sector_df$contig), xlim = xlim_mat)
    plot_ideogram_track(sector_df)

    if (!is.null(gc_skew_df) && nrow(gc_skew_df) > 0) {
      plot_gc_skew_track(gc_skew_df, track_height = 0.10)
    }
    plot_feature_track(features, cmap, track_height = 0.09, label_cex = 0.45, show_links = FALSE)

    text(
      x = 0,
      y = 0,
      labels = genome,
      cex = 1.2,
      font = 2
    )
    circos.clear()
  }

  draw_once(function() pdf(pdf_path, width = pdf_width, height = pdf_height, useDingbats = FALSE))
  draw_once(function() {
    png_width_px <- max(1L, as.integer(round(pdf_width * 300)))
    png_height_px <- max(1L, as.integer(round(pdf_height * 300)))
    png(png_path, width = png_width_px, height = png_height_px, res = 300)
  })

  invisible(list(pdf = pdf_path, png = png_path))
}

run <- function(opt) {
  genome <- opt$sample
  if (is.null(genome) || !nzchar(genome)) stop("[ERROR] --sample is required", call. = FALSE)
  if (is.null(opt$hmmscan) || !file.exists(opt$hmmscan)) stop(sprintf("[ERROR] hmmscan file not found: %s", opt$hmmscan), call. = FALSE)
  if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("[ERROR] --outdir is required", call. = FALSE)

  pathway_map_path <- if (!is.null(opt$pathway_map) && nzchar(opt$pathway_map)) opt$pathway_map else infer_default_pathway_map()
  if (!file.exists(pathway_map_path)) stop(sprintf("[ERROR] Pathway map file not found: %s", pathway_map_path), call. = FALSE)
  message(sprintf("[info] Using pathway map: %s", pathway_map_path))

  sample_dir <- file.path(opt$outdir, paste0(genome, "_circos_plot"))
  ensure_dir(sample_dir)

  pathway_df <- load_pathway_map(pathway_map_path)

  contig_lengths_path <- ensure_contig_lengths(opt, sample_dir)
  contig_lengths <- read_contig_lengths(contig_lengths_path, genome)
  if (nrow(contig_lengths) == 0) stop(sprintf("[ERROR] No contig lengths found for genome: %s", genome), call. = FALSE)

  hmmscan_df <- read_csv(opt$hmmscan, show_col_types = FALSE)
  required_cols <- c("sample", "hmm", "hits", "hit_headers")
  missing <- setdiff(required_cols, colnames(hmmscan_df))
  if (length(missing) > 0) {
    stop(sprintf("[ERROR] hmmscan CSV missing required columns: %s", paste(sort(missing), collapse = ", ")), call. = FALSE)
  }

  selected_rows <- hmmscan_df |>
    filter(sample == genome, as.integer(hits) >= 1)

  features_list <- pmap(
    list(selected_rows$sample, selected_rows$hmm, selected_rows$hit_headers),
    parse_hit_headers
  )
  features <- bind_rows(features_list)

  hmms_with_hits <- sort(unique(features$hmm))

  cmap <- if (length(hmms_with_hits) > 0) build_color_map(hmms_with_hits, pathway_df) else tibble(
    hmm = character(), color_name = character(), rgb = character(), color_hex = character(), pathway = character()
  )

  gc_skew_df <- compute_gc_skew_windows(
    dna_path = opt$dna,
    contigs_to_keep = contig_lengths$contig,
    window_size = 5000L,
    step_size = 5000L
  )

  write_hmm_colors_tsv(file.path(sample_dir, "hmm_colors.tsv"), cmap)
  write_karyotype_tsv(file.path(sample_dir, "karyotype.tsv"), contig_lengths)
  write_gene_hits_tsv(file.path(sample_dir, "gene_hits.tsv"), features, cmap)
  write_gene_labels_tsv(file.path(sample_dir, "gene_labels.tsv"), features)
  write_genbank_like_hits_tsv(file.path(sample_dir, "gene_hits_export.tsv"), features, contig_lengths)
  write_gc_skew_tsv(file.path(sample_dir, "gc_skew_5kb.tsv"), gc_skew_df)

  plot_paths <- make_circlize_plot(
    sample_dir = sample_dir,
    genome = genome,
    contig_lengths = contig_lengths,
    features = features,
    cmap = cmap,
    gc_skew_df = gc_skew_df,
    pdf_width = opt$pdf_width,
    pdf_height = opt$pdf_height
  )

  message(sprintf("[OK] circlize project for %s -> %s", genome, sample_dir))
  message(sprintf("[OK] Plot PDF -> %s", plot_paths$pdf))
  message(sprintf("[OK] Plot PNG -> %s", plot_paths$png))

  invisible(sample_dir)
}

option_list <- list(
  make_option("--hmmscan", type = "character", help = "Path to btex_hmm_summary.csv or similar hmmscan summary CSV"),
  make_option("-o", "--outdir", type = "character", help = "Output directory"),
  make_option("--dna", type = "character", default = NULL, help = "Genome FASTA file for the selected sample"),
  make_option("-s", "--sample", type = "character", help = "Sample/genome name to plot"),
  make_option("--contig-lengths", dest = "contig_lengths", type = "character", default = NULL, help = "Optional TSV with columns: sample, contig, length"),
  make_option("--pathway-map", dest = "pathway_map", type = "character", default = NULL, help = "Optional override for pathway_map.tsv"),
  make_option("--pdf-width", dest = "pdf_width", type = "double", default = 10, help = "Output PDF/PNG width in inches [default %default]"),
  make_option("--pdf-height", dest = "pdf_height", type = "double", default = 10, help = "Output PDF/PNG height in inches [default %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$hmmscan) || is.null(opt$outdir) || is.null(opt$sample)) {
  print_help(parser)
  quit(status = 1)
}

run(opt)
