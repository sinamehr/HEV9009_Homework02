#!/usr/bin/env Rscript

active <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else Sys.getenv("MMC_GRAPH_ACTIVE_SCRIPT")
})
if (!nzchar(active)) stop("Cannot determine the sunburst graph script path.", call. = FALSE)
source(file.path(dirname(normalizePath(active)), "0103030408_mmc_v01_graph_helpers.inc"), local = TRUE)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Required R package is missing: ggplot2", call. = FALSE)

publication_font <- "Arial"

mg_mix_white <- function(color, amount) {
  rgb <- grDevices::col2rgb(color) / 255
  grDevices::rgb(rgb[1L, ] + (1 - rgb[1L, ]) * amount,
                 rgb[2L, ] + (1 - rgb[2L, ]) * amount,
                 rgb[3L, ] + (1 - rgb[3L, ]) * amount)
}

dimension_colors <- mg_dimension_colors
indicator_order <- list(
  Technology = paste0("w_technology_ai_", c("tml", "ttm", "tpa", "tir", "tsr", "tnlg", "tar")),
  Purpose = paste0("w_purpose_ai_", c("pms", "ppp", "plog", "pits", "pbam_harmonised", "pba", "pme", "phr", "pbam", "pfin", "prdi")),
  Acquisition = paste0("w_acquisition_ai_", c("adown", "amown", "aos", "ardy", "aext"))
)
indicator_colors <- unlist(lapply(names(indicator_order), function(dimension) {
  values <- indicator_order[[dimension]]
  setNames(mg_mix_white(dimension_colors[[dimension]], seq(0.08, 0.55, length.out = length(values))), values)
}), use.names = TRUE)
fill_colors <- c(setNames(dimension_colors, paste0("dimension_", tolower(names(dimension_colors)))), indicator_colors)

dimension_codes <- c(dimension_technology = "T", dimension_purpose = "P", dimension_acquisition = "A")
indicator_codes <- unlist(lapply(names(indicator_order), function(dimension) {
  values <- indicator_order[[dimension]]
  setNames(paste0(substr(dimension, 1L, 1L), seq_along(values)), values)
}), use.names = TRUE)
fill_codes <- c(dimension_codes, indicator_codes)

indicator_label <- function(weight) {
  label <- toupper(sub("^w_(?:technology|purpose|acquisition)_ai_", "", weight))
  ifelse(label == "PBAM_HARMONISED", "PBAM-H", label)
}
fill_labels <- c(
  setNames(unname(mg_dimension_legend_labels()[mg_dimension_levels]), names(dimension_codes)),
  mg_indicator_legend_labels(names(indicator_codes), unname(indicator_codes))
)

mg_relative_luminance <- function(colors) {
  rgb <- grDevices::col2rgb(colors) / 255
  rgb <- ifelse(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  as.numeric(c(0.2126, 0.7152, 0.0722) %*% rgb)
}

mg_text_color <- function(colors) {
  luminance <- mg_relative_luminance(colors)
  black_contrast <- (luminance + 0.05) / 0.05
  white_contrast <- 1.05 / (luminance + 0.05)
  ifelse(black_contrast >= white_contrast, "black", "white")
}

mg_validate_palette <- function(colors) {
  if (anyNA(grDevices::col2rgb(colors))) stop("The publication palette contains an invalid color.", call. = FALSE)
  luminance <- mg_relative_luminance(colors)
  best_contrast <- pmax((luminance + 0.05) / 0.05, 1.05 / (luminance + 0.05))
  if (any(best_contrast < 4.5)) stop("The publication palette fails the text-contrast preflight.", call. = FALSE)
  invisible(TRUE)
}

build_geometry <- function(record) {
  rows <- list(); centers <- list(); index <- 0L
  data <- record$data
  for (row_number in seq_len(nrow(data))) {
    group_raw <- if (record$grouped) as.character(data[[record$group_column]][[row_number]]) else ""
    group_title <- if (!nzchar(group_raw)) "" else if (grepl("group", group_raw, ignore.case = TRUE)) group_raw else paste("Group", group_raw)
    fit_title <- paste(c(record$year, record$scope, record$method, group_title)[nzchar(c(record$year, record$scope, record$method, group_title))], collapse = " | ")
    panel_title <- if (record$grouped) paste0(LETTERS[[row_number]], ". ", group_title) else ""
    if (identical(record$structure, "flat")) {
      weights <- record$spec$flat
      values <- mg_normalize(unlist(data[row_number, weights], use.names = FALSE), paste(fit_title, "flat weights"))
      ends <- cumsum(values)
      starts <- c(0, head(ends, -1L))
      for (item in seq_along(weights)) {
        index <- index + 1L
        rows[[index]] <- data.frame(
          panel = panel_title, level = "Weight", label = toupper(sub("^w_", "", weights[[item]])),
          fill_key = weights[[item]], ymin = starts[[item]], ymax = ends[[item]],
          xmin = 0.78, xmax = 1.68, stringsAsFactors = FALSE
        )
      }
      centers[[row_number]] <- data.frame(panel = panel_title, x = 0, y = 0, label = "AI\nAdoption")
      next
    }
    dimension_weights <- mg_normalize(unlist(data[row_number, record$spec$dimensions], use.names = FALSE), paste(fit_title, "dimension weights"))
    names(dimension_weights) <- mg_dimension_levels
    dimension_end <- cumsum(dimension_weights)
    dimension_start <- c(0, head(dimension_end, -1L))
    names(dimension_start) <- names(dimension_weights)
    for (dimension in mg_dimension_levels) {
      index <- index + 1L
      rows[[index]] <- data.frame(
        panel = panel_title, level = "Dimension", label = dimension,
        fill_key = paste0("dimension_", tolower(dimension)),
        ymin = dimension_start[[dimension]], ymax = dimension_end[[dimension]],
        xmin = 0.62, xmax = 1.42, stringsAsFactors = FALSE
      )
      weights <- record$spec[[tolower(dimension)]]
      within <- mg_normalize(unlist(data[row_number, weights], use.names = FALSE), paste(fit_title, dimension, "indicator weights"))
      outer_end <- dimension_start[[dimension]] + cumsum(within) * dimension_weights[[dimension]]
      outer_start <- c(dimension_start[[dimension]], head(outer_end, -1L))
      for (item in seq_along(weights)) {
        index <- index + 1L
        rows[[index]] <- data.frame(
          panel = panel_title, level = "Indicator", label = indicator_label(weights[[item]]),
          fill_key = weights[[item]], ymin = outer_start[[item]], ymax = outer_end[[item]],
          xmin = 1.42, xmax = 2.32, stringsAsFactors = FALSE
        )
      }
    }
    centers[[row_number]] <- data.frame(panel = panel_title, x = 0, y = 0, label = "AI\nAdoption")
  }
  list(rectangles = do.call(rbind, rows), centers = do.call(rbind, centers))
}

plot_record <- function(record) {
  geometry <- build_geometry(record)
  panels <- unique(geometry$rectangles$panel)
  if (identical(record$structure, "hierarchical")) {
    active_colors <- setNames(dimension_colors, paste0("dimension_", tolower(mg_dimension_levels)))
    active_labels <- setNames(unname(mg_dimension_legend_labels()[mg_dimension_levels]), names(active_colors))
    active_codes <- dimension_codes
    indicator_counts <- integer(length(mg_dimension_levels))
    names(indicator_counts) <- mg_dimension_levels
    for (dimension in mg_dimension_levels) {
      weights <- record$spec[[tolower(dimension)]]
      indicator_counts[[dimension]] <- length(weights)
      codes <- paste0(substr(dimension, 1L, 1L), seq_along(weights))
      shades <- mg_mix_white(
        dimension_colors[[dimension]],
        seq(0.08, 0.55, length.out = length(weights))
      )
      if (length(shades) > 1L && any(diff(mg_relative_luminance(shades)) <= 0)) {
        stop(dimension, " indicator shades must progress from darkest to lightest.", call. = FALSE)
      }
      active_colors <- c(active_colors, setNames(shades, weights))
      active_labels <- c(active_labels, mg_indicator_legend_labels(weights, codes))
      active_codes <- c(active_codes, setNames(codes, weights))
    }
    max_indicators <- max(indicator_counts)
    dimension_keys <- paste0("dimension_", tolower(mg_dimension_levels))
    gap_keys <- paste0("legend_gap_", seq_along(mg_dimension_levels))
    indicator_rows <- lapply(seq_len(max_indicators), function(item) {
      vapply(mg_dimension_levels, function(dimension) {
        weights <- record$spec[[tolower(dimension)]]
        if (item <= length(weights)) weights[[item]] else paste0("legend_pad_", tolower(dimension), "_", item)
      }, character(1L))
    })
    legend_keys <- c(dimension_keys, gap_keys, unlist(indicator_rows, use.names = FALSE))
    spacer_keys <- setdiff(legend_keys, names(active_colors))
    legend_colors <- c(active_colors, setNames(rep("transparent", length(spacer_keys)), spacer_keys))
    legend_labels <- c(active_labels, setNames(rep("", length(spacer_keys)), spacer_keys))
    geometry$rectangles$sector_code <- unname(active_codes[geometry$rectangles$fill_key])
    legend_ncol <- length(mg_dimension_levels)
    legend_byrow <- TRUE
  } else {
    active_keys <- unique(geometry$rectangles$fill_key)
    active_colors <- fill_colors[active_keys]
    active_labels <- fill_labels[active_keys]
    legend_keys <- active_keys
    legend_colors <- active_colors
    legend_labels <- active_labels
    geometry$rectangles$sector_code <- unname(fill_codes[geometry$rectangles$fill_key])
    legend_ncol <- min(3L, length(active_keys))
    legend_byrow <- TRUE
  }
  missing_keys <- setdiff(unique(geometry$rectangles$fill_key), names(active_colors))
  if (length(missing_keys)) {
    stop("No reviewed publication color/code is registered for: ", paste(missing_keys, collapse = ", "), call. = FALSE)
  }
  mg_validate_palette(active_colors)
  geometry$rectangles$text_color <- mg_text_color(unname(active_colors[geometry$rectangles$fill_key]))
  geometry$rectangles$label_x <- (geometry$rectangles$xmin + geometry$rectangles$xmax) / 2
  geometry$rectangles$label_y <- (geometry$rectangles$ymin + geometry$rectangles$ymax) / 2
  label_data <- mg_stack_overlapping_polar_labels(
    geometry$rectangles, c("panel", "level")
  )
  plot <- ggplot2::ggplot(geometry$rectangles) +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_key),
                       color = "white", linewidth = 0.25) +
    ggplot2::geom_text(data = label_data,
                       ggplot2::aes(x = label_x, y = label_y, label = sector_code, color = text_color),
                       inherit.aes = FALSE, family = publication_font, size = 2.47,
                       lineheight = 0.9, show.legend = FALSE) +
    ggplot2::geom_text(data = geometry$centers,
                       ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE,
                       family = publication_font, fontface = "bold", size = 2.7, lineheight = 0.9) +
    ggplot2::coord_polar(theta = "y", start = 0, direction = 1, clip = "off") +
    ggplot2::xlim(0, 2.35) +
    ggplot2::scale_fill_manual(
      values = legend_colors, limits = legend_keys, breaks = legend_keys,
      labels = unname(mg_wrap_legend_labels(legend_labels[legend_keys])), drop = FALSE
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::guides(fill = ggplot2::guide_legend(
      title = NULL, ncol = legend_ncol, byrow = legend_byrow,
      keywidth = grid::unit(3.8, "mm"), keyheight = grid::unit(3.8, "mm")
    )) +
    ggplot2::theme_void(base_family = publication_font, base_size = 8) +
    ggplot2::theme(
      legend.position = "bottom", legend.justification = "center", legend.box.just = "center",
      legend.text = ggplot2::element_text(size = 7, angle = 0),
      legend.key.width = grid::unit(3.8, "mm"), legend.key.height = grid::unit(3.8, "mm"),
      legend.box = "vertical", legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 8, face = "bold"), aspect.ratio = 1
    ) +
    mg_graph_frame(8)
  if (record$grouped) {
    columns <- if (length(panels) == 6L) 3L else min(3L, length(panels))
    plot <- plot + ggplot2::facet_wrap(ggplot2::vars(panel), ncol = columns, nrow = ceiling(length(panels) / columns))
  }
  plot
}

write_plot <- function(record, output_dir, settings) {
  plot <- mg_single_method_badge(plot_record(record), record$method)
  output_stem <- file.path(output_dir, paste(
    "sun", tolower(record$method), mg_short_scope(record$scope),
    mg_short_dataset(record$dataset), sep = "_"
  ))
  width <- 20 + if (record$grouped) as.numeric(settings$grouped_width_mm) else as.numeric(settings$ungrouped_width_mm)
  height <- 20 + if (record$grouped) as.numeric(settings$grouped_height_mm) else as.numeric(settings$ungrouped_height_mm)
  formats <- unlist(settings$formats, use.names = FALSE)
  written <- character()
  caption <- paste(
    "Sunburst fitted weights:", record$year, record$scope, record$method,
    if (record$grouped) paste0("groups: ", paste(as.character(record$data[[record$group_column]]), collapse = ", ")) else "Global",
    paste0("Source: ", basename(record$source)), sep = " | "
  )
  writeLines(caption, paste0(output_stem, "_caption.txt"), useBytes = TRUE)
  if ("pdf" %in% formats) {
    path <- paste0(output_stem, ".pdf")
    if (!capabilities("cairo")) stop("Cairo is required for publication PDF output.", call. = FALSE)
    ggplot2::ggsave(path, plot = plot, device = grDevices::cairo_pdf,
                    width = width, height = height, units = "mm", bg = "white")
    written <- c(written, path)
  }
  if ("tiff" %in% formats) {
    path <- paste0(output_stem, ".tiff")
    ggplot2::ggsave(path, plot = plot, device = "tiff", width = width, height = height, units = "mm",
                    dpi = as.numeric(settings$tiff_dpi), compression = "lzw", type = "cairo", bg = "white")
    written <- c(written, path)
  }
  if (!length(written)) stop("No supported graph output format is configured.", call. = FALSE)
  written
}

config <- mg_config()
inputs <- mg_input_paths(config)
output_dir <- mg_script_output_dir(config, active)
candidates <- mg_candidates(inputs)
if (!length(candidates)) stop("No safe aggregate fitted-weight candidates were found.", call. = FALSE)
records <- unlist(lapply(candidates, mg_records), recursive = FALSE)
if (!length(records)) stop("No input contained a complete supported fitted-weight set.", call. = FALSE)

keys <- vapply(records, function(item) item$fit_key, character(1L))
selected <- lapply(split(seq_along(records), keys), function(indices) {
  priorities <- vapply(records[indices], function(item) item$priority, integer(1L))
  records[[indices[[which.max(priorities)]]]]
})
records <- unname(selected)

created <- unlist(lapply(records, write_plot, output_dir = output_dir,
                         settings = config$value$graph_generation), use.names = FALSE)
message("Created ", length(created), " sunburst figure file(s) from ", length(records), " fitted weight input(s).")
