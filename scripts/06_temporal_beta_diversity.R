# ============================================================================ #
# 06_temporal_beta_diversity.R
# ============================================================================ #
#
# Purpose:
# Reproduce the temporal beta-diversity analyses used in Figure 5:
#
#   - Figure 5A: weighted UniFrac similarity to the initial community (T0)
#   - Figure 5B: Bray-Curtis similarity to the initial community (T0)
#   - Figure 5C: pairwise temporal distance vs weighted UniFrac distance
#   - Figure 5D: pairwise temporal distance vs Bray-Curtis dissimilarity
#
# Notes:
# - Analyses use the non-rarefied final phyloseq object.
# - Baseline similarity is computed separately for each permanent sampling
#   point using that point's earliest available sample as T0.
# - Similarity is calculated as (1 - distance) * 100.
# - Pairwise temporal turnover is assessed within each basin.
# - Mantel tests use Spearman correlation and 9,999 permutations.
# - GLMs are used only for the fitted visual trend and slope annotation,
#   matching the original analysis code.
# - The wet-period shading in the T0 plots follows the original plotting
#   workflow: 0-4 days in HCP and 0-28 days in CHF and YSB, relative to the
#   first sampled date within each sampling-point time series.
#
# Input:
# - output_phyloseq_tree/ps_mdecon_v2_FTree.rds
#
# Main outputs:
# - Figure5A_similarity_T0_wunifrac.*
# - Figure5B_similarity_T0_bray.*
# - Figure5C_temporal_distance_wunifrac.*
# - Figure5D_temporal_distance_bray.*
# - baseline-similarity tables
# - pairwise temporal-distance tables
# - Mantel-test tables
# - GLM slope tables
#
# ============================================================================ #


# ============================================================================ #
# Packages
# ============================================================================ #

library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(cowplot)
library(broom)
library(scales)
library(openxlsx)


# ============================================================================ #
# Reproducibility
# ============================================================================ #

set.seed(2026)


# ============================================================================ #
# User paths
# ============================================================================ #

phyloseq_path <- file.path(
  "output_phyloseq_tree",
  "ps_mdecon_v2_FTree.rds"
)

output_dir <- "output_temporal_beta_diversity"
figure_dir <- file.path(
  output_dir,
  "figures"
)
table_dir <- file.path(
  output_dir,
  "tables"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  figure_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  table_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================================ #
# Shared settings
# ============================================================================ #

basin_levels <- c(
  "Herradura Clay Pan",
  "Ckoirama Halite Field",
  "Yungay Station Basin"
)

sampling_point_levels <- c(
  "HCP1",
  "HCP2",
  "HCP3",
  "CHF1",
  "CHF2",
  "CHF3",
  "YSB1",
  "YSB2",
  "YSB3"
)

sampling_point_colors <- c(
  "HCP1" = "#006d2c",
  "HCP2" = "#1b9e77",
  "HCP3" = "#a1d99b",
  "CHF1" = "#a63603",
  "CHF2" = "#d95f02",
  "CHF3" = "#fdae6b",
  "YSB1" = "#54278f",
  "YSB2" = "#7570b3",
  "YSB3" = "#bcbddc"
)

basin_colors <- c(
  "Herradura Clay Pan" = "#1b9e77",
  "Ckoirama Halite Field" = "#d95f02",
  "Yungay Station Basin" = "#7570b3"
)

distance_methods <- c(
  "wunifrac",
  "bray"
)


# ============================================================================ #
# Load final phyloseq object
# ============================================================================ #

ps <- readRDS(
  phyloseq_path
)

sample_data(ps)$Ephemeral_Basin <- factor(
  sample_data(ps)$Ephemeral_Basin,
  levels = basin_levels
)

sample_data(ps)$Sampling_Point <- factor(
  sample_data(ps)$Sampling_Point,
  levels = sampling_point_levels
)

ps_ra <- transform_sample_counts(
  ps,
  function(x) {

    if (
      sum(
        x,
        na.rm = TRUE
      ) == 0
    ) {
      x
    } else {
      100 *
        x /
        sum(
          x,
          na.rm = TRUE
        )
    }
  }
)


# ============================================================================ #
# Distance helpers
# ============================================================================ #

get_distance_label <- function(
    dist_method
) {

  dist_method <- match.arg(
    dist_method,
    c(
      "bray",
      "wunifrac"
    )
  )

  switch(
    dist_method,
    bray = "Bray-Curtis dissimilarity",
    wunifrac = "Weighted UniFrac distance"
  )
}


compute_selected_distance <- function(
    ps_obj,
    ps_ra_obj,
    dist_method
) {

  dist_method <- match.arg(
    dist_method,
    c(
      "bray",
      "wunifrac"
    )
  )

  switch(
    dist_method,
    bray = phyloseq::distance(
      ps_ra_obj,
      method = "bray"
    ),
    wunifrac = phyloseq::distance(
      ps_obj,
      method = "wunifrac"
    )
  )
}


format_p_plot <- function(
    p
) {

  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "< 0.001",
    p < 0.01 ~ "< 0.01",
    p < 0.05 ~ "< 0.05",
    TRUE ~ paste0(
      "= ",
      format(
        round(
          p,
          3
        ),
        nsmall = 3
      )
    )
  )
}


# ============================================================================ #
# Baseline similarity helper
# ============================================================================ #

compute_similarity_to_baseline <- function(
    ps_obj,
    ps_ra_obj,
    dist_method
) {

  metadata <- data.frame(
    sample_data(
      ps_obj
    ),
    check.names = FALSE
  ) |>
    tibble::rownames_to_column(
      "SampleID"
    ) |>
    dplyr::mutate(
      Date = as.Date(
        Date
      )
    ) |>
    dplyr::arrange(
      Ephemeral_Basin,
      Sampling_Point,
      Date
    )

  output_list <- list()

  for (
    basin in basin_levels
  ) {

    points_basin <- metadata |>
      dplyr::filter(
        Ephemeral_Basin == basin
      ) |>
      dplyr::pull(
        Sampling_Point
      ) |>
      unique()

    for (
      sampling_point in points_basin
    ) {

      meta_sub <- metadata |>
        dplyr::filter(
          Ephemeral_Basin == basin,
          Sampling_Point == sampling_point
        ) |>
        dplyr::arrange(
          Date
        )

      if (
        nrow(
          meta_sub
        ) < 2
      ) {
        next
      }

      samples_sub <- meta_sub$SampleID

      ps_sub <- prune_samples(
        samples_sub,
        ps_obj
      )

      ps_ra_sub <- prune_samples(
        samples_sub,
        ps_ra_obj
      )

      dist_mat <- as.matrix(
        compute_selected_distance(
          ps_obj = ps_sub,
          ps_ra_obj = ps_ra_sub,
          dist_method = dist_method
        )
      )

      baseline_sample <- meta_sub$SampleID[1]
      baseline_date <- meta_sub$Date[1]

      output_list[[
        paste(
          basin,
          sampling_point,
          "baseline",
          sep = "_"
        )
      ]] <- data.frame(
        Ephemeral_Basin = basin,
        Sampling_Point = sampling_point,
        SampleID = baseline_sample,
        Date = baseline_date,
        Days = 0,
        Dissimilarity = 0,
        stringsAsFactors = FALSE
      )

      if (
        nrow(
          meta_sub
        ) >= 2
      ) {

        for (
          i in 2:nrow(
            meta_sub
          )
        ) {

          sample_i <- meta_sub$SampleID[i]
          date_i <- meta_sub$Date[i]

          output_list[[
            paste(
              basin,
              sampling_point,
              sample_i,
              sep = "_"
            )
          ]] <- data.frame(
            Ephemeral_Basin = basin,
            Sampling_Point = sampling_point,
            SampleID = sample_i,
            Date = date_i,
            Days = as.numeric(
              date_i -
                baseline_date
            ),
            Dissimilarity = dist_mat[
              baseline_sample,
              sample_i
            ],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  dplyr::bind_rows(
    output_list
  ) |>
    dplyr::mutate(
      Similarity = (
        1 -
          Dissimilarity
      ) *
        100,
      Ephemeral_Basin = factor(
        Ephemeral_Basin,
        levels = basin_levels
      ),
      Sampling_Point = factor(
        Sampling_Point,
        levels = sampling_point_levels
      ),
      PointDays = Days +
        dplyr::case_when(
          Sampling_Point %in% c(
            "HCP1",
            "CHF1",
            "YSB1"
          ) ~ -1.5,
          Sampling_Point %in% c(
            "HCP2",
            "CHF2",
            "YSB2"
          ) ~ 0,
          Sampling_Point %in% c(
            "HCP3",
            "CHF3",
            "YSB3"
          ) ~ 1.5,
          TRUE ~ 0
        )
    )
}


# ============================================================================ #
# Baseline similarity plot helper
# ============================================================================ #

make_baseline_similarity_plot <- function(
    similarity_data
) {

  shaded_periods <- data.frame(
    Ephemeral_Basin = factor(
      basin_levels,
      levels = basin_levels
    ),
    xmin = c(
      0,
      0,
      0
    ),
    xmax = c(
      4,
      28,
      28
    )
  )

  water_end_lines <- shaded_periods |>
    dplyr::transmute(
      Ephemeral_Basin = Ephemeral_Basin,
      xintercept = xmax
    )

  ggplot(
    similarity_data,
    aes(
      x = Days,
      y = Similarity,
      colour = Sampling_Point
    )
  ) +
    geom_rect(
      data = shaded_periods,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "grey70",
      alpha = 0.5
    ) +
    geom_vline(
      data = water_end_lines,
      aes(
        xintercept = xintercept
      ),
      inherit.aes = FALSE,
      color = "grey15",
      linetype = "dotted",
      linewidth = 0.5
    ) +
    geom_smooth(
      aes(
        group = Sampling_Point
      ),
      method = "loess",
      se = FALSE,
      span = 0.14,
      linewidth = 0.8
    ) +
    geom_point(
      aes(
        x = PointDays,
        fill = Sampling_Point
      ),
      shape = 21,
      size = 2,
      stroke = 0.35,
      alpha = 0.9,
      color = "black"
    ) +
    facet_wrap(
      ~ Ephemeral_Basin,
      ncol = 1,
      scales = "fixed"
    ) +
    labs(
      x = "Time interval (days)",
      y = "Community similarity to baseline (T0) (%)"
    ) +
    theme_cowplot() +
    scale_x_continuous(
      limits = c(
        -3,
        460
      ),
      expand = c(
        0.01,
        0
      ),
      breaks = c(
        0,
        50,
        100,
        150,
        200,
        250,
        300,
        350,
        400,
        450
      )
    ) +
    scale_y_continuous(
      limits = c(
        0,
        100
      )
    ) +
    scale_color_manual(
      values = sampling_point_colors
    ) +
    scale_fill_manual(
      values = sampling_point_colors
    ) +
    theme(
      axis.title.x = element_text(
        size = 14
      ),
      axis.title.y = element_text(
        size = 14
      ),
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.5
      ),
      panel.grid.major.x = element_line(
        color = "grey15",
        linetype = "dotted"
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 12,
        face = "bold"
      )
    )
}


# ============================================================================ #
# Pairwise temporal-distance helper
# ============================================================================ #

compute_pairwise_temporal_data <- function(
    ps_obj,
    ps_ra_obj,
    dist_method
) {

  dist_obj <- compute_selected_distance(
    ps_obj = ps_obj,
    ps_ra_obj = ps_ra_obj,
    dist_method = dist_method
  )

  metadata <- data.frame(
    sample_data(
      ps_obj
    ),
    check.names = FALSE
  ) |>
    tibble::rownames_to_column(
      "SampleID"
    ) |>
    dplyr::mutate(
      Date = as.Date(
        Date
      )
    )

  dist_long <- as.data.frame(
    as.matrix(
      dist_obj
    )
  ) |>
    tibble::rownames_to_column(
      "Sample1"
    ) |>
    tidyr::pivot_longer(
      -Sample1,
      names_to = "Sample2",
      values_to = "Dissimilarity"
    ) |>
    dplyr::filter(
      Sample1 <
        Sample2
    ) |>
    dplyr::left_join(
      metadata |>
        dplyr::select(
          Sample1 = SampleID,
          Date1 = Date,
          Basin1 = Ephemeral_Basin
        ),
      by = "Sample1"
    ) |>
    dplyr::left_join(
      metadata |>
        dplyr::select(
          Sample2 = SampleID,
          Date2 = Date,
          Basin2 = Ephemeral_Basin
        ),
      by = "Sample2"
    ) |>
    dplyr::filter(
      Basin1 ==
        Basin2
    ) |>
    dplyr::mutate(
      Ephemeral_Basin = factor(
        Basin1,
        levels = basin_levels
      ),
      DiffDate = abs(
        as.numeric(
          Date2 -
            Date1
        )
      )
    )

  dist_long
}


# ============================================================================ #
# Mantel helper
# ============================================================================ #

compute_temporal_mantel_by_basin <- function(
    ps_obj,
    ps_ra_obj,
    dist_method,
    n_perm = 9999
) {

  metadata <- data.frame(
    sample_data(
      ps_obj
    ),
    check.names = FALSE
  ) |>
    tibble::rownames_to_column(
      "SampleID"
    ) |>
    dplyr::mutate(
      Date = as.Date(
        Date
      )
    ) |>
    dplyr::filter(
      !is.na(
        Date
      ),
      !is.na(
        Ephemeral_Basin
      )
    )

  output_list <- list()

  for (
    basin in basin_levels
  ) {

    meta_basin <- metadata |>
      dplyr::filter(
        Ephemeral_Basin == basin
      ) |>
      dplyr::arrange(
        Date
      )

    samples_basin <- meta_basin$SampleID

    if (
      length(
        samples_basin
      ) < 3
    ) {
      next
    }

    ps_sub <- prune_samples(
      samples_basin,
      ps_obj
    )

    ps_ra_sub <- prune_samples(
      samples_basin,
      ps_ra_obj
    )

    comm_dist <- compute_selected_distance(
      ps_obj = ps_sub,
      ps_ra_obj = ps_ra_sub,
      dist_method = dist_method
    )

    comm_labels <- attr(
      comm_dist,
      "Labels"
    )

    date_vec <- meta_basin$Date[
      match(
        comm_labels,
        meta_basin$SampleID
      )
    ]

    time_dist <- dist(
      as.numeric(
        date_vec
      ),
      method = "manhattan"
    )

    attr(
      time_dist,
      "Labels"
    ) <- comm_labels

    mantel_obj <- vegan::mantel(
      xdis = comm_dist,
      ydis = time_dist,
      method = "spearman",
      permutations = n_perm,
      na.rm = TRUE
    )

    output_list[[
      basin
    ]] <- data.frame(
      Ephemeral_Basin = basin,
      Mantel_method = "spearman",
      Mantel_r = as.numeric(
        mantel_obj$statistic
      ),
      P_value = mantel_obj$signif,
      Distance = dist_method,
      Permutations = n_perm,
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(
    output_list
  ) |>
    dplyr::mutate(
      Ephemeral_Basin = factor(
        Ephemeral_Basin,
        levels = basin_levels
      )
    )
}


# ============================================================================ #
# GLM slope helper for temporal-distance plots
# ============================================================================ #

extract_glm_slopes <- function(
    pairwise_data
) {

  pairwise_data |>
    dplyr::group_by(
      Ephemeral_Basin
    ) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(
        data,
        ~ glm(
          Dissimilarity ~ DiffDate,
          data = .x
        )
      ),
      tidy = purrr::map(
        model,
        broom::tidy
      )
    ) |>
    tidyr::unnest(
      tidy
    ) |>
    dplyr::filter(
      term == "DiffDate"
    ) |>
    dplyr::transmute(
      Ephemeral_Basin,
      Model = "GLM",
      Slope = estimate,
      P_value = p.value
    )
}


# ============================================================================ #
# Pairwise temporal-distance plot helper
# ============================================================================ #

make_temporal_distance_plot <- function(
    pairwise_data,
    mantel_results,
    glm_results,
    dist_method
) {

  distance_label <- get_distance_label(
    dist_method
  )

  strip_labels <- setNames(
    paste0(
      mantel_results$Ephemeral_Basin,
      "\nMantel test: r = ",
      signif(
        mantel_results$Mantel_r,
        3
      ),
      ", p-value ",
      format_p_plot(
        mantel_results$P_value
      )
    ),
    mantel_results$Ephemeral_Basin
  )

  glm_labels <- glm_results |>
    dplyr::mutate(
      label = paste0(
        "Slope = ",
        formatC(
          Slope,
          format = "e",
          digits = 2
        ),
        "\np-value ",
        format_p_plot(
          P_value
        )
      )
    )

  ggplot(
    pairwise_data,
    aes(
      x = DiffDate,
      y = Dissimilarity
    )
  ) +
    geom_point(
      shape = 21,
      size = 1.35,
      stroke = 0.25,
      fill = scales::alpha(
        "black",
        0.22
      ),
      colour = scales::alpha(
        "black",
        0.55
      ),
      position = position_jitter(
        width = 1.5,
        height = 0
      )
    ) +
    geom_smooth(
      aes(
        colour = Ephemeral_Basin
      ),
      method = "glm",
      se = TRUE,
      fill = "grey50",
      linewidth = 1.15,
      alpha = 0.4
    ) +
    facet_wrap(
      ~ Ephemeral_Basin,
      labeller = labeller(
        Ephemeral_Basin = strip_labels
      )
    ) +
    scale_color_manual(
      values = basin_colors,
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = seq(
        0,
        450,
        by = 100
      )
    ) +
    geom_text(
      data = glm_labels,
      aes(
        x = Inf,
        y = -Inf,
        label = label
      ),
      hjust = 1.03,
      vjust = -0.35,
      size = 3,
      colour = "black",
      lineheight = 1.05,
      inherit.aes = FALSE
    ) +
    labs(
      x = "Time difference (days)",
      y = distance_label
    ) +
    coord_cartesian(
      xlim = c(
        0,
        450
      ),
      clip = "off"
    ) +
    theme_cowplot() +
    theme(
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.7
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 11,
        face = "bold",
        margin = margin(
          b = 5
        )
      ),
      panel.spacing.x = grid::unit(
        0.35,
        "cm"
      ),
      axis.title = element_text(
        size = 15
      ),
      axis.text = element_text(
        size = 11
      ),
      plot.margin = margin(
        t = 6,
        r = 10,
        b = 8,
        l = 8
      )
    )
}


# ============================================================================ #
# Figure 5A-B | Similarity to initial community (T0)
# ============================================================================ #

baseline_results <- list()
baseline_plots <- list()

for (
  dist_method in distance_methods
) {

  similarity_data <- compute_similarity_to_baseline(
    ps_obj = ps,
    ps_ra_obj = ps_ra,
    dist_method = dist_method
  )

  baseline_results[[
    dist_method
  ]] <- similarity_data

  baseline_plots[[
    dist_method
  ]] <- make_baseline_similarity_plot(
    similarity_data
  )

  write.csv(
    similarity_data,
    file.path(
      table_dir,
      paste0(
        "similarity_to_T0_",
        dist_method,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

ggsave(
  file.path(
    figure_dir,
    "Figure5A_similarity_T0_wunifrac.tiff"
  ),
  plot = baseline_plots$wunifrac,
  width = 8,
  height = 6,
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure5A_similarity_T0_wunifrac.svg"
  ),
  plot = baseline_plots$wunifrac,
  width = 8,
  height = 6,
  device = "svg"
)

ggsave(
  file.path(
    figure_dir,
    "Figure5B_similarity_T0_bray.tiff"
  ),
  plot = baseline_plots$bray,
  width = 8,
  height = 6,
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure5B_similarity_T0_bray.svg"
  ),
  plot = baseline_plots$bray,
  width = 8,
  height = 6,
  device = "svg"
)


# ============================================================================ #
# Figure 5C-D | Pairwise temporal distance-decay
# ============================================================================ #

pairwise_results <- list()
mantel_results_all <- list()
glm_results_all <- list()
temporal_plots <- list()

for (
  dist_method in distance_methods
) {

  pairwise_data <- compute_pairwise_temporal_data(
    ps_obj = ps,
    ps_ra_obj = ps_ra,
    dist_method = dist_method
  )

  mantel_results <- compute_temporal_mantel_by_basin(
    ps_obj = ps,
    ps_ra_obj = ps_ra,
    dist_method = dist_method,
    n_perm = 9999
  )

  glm_results <- extract_glm_slopes(
    pairwise_data
  )

  temporal_plot <- make_temporal_distance_plot(
    pairwise_data = pairwise_data,
    mantel_results = mantel_results,
    glm_results = glm_results,
    dist_method = dist_method
  )

  pairwise_results[[
    dist_method
  ]] <- pairwise_data

  mantel_results_all[[
    dist_method
  ]] <- mantel_results

  glm_results_all[[
    dist_method
  ]] <- glm_results

  temporal_plots[[
    dist_method
  ]] <- temporal_plot

  write.csv(
    pairwise_data,
    file.path(
      table_dir,
      paste0(
        "pairwise_temporal_distance_",
        dist_method,
        ".csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    mantel_results,
    file.path(
      table_dir,
      paste0(
        "mantel_temporal_",
        dist_method,
        "_spearman_9999.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    glm_results,
    file.path(
      table_dir,
      paste0(
        "glm_temporal_",
        dist_method,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

ggsave(
  file.path(
    figure_dir,
    "Figure5C_temporal_distance_wunifrac.tiff"
  ),
  plot = temporal_plots$wunifrac,
  width = 7,
  height = 3,
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure5C_temporal_distance_wunifrac.svg"
  ),
  plot = temporal_plots$wunifrac,
  width = 7,
  height = 3,
  device = "svg"
)

ggsave(
  file.path(
    figure_dir,
    "Figure5D_temporal_distance_bray.tiff"
  ),
  plot = temporal_plots$bray,
  width = 7,
  height = 3,
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure5D_temporal_distance_bray.svg"
  ),
  plot = temporal_plots$bray,
  width = 7,
  height = 3,
  device = "svg"
)


# ============================================================================ #
# Combined statistics workbook
# ============================================================================ #

mantel_export <- dplyr::bind_rows(
  lapply(
    names(
      mantel_results_all
    ),
    function(metric) {

      mantel_results_all[[
        metric
      ]] |>
        dplyr::mutate(
          Distance = metric
        )
    }
  )
)

glm_export <- dplyr::bind_rows(
  lapply(
    names(
      glm_results_all
    ),
    function(metric) {

      glm_results_all[[
        metric
      ]] |>
        dplyr::mutate(
          Distance = metric
        )
    }
  )
)

openxlsx::write.xlsx(
  list(
    Mantel_tests = mantel_export,
    GLM_slopes = glm_export
  ),
  file = file.path(
    table_dir,
    "Figure5_temporal_turnover_statistics.xlsx"
  )
)


# ============================================================================ #
# Final summary
# ============================================================================ #

cat(
  "\nTemporal beta-diversity analysis completed.\n"
)

cat(
  "Samples: ",
  nsamples(
    ps
  ),
  "\n",
  sep = ""
)

cat(
  "ASVs: ",
  ntaxa(
    ps
  ),
  "\n",
  sep = ""
)

cat(
  "Distance metrics: weighted UniFrac and Bray-Curtis\n"
)

cat(
  "Mantel permutations: 9999\n"
)
