# ============================================================================ #
# 16_Mantel_community_environment.R
# ============================================================================ #
#
# Purpose:
# Reproduce the final basin-specific Mantel analyses supporting
# Supplementary Figure 7 and Supplementary Table 5.
#
# Analysis blocks:
#
#   A. Core physicochemical Mantel analysis
#      Wc, Salinity, pH, ORP
#      3 basins × 2 community metrics × 4 variables = 24 tests
#
#   B. Geochemical Mantel analysis
#      TOC, Na, Mo, Sr, Se, Pb, K, Li, Fe
#      3 basins × 2 community metrics × 9 variables = 54 tests
#
# Community metrics:
#   - Bray-Curtis
#   - Weighted UniFrac
#
# Multiple-testing correction:
#   - global BH across the 24 core tests
#   - global BH across the 54 geochemical tests
#   - basin-wise BH retained only as a sensitivity reference
#
# Robustness analyses:
#   - temporally restricted cyclic permutations for the 24 core tests
#   - leave-one-time-point-out (LOTO) directional stability for core tests
#   - LOTO directional stability for geochemical tests
#
# Important:
# - Uses the NON-RAREFIED final phyloseq object.
# - Core analysis works at the temporal level. For each pair of time points,
#   community distances are averaged across matching permanent sampling
#   points within each basin. Environmental values are basin-time medians.
# - Wc and Salinity are log1p-transformed before standardization in the core
#   environmental distance matrices.
# - The geochemical analysis uses the complete 29-sample subset for all nine
#   selected variables. Concentrations are log1p-transformed and standardized.
# - Mantel tests use Spearman correlation and 9,999 permutations.
#
# Input:
# - output_phyloseq_tree/ps_mdecon_v2_FTree.rds
#
# Main outputs:
# - Mantel_core_results.tsv
# - Mantel_geochemistry_results.tsv
# - Mantel_core_final.*
# - Mantel_geochemistry_matrix_links.*
# - Table_Supplementary5_Mantel.xlsx
# - temporal sensitivity / LOTO tables
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
library(ggplot2)
library(cowplot)
library(patchwork)
library(openxlsx)


# ============================================================================ #
# Paths
# ============================================================================ #

phyloseq_path <- file.path(
  "output_phyloseq_tree",
  "ps_mdecon_v2_FTree.rds"
)

output_dir <- "output_mantel_community_environment"

figures_dir <- file.path(
  output_dir,
  "figures"
)

tables_dir <- file.path(
  output_dir,
  "tables"
)

objects_dir <- file.path(
  output_dir,
  "objects"
)

dir.create(
  figures_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  tables_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  objects_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================================ #
# Load final phyloseq object
# ============================================================================ #

ps <- readRDS(
  phyloseq_path
)

ps <- prune_taxa(
  taxa_sums(
    ps
  ) >
    0,
  ps
)

meta_mantel <- data.frame(
  sample_data(
    ps
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta_mantel$SampleID <- rownames(
  meta_mantel
)


# ============================================================================ #
# Shared settings
# ============================================================================ #

basin_levels <- c(
  "HCP",
  "CHF",
  "YSB"
)

basin_names <- c(
  "HCP" = "Herradura Clay Pan",
  "CHF" = "Ckoirama Halite Field",
  "YSB" = "Yungay Station Basin"
)

community_metrics <- c(
  "Bray-Curtis",
  "Weighted UniFrac"
)

vars_core <- c(
  "Wc",
  "Salinity",
  "pH",
  "ORP"
)

vars_geo_29 <- c(
  "TOC",
  "Na",
  "Mo",
  "Sr",
  "Se",
  "Pb",
  "K",
  "Li",
  "Fe"
)

N_PERM <- 9999


# ============================================================================ #
# Core analysis | Temporal community distance
# ============================================================================ #

build_temporal_comm_dist <- function(
    ps_obj,
    meta_df,
    basin_name,
    metric
) {

  meta_b <- meta_df |>
    filter(
      as.character(
        Basin
      ) ==
        basin_name
    ) |>
    arrange(
      Time_Point,
      Sampling_Point
    )

  ps_b <- prune_samples(
    sample_names(
      ps_obj
    ) %in%
      meta_b$SampleID,
    ps_obj
  )

  ps_b <- prune_taxa(
    taxa_sums(
      ps_b
    ) >
      0,
    ps_b
  )

  ps_b_rel <- transform_sample_counts(
    ps_b,
    function(x) {
      x /
        sum(
          x
        )
    }
  )

  if (
    metric ==
      "Bray-Curtis"
  ) {

    d_samples <- as.matrix(
      phyloseq::distance(
        ps_b_rel,
        method = "bray"
      )
    )

  } else if (
    metric ==
      "Weighted UniFrac"
  ) {

    d_samples <- as.matrix(
      phyloseq::UniFrac(
        ps_b_rel,
        weighted = TRUE,
        normalized = TRUE,
        parallel = FALSE,
        fast = TRUE
      )
    )

  } else {

    stop(
      "Unknown community metric."
    )
  }

  times <- sort(
    unique(
      meta_b$Time_Point
    )
  )

  D_time <- matrix(
    0,
    nrow = length(
      times
    ),
    ncol = length(
      times
    ),
    dimnames = list(
      as.character(
        times
      ),
      as.character(
        times
      )
    )
  )

  sampling_points <- unique(
    as.character(
      meta_b$Sampling_Point
    )
  )

  for (
    i in seq_len(
      length(
        times
      ) -
        1
    )
  ) {

    for (
      j in (
        i +
          1
      ):length(
        times
      )
    ) {

      vals <- numeric(
        0
      )

      for (
        pt in
          sampling_points
      ) {

        id_i <- meta_b$SampleID[
          meta_b$Time_Point ==
            times[
              i
            ] &
            as.character(
              meta_b$Sampling_Point
            ) ==
            pt
        ]

        id_j <- meta_b$SampleID[
          meta_b$Time_Point ==
            times[
              j
            ] &
            as.character(
              meta_b$Sampling_Point
            ) ==
            pt
        ]

        if (
          length(
            id_i
          ) ==
            1 &&
            length(
              id_j
            ) ==
            1 &&
            id_i %in%
            rownames(
              d_samples
            ) &&
            id_j %in%
            rownames(
              d_samples
            )
        ) {

          vals <- c(
            vals,
            d_samples[
              id_i,
              id_j
            ]
          )
        }
      }

      if (
        length(
          vals
        ) ==
          0
      ) {

        D_time[
          i,
          j
        ] <- NA_real_

        D_time[
          j,
          i
        ] <- NA_real_

      } else {

        D_time[
          i,
          j
        ] <- mean(
          vals
        )

        D_time[
          j,
          i
        ] <- mean(
          vals
        )
      }
    }
  }

  if (
    anyNA(
      D_time
    )
  ) {
    stop(
      "Missing temporal microbial distances in ",
      basin_name,
      " / ",
      metric,
      "."
    )
  }

  as.dist(
    D_time
  )
}


# ============================================================================ #
# Core analysis | Temporal environmental distance
# ============================================================================ #

build_temporal_env_dist <- function(
    meta_df,
    basin_name,
    variable
) {

  env_b <- meta_df |>
    filter(
      as.character(
        Basin
      ) ==
        basin_name
    ) |>
    group_by(
      Time_Point
    ) |>
    summarise(
      value = median(
        .data[
          [
            variable
          ]
        ],
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    arrange(
      Time_Point
    )

  if (
    variable %in%
      c(
        "Wc",
        "Salinity"
      )
  ) {
    env_b$value <- log1p(
      env_b$value
    )
  }

  env_b$value <- as.numeric(
    scale(
      env_b$value
    )
  )

  env_matrix <- matrix(
    env_b$value,
    ncol = 1,
    dimnames = list(
      as.character(
        env_b$Time_Point
      ),
      variable
    )
  )

  dist(
    env_matrix
  )
}


# ============================================================================ #
# Core analysis | Cache community-distance objects
# ============================================================================ #

core_dist_file <- file.path(
  objects_dir,
  "Mantel_core_community_distances.rds"
)

if (
  file.exists(
    core_dist_file
  )
) {

  core_comm_dist <- readRDS(
    core_dist_file
  )

} else {

  core_comm_dist <- list()

  for (
    basin_i in
      basin_levels
  ) {

    for (
      metric_i in
        community_metrics
    ) {

      key_i <- paste(
        basin_i,
        metric_i,
        sep = "__"
      )

      core_comm_dist[
        [
          key_i
        ]
      ] <- build_temporal_comm_dist(
        ps_obj = ps,
        meta_df = meta_mantel,
        basin_name = basin_i,
        metric = metric_i
      )
    }
  }

  saveRDS(
    core_comm_dist,
    core_dist_file
  )
}


# ============================================================================ #
# Core analysis | Mantel tests
# ============================================================================ #

set.seed(
  2026
)

mantel_results_list <- list()
counter <- 1

for (
  basin_i in
    basin_levels
) {

  for (
    metric_i in
      community_metrics
  ) {

    key_i <- paste(
      basin_i,
      metric_i,
      sep = "__"
    )

    d_comm <- core_comm_dist[
      [
        key_i
      ]
    ]

    for (
      env_i in
        vars_core
    ) {

      d_env <- build_temporal_env_dist(
        meta_df = meta_mantel,
        basin_name = basin_i,
        variable = env_i
      )

      if (
        !identical(
          attr(
            d_comm,
            "Labels"
          ),
          attr(
            d_env,
            "Labels"
          )
        )
      ) {
        stop(
          "Temporal labels do not match for ",
          basin_i,
          " / ",
          metric_i,
          " / ",
          env_i,
          "."
        )
      }

      mantel_i <- vegan::mantel(
        d_comm,
        d_env,
        method = "spearman",
        permutations = N_PERM
      )

      mantel_results_list[
        [
          counter
        ]
      ] <- data.frame(
        Basin = basin_i,
        Community_metric = metric_i,
        Environment = env_i,
        Mantel_r = unname(
          mantel_i$statistic
        ),
        p_value =
          mantel_i$signif,
        stringsAsFactors = FALSE
      )

      counter <- counter +
        1
    }
  }
}

mantel_results <- bind_rows(
  mantel_results_list
)

mantel_results$p_adj_BH_global <- p.adjust(
  mantel_results$p_value,
  method = "BH"
)

mantel_results <- mantel_results |>
  group_by(
    Basin
  ) |>
  mutate(
    p_adj_BH_basin = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  ungroup() |>
  mutate(
    Significance = case_when(
      p_adj_BH_global <
        0.001 ~
        "FDR < 0.001",
      p_adj_BH_global <
        0.01 ~
        "FDR < 0.01",
      p_adj_BH_global <
        0.05 ~
        "FDR < 0.05",
      TRUE ~
        "ns"
    )
  )

write.table(
  mantel_results,
  file.path(
    tables_dir,
    "Mantel_core_results.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================================ #
# Core analysis | Environmental Spearman correlations
# ============================================================================ #

env_time <- meta_mantel |>
  group_by(
    Basin,
    Time_Point
  ) |>
  summarise(
    Wc = median(
      Wc,
      na.rm = TRUE
    ),
    Salinity = median(
      Salinity,
      na.rm = TRUE
    ),
    pH = median(
      pH,
      na.rm = TRUE
    ),
    ORP = median(
      ORP,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

env_cor_list <- list()

for (
  basin_i in
    basin_levels
) {

  x <- env_time |>
    filter(
      as.character(
        Basin
      ) ==
        basin_i
    ) |>
    select(
      all_of(
        vars_core
      )
    )

  cor_i <- cor(
    x,
    method = "spearman",
    use = "pairwise.complete.obs"
  )

  cor_df <- as.data.frame(
    as.table(
      cor_i
    )
  )

  colnames(
    cor_df
  ) <- c(
    "Variable1",
    "Variable2",
    "Spearman_r"
  )

  cor_df$Basin <-
    basin_i

  env_cor_list[
    [
      basin_i
    ]
  ] <- cor_df
}

env_cor_results <- bind_rows(
  env_cor_list
)

env_cor_results$Variable1 <- factor(
  env_cor_results$Variable1,
  levels = vars_core
)

env_cor_results$Variable2 <- factor(
  env_cor_results$Variable2,
  levels = vars_core
)

env_cor_results <- env_cor_results |>
  mutate(
    i = as.numeric(
      Variable1
    ),
    j = as.numeric(
      Variable2
    )
  ) |>
  filter(
    i <
      j
  ) |>
  select(
    Basin,
    Variable1,
    Variable2,
    Spearman_r
  )

write.table(
  env_cor_results,
  file.path(
    tables_dir,
    "Mantel_environmental_correlations.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================================ #
# Core analysis | Figure component
# ============================================================================ #

env_order <- vars_core

cor_plot_df <- env_cor_results |>
  mutate(
    Variable1 = factor(
      Variable1,
      levels = env_order
    ),
    Variable2 = factor(
      Variable2,
      levels = env_order
    ),
    i = as.numeric(
      Variable1
    ),
    j = as.numeric(
      Variable2
    )
  ) |>
  mutate(
    x = j,
    y = 5 -
      i,
    cell_size =
      0.12 +
      0.70 *
      abs(
        Spearman_r
      )
  )

diag_df <- expand.grid(
  Basin = basin_levels,
  Environment = env_order,
  stringsAsFactors = FALSE
) |>
  mutate(
    x = match(
      Environment,
      env_order
    ),
    y = 5 -
      x
  )

top_env_labels <- expand.grid(
  Basin = basin_levels,
  Environment = env_order,
  stringsAsFactors = FALSE
) |>
  mutate(
    x = match(
      Environment,
      env_order
    ),
    y = 4.95
  )

right_env_labels <- expand.grid(
  Basin = basin_levels,
  Environment = env_order,
  stringsAsFactors = FALSE
) |>
  mutate(
    x = 4.55,
    y = 5 -
      match(
        Environment,
        env_order
      )
  )

bio_nodes <- expand.grid(
  Basin = basin_levels,
  Community_metric =
    community_metrics,
  stringsAsFactors = FALSE
) |>
  mutate(
    x = -1.25,
    y = ifelse(
      Community_metric ==
        "Bray-Curtis",
      3.15,
      1.85
    )
  )

mantel_edges <- mantel_results |>
  mutate(
    x = -1.25,
    y = ifelse(
      Community_metric ==
        "Bray-Curtis",
      3.15,
      1.85
    ),
    xend = match(
      Environment,
      env_order
    ),
    yend = 5 -
      xend,
    Plot_significance = ifelse(
      p_adj_BH_global <
        0.05,
      "BH-FDR < 0.05",
      "ns"
    )
  )

for (
  object_name in
    c(
      "cor_plot_df",
      "diag_df",
      "top_env_labels",
      "right_env_labels",
      "bio_nodes",
      "mantel_edges"
    )
) {

  x <- get(
    object_name
  )

  x$Basin <- factor(
    x$Basin,
    levels = basin_levels
  )

  assign(
    object_name,
    x
  )
}

p_mantel_core <- ggplot() +
  geom_curve(
    data = filter(
      mantel_edges,
      Community_metric ==
        "Bray-Curtis"
    ),
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      colour =
        Plot_significance,
      linewidth =
        abs(
          Mantel_r
        )
    ),
    curvature = 0.08,
    lineend = "round",
    alpha = 0.90
  ) +
  geom_curve(
    data = filter(
      mantel_edges,
      Community_metric ==
        "Weighted UniFrac"
    ),
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      colour =
        Plot_significance,
      linewidth =
        abs(
          Mantel_r
        )
    ),
    curvature = -0.08,
    lineend = "round",
    alpha = 0.90
  ) +
  geom_tile(
    data = cor_plot_df,
    aes(
      x = x,
      y = y
    ),
    width = 0.88,
    height = 0.88,
    fill = NA,
    colour = "#bdbdbd",
    linewidth = 0.35
  ) +
  geom_tile(
    data = cor_plot_df,
    aes(
      x = x,
      y = y,
      width =
        cell_size,
      height =
        cell_size,
      fill =
        Spearman_r
    ),
    colour = "#7f7f7f",
    linewidth = 0.30
  ) +
  geom_point(
    data = diag_df,
    aes(
      x = x,
      y = y
    ),
    shape = 21,
    size = 2.8,
    fill = "white",
    colour = "black",
    stroke = 0.7
  ) +
  geom_text(
    data = top_env_labels,
    aes(
      x = x,
      y = y,
      label =
        Environment
    ),
    angle = 45,
    hjust = 0,
    vjust = 0.5,
    size = 3.1
  ) +
  geom_text(
    data = right_env_labels,
    aes(
      x = x,
      y = y,
      label =
        Environment
    ),
    hjust = 0,
    size = 3.1
  ) +
  geom_point(
    data = bio_nodes,
    aes(
      x = x,
      y = y,
      shape =
        Community_metric
    ),
    size = 4,
    fill = "white",
    colour = "black",
    stroke = 0.8,
    show.legend = FALSE
  ) +
  geom_text(
    data = bio_nodes,
    aes(
      x = x -
        0.18,
      y = y,
      label =
        Community_metric
    ),
    hjust = 1,
    size = 3.1
  ) +
  scale_fill_gradient2(
    low = "#8c510a",
    mid = "white",
    high = "#01665e",
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    name = "Spearman \u03c1"
  ) +
  scale_colour_manual(
    values = c(
      "BH-FDR < 0.05" =
        "#238b45",
      "ns" =
        "#dddddd"
    ),
    name = "Mantel BH-FDR"
  ) +
  scale_linewidth_continuous(
    range = c(
      0.25,
      2.25
    ),
    limits = c(
      0,
      0.55
    ),
    breaks = c(
      0.1,
      0.2,
      0.3,
      0.4,
      0.5
    ),
    name = "Mantel |r|"
  ) +
  scale_shape_manual(
    values = c(
      "Bray-Curtis" = 21,
      "Weighted UniFrac" = 22
    )
  ) +
  facet_wrap(
    ~ Basin,
    nrow = 1,
    labeller =
      as_labeller(
        basin_names
      )
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(
      -3.0,
      5.55
    ),
    ylim = c(
      0.30,
      5.35
    ),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_void(
    base_size = 11
  ) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 12,
      margin = margin(
        b = 10
      )
    ),
    panel.spacing = grid::unit(
      1,
      "lines"
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    legend.text = element_text(
      size = 9
    ),
    panel.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    plot.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    legend.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    plot.margin = margin(
      25,
      30,
      10,
      55
    )
  )

ggsave(
  file.path(
    figures_dir,
    "Mantel_core_final.png"
  ),
  p_mantel_core,
  width = 13,
  height = 5.2,
  dpi = 300,
  bg = "transparent"
)

ggsave(
  file.path(
    figures_dir,
    "Mantel_core_final.svg"
  ),
  p_mantel_core,
  width = 13,
  height = 5.2,
  bg = "transparent"
)


# ============================================================================ #
# Geochemistry | Complete 29-sample subset
# ============================================================================ #

meta_geo_29 <- meta_mantel |>
  filter(
    if_all(
      all_of(
        vars_geo_29
      ),
      ~ !is.na(
        .x
      )
    )
  ) |>
  arrange(
    Basin,
    Time_Point
  )

ps_geo_29 <- prune_samples(
  sample_names(
    ps
  ) %in%
    meta_geo_29$SampleID,
  ps
)

ps_geo_29 <- prune_taxa(
  taxa_sums(
    ps_geo_29
  ) >
    0,
  ps_geo_29
)


# ============================================================================ #
# Geochemistry | Community distance
# ============================================================================ #

build_geo_comm_dist <- function(
    ps_obj,
    meta_df,
    basin_name,
    metric
) {

  meta_b <- meta_df |>
    filter(
      as.character(
        Basin
      ) ==
        basin_name
    ) |>
    arrange(
      Time_Point
    )

  ps_b <- prune_samples(
    sample_names(
      ps_obj
    ) %in%
      meta_b$SampleID,
    ps_obj
  )

  ps_b <- prune_taxa(
    taxa_sums(
      ps_b
    ) >
      0,
    ps_b
  )

  ps_b <- prune_samples(
    meta_b$SampleID,
    ps_b
  )

  ps_b_rel <- transform_sample_counts(
    ps_b,
    function(x) {
      x /
        sum(
          x
        )
    }
  )

  if (
    metric ==
      "Bray-Curtis"
  ) {

    d_comm <- phyloseq::distance(
      ps_b_rel,
      method = "bray"
    )

  } else if (
    metric ==
      "Weighted UniFrac"
  ) {

    d_comm <- phyloseq::UniFrac(
      ps_b_rel,
      weighted = TRUE,
      normalized = TRUE,
      parallel = FALSE,
      fast = TRUE
    )

  } else {

    stop(
      "Unknown community metric."
    )
  }

  d_comm
}


# ============================================================================ #
# Geochemistry | Environmental distance
# ============================================================================ #

build_geo_env_dist <- function(
    meta_df,
    basin_name,
    variable
) {

  env_b <- meta_df |>
    filter(
      as.character(
        Basin
      ) ==
        basin_name
    ) |>
    arrange(
      Time_Point
    )

  x <- env_b[
    [
      variable
    ]
  ]

  if (
    any(
      x <
        0,
      na.rm = TRUE
    )
  ) {
    stop(
      "Negative values detected in ",
      variable,
      "."
    )
  }

  x <- log1p(
    x
  )

  x <- as.numeric(
    scale(
      x
    )
  )

  env_matrix <- matrix(
    x,
    ncol = 1,
    dimnames = list(
      env_b$SampleID,
      variable
    )
  )

  dist(
    env_matrix
  )
}


# ============================================================================ #
# Geochemistry | Mantel tests
# ============================================================================ #

set.seed(
  20260826
)

mantel_geo_list <- list()
counter <- 1

for (
  basin_i in
    basin_levels
) {

  for (
    metric_i in
      community_metrics
  ) {

    d_comm <- build_geo_comm_dist(
      ps_obj = ps_geo_29,
      meta_df = meta_geo_29,
      basin_name = basin_i,
      metric = metric_i
    )

    for (
      geo_i in
        vars_geo_29
    ) {

      d_geo <- build_geo_env_dist(
        meta_df = meta_geo_29,
        basin_name = basin_i,
        variable = geo_i
      )

      if (
        !identical(
          attr(
            d_comm,
            "Labels"
          ),
          attr(
            d_geo,
            "Labels"
          )
        )
      ) {
        stop(
          "Sample order mismatch for ",
          basin_i,
          " / ",
          metric_i,
          " / ",
          geo_i,
          "."
        )
      }

      mantel_i <- vegan::mantel(
        d_comm,
        d_geo,
        method = "spearman",
        permutations = N_PERM
      )

      mantel_geo_list[
        [
          counter
        ]
      ] <- data.frame(
        Basin = basin_i,
        Community_metric =
          metric_i,
        Geochemistry =
          geo_i,
        n_samples = nrow(
          filter(
            meta_geo_29,
            as.character(
              Basin
            ) ==
              basin_i
          )
        ),
        Mantel_r = unname(
          mantel_i$statistic
        ),
        p_value =
          mantel_i$signif,
        stringsAsFactors = FALSE
      )

      counter <- counter +
        1
    }
  }
}

mantel_geo_results <- bind_rows(
  mantel_geo_list
)

mantel_geo_results$p_adj_BH_global <- p.adjust(
  mantel_geo_results$p_value,
  method = "BH"
)

mantel_geo_results <- mantel_geo_results |>
  group_by(
    Basin
  ) |>
  mutate(
    p_adj_BH_basin = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  ungroup() |>
  mutate(
    Significance = ifelse(
      p_adj_BH_global <
        0.05,
      "BH-FDR < 0.05",
      "ns"
    )
  )

write.table(
  mantel_geo_results,
  file.path(
    tables_dir,
    "Mantel_geochemistry_results.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================================ #
# Geochemistry | Spearman correlation matrices
# ============================================================================ #

geo_cor_list <- list()

for (
  basin_i in
    basin_levels
) {

  geo_b <- meta_geo_29 |>
    filter(
      as.character(
        Basin
      ) ==
        basin_i
    ) |>
    select(
      all_of(
        vars_geo_29
      )
    )

  cor_i <- cor(
    geo_b,
    method = "spearman",
    use = "pairwise.complete.obs"
  )

  cor_df <- as.data.frame(
    as.table(
      cor_i
    )
  )

  colnames(
    cor_df
  ) <- c(
    "Variable1",
    "Variable2",
    "Spearman_r"
  )

  cor_df$Basin <-
    basin_i

  geo_cor_list[
    [
      basin_i
    ]
  ] <- cor_df
}

geo_cor_results <- bind_rows(
  geo_cor_list
)

write.table(
  geo_cor_results,
  file.path(
    tables_dir,
    "Mantel_geochemistry_Spearman_correlations.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================================ #
# Geochemistry | Figure component
# ============================================================================ #

geo_order <- vars_geo_29

geo_cor_plot <- geo_cor_results |>
  mutate(
    Variable1 = factor(
      Variable1,
      levels = geo_order
    ),
    Variable2 = factor(
      Variable2,
      levels = geo_order
    ),
    i = as.numeric(
      Variable1
    ),
    j = as.numeric(
      Variable2
    )
  ) |>
  filter(
    i <
      j
  ) |>
  mutate(
    x = j,
    y = 10 -
      i,
    cell_size =
      0.12 +
      0.70 *
      abs(
        Spearman_r
      )
  )

diag_geo <- expand.grid(
  Basin = basin_levels,
  Geochemistry = geo_order,
  stringsAsFactors = FALSE
) |>
  mutate(
    x = match(
      Geochemistry,
      geo_order
    ),
    y = 10 -
      x
  )

top_geo_labels <- expand.grid(
  Basin = basin_levels,
  Geochemistry = geo_order,
  stringsAsFactors = FALSE
) |>
  mutate(
    x = match(
      Geochemistry,
      geo_order
    ),
    y = 9.95
  )

right_geo_labels <- expand.grid(
  Basin = basin_levels,
  Geochemistry = geo_order,
  stringsAsFactors = FALSE
) |>
  mutate(
    x = 9.55,
    y = 10 -
      match(
        Geochemistry,
        geo_order
      )
  )

bio_geo_nodes <- expand.grid(
  Basin = basin_levels,
  Community_metric =
    community_metrics,
  stringsAsFactors = FALSE
) |>
  mutate(
    x = -1.55,
    y = ifelse(
      Community_metric ==
        "Bray-Curtis",
      7.2,
      4.8
    )
  )

geo_edges <- mantel_geo_results |>
  mutate(
    x = -1.55,
    y = ifelse(
      Community_metric ==
        "Bray-Curtis",
      7.2,
      4.8
    ),
    xend = match(
      Geochemistry,
      geo_order
    ),
    yend = 10 -
      xend
  )

for (
  object_name in
    c(
      "geo_cor_plot",
      "diag_geo",
      "top_geo_labels",
      "right_geo_labels",
      "bio_geo_nodes",
      "geo_edges"
    )
) {

  x <- get(
    object_name
  )

  x$Basin <- factor(
    x$Basin,
    levels = basin_levels
  )

  assign(
    object_name,
    x
  )
}

p_mantel_geo <- ggplot() +
  geom_curve(
    data = filter(
      geo_edges,
      Community_metric ==
        "Bray-Curtis"
    ),
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      colour =
        Significance,
      linewidth =
        abs(
          Mantel_r
        )
    ),
    curvature = 0.08,
    lineend = "round",
    alpha = 0.90
  ) +
  geom_curve(
    data = filter(
      geo_edges,
      Community_metric ==
        "Weighted UniFrac"
    ),
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      colour =
        Significance,
      linewidth =
        abs(
          Mantel_r
        )
    ),
    curvature = -0.08,
    lineend = "round",
    alpha = 0.90
  ) +
  geom_tile(
    data = geo_cor_plot,
    aes(
      x = x,
      y = y
    ),
    width = 0.88,
    height = 0.88,
    fill = NA,
    colour = "#bdbdbd",
    linewidth = 0.35
  ) +
  geom_tile(
    data = geo_cor_plot,
    aes(
      x = x,
      y = y,
      width =
        cell_size,
      height =
        cell_size,
      fill =
        Spearman_r
    ),
    colour = "#7f7f7f",
    linewidth = 0.30
  ) +
  geom_point(
    data = diag_geo,
    aes(
      x = x,
      y = y
    ),
    shape = 21,
    size = 2.8,
    fill = "white",
    colour = "black",
    stroke = 0.7
  ) +
  geom_text(
    data = top_geo_labels,
    aes(
      x = x,
      y = y,
      label =
        Geochemistry
    ),
    angle = 45,
    hjust = 0,
    vjust = 0.5,
    size = 3.1
  ) +
  geom_text(
    data = right_geo_labels,
    aes(
      x = x,
      y = y,
      label =
        Geochemistry
    ),
    hjust = 0,
    size = 3.1
  ) +
  geom_point(
    data = bio_geo_nodes,
    aes(
      x = x,
      y = y,
      shape =
        Community_metric
    ),
    size = 4,
    fill = "white",
    colour = "black",
    stroke = 0.8,
    show.legend = FALSE
  ) +
  geom_text(
    data = bio_geo_nodes,
    aes(
      x = x -
        0.18,
      y = y,
      label =
        Community_metric
    ),
    hjust = 1,
    size = 3.1
  ) +
  scale_fill_gradient2(
    low = "#8c510a",
    mid = "white",
    high = "#01665e",
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    name = "Spearman \u03c1"
  ) +
  scale_colour_manual(
    values = c(
      "BH-FDR < 0.05" =
        "#238b45",
      "ns" =
        "#dddddd"
    ),
    name = "Mantel BH-FDR"
  ) +
  scale_linewidth_continuous(
    range = c(
      0.25,
      2.25
    ),
    limits = c(
      0,
      0.55
    ),
    breaks = c(
      0.1,
      0.2,
      0.3,
      0.4,
      0.5
    ),
    name = "Mantel |r|"
  ) +
  scale_shape_manual(
    values = c(
      "Bray-Curtis" = 21,
      "Weighted UniFrac" = 22
    )
  ) +
  facet_wrap(
    ~ Basin,
    nrow = 1,
    labeller =
      as_labeller(
        basin_names
      )
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(
      -3.3,
      11.2
    ),
    ylim = c(
      0.3,
      10.6
    ),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_void(
    base_size = 11
  ) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 12,
      margin = margin(
        b = 10
      )
    ),
    panel.spacing = grid::unit(
      1,
      "lines"
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    legend.text = element_text(
      size = 9
    ),
    panel.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    plot.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    legend.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    plot.margin = margin(
      25,
      30,
      10,
      55
    )
  )

ggsave(
  file.path(
    figures_dir,
    "Mantel_geochemistry_matrix_links.png"
  ),
  p_mantel_geo,
  width = 16,
  height = 7.5,
  dpi = 300,
  bg = "transparent"
)

ggsave(
  file.path(
    figures_dir,
    "Mantel_geochemistry_matrix_links.svg"
  ),
  p_mantel_geo,
  width = 16,
  height = 7.5,
  bg = "transparent"
)


# ============================================================================ #
# Core robustness | Restricted cyclic permutations
# ============================================================================ #

make_cyclic_permutations <- function(
    n,
    mirror = TRUE
) {

  base_order <- seq_len(
    n
  )

  shift_vector <- function(
      x,
      k
  ) {

    if (
      k ==
        0
    ) {
      return(
        x
      )
    }

    c(
      x[
        (
          k +
            1
        ):length(
          x
        )
      ],
      x[
        seq_len(
          k
        )
      ]
    )
  }

  perms <- lapply(
    seq_len(
      n -
        1
    ),
    function(k) {
      shift_vector(
        base_order,
        k
      )
    }
  )

  if (
    mirror
  ) {

    reversed_order <- rev(
      base_order
    )

    perms <- c(
      perms,
      lapply(
        0:(
          n -
            1
        ),
        function(k) {
          shift_vector(
            reversed_order,
            k
          )
        }
      )
    )
  }

  perm_matrix <- unique(
    do.call(
      rbind,
      perms
    )
  )

  keep <- apply(
    perm_matrix,
    1,
    function(x) {
      !all(
        x ==
          base_order
      )
    }
  )

  perm_matrix[
    keep,
    ,
    drop = FALSE
  ]
}


core_series_list <- list()
counter <- 1

for (
  basin_i in
    basin_levels
) {

  for (
    metric_i in
      community_metrics
  ) {

    key_i <- paste(
      basin_i,
      metric_i,
      sep = "__"
    )

    d_comm <- core_comm_dist[
      [
        key_i
      ]
    ]

    n_time <- attr(
      d_comm,
      "Size"
    )

    perm_series <- make_cyclic_permutations(
      n = n_time,
      mirror = TRUE
    )

    for (
      env_i in
        vars_core
    ) {

      d_env <- build_temporal_env_dist(
        meta_df = meta_mantel,
        basin_name = basin_i,
        variable = env_i
      )

      if (
        !identical(
          attr(
            d_comm,
            "Labels"
          ),
          attr(
            d_env,
            "Labels"
          )
        )
      ) {
        stop(
          "Temporal label mismatch in cyclic test."
        )
      }

      mantel_series_i <- vegan::mantel(
        d_comm,
        d_env,
        method = "spearman",
        permutations =
          perm_series
      )

      core_series_list[
        [
          counter
        ]
      ] <- data.frame(
        Basin = basin_i,
        Community_metric =
          metric_i,
        Environment =
          env_i,
        n_timepoints =
          n_time,
        n_restricted_permutations =
          nrow(
            perm_series
          ),
        Mantel_r_series = unname(
          mantel_series_i$statistic
        ),
        p_series =
          mantel_series_i$signif,
        stringsAsFactors = FALSE
      )

      counter <- counter +
        1
    }
  }
}

core_series_sensitivity <- bind_rows(
  core_series_list
)

core_series_sensitivity$p_adj_BH_series <- p.adjust(
  core_series_sensitivity$p_series,
  method = "BH"
)

core_series_sensitivity <- core_series_sensitivity |>
  left_join(
    mantel_results |>
      select(
        Basin,
        Community_metric,
        Environment,
        Mantel_r,
        p_value,
        p_adj_BH_global
      ),
    by = c(
      "Basin",
      "Community_metric",
      "Environment"
    )
  ) |>
  rename(
    Mantel_r_original =
      Mantel_r,
    p_original =
      p_value,
    p_adj_BH_original =
      p_adj_BH_global
  )

write.table(
  core_series_sensitivity,
  file.path(
    tables_dir,
    "Mantel_core_temporal_sensitivity.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================================ #
# Geochemistry robustness | Leave-one-time-point-out
# ============================================================================ #

geo_loto_list <- list()
counter <- 1

for (
  basin_i in
    basin_levels
) {

  meta_geo_b <- meta_geo_29 |>
    filter(
      as.character(
        Basin
      ) ==
        basin_i
    ) |>
    arrange(
      Time_Point
    )

  for (
    metric_i in
      community_metrics
  ) {

    d_comm_full <- build_geo_comm_dist(
      ps_obj = ps_geo_29,
      meta_df = meta_geo_29,
      basin_name = basin_i,
      metric = metric_i
    )

    comm_matrix_full <- as.matrix(
      d_comm_full
    )

    comm_matrix_full <- comm_matrix_full[
      meta_geo_b$SampleID,
      meta_geo_b$SampleID,
      drop = FALSE
    ]

    for (
      drop_i in
        seq_len(
          nrow(
            meta_geo_b
          )
        )
    ) {

      drop_sample <- meta_geo_b$SampleID[
        drop_i
      ]

      drop_time <- meta_geo_b$Time_Point[
        drop_i
      ]

      keep_samples <- setdiff(
        meta_geo_b$SampleID,
        drop_sample
      )

      d_comm_loto <- as.dist(
        comm_matrix_full[
          keep_samples,
          keep_samples,
          drop = FALSE
        ]
      )

      meta_loto <- meta_geo_29 |>
        filter(
          SampleID !=
            drop_sample
        )

      for (
        geo_i in
          vars_geo_29
      ) {

        d_geo_loto <- build_geo_env_dist(
          meta_df = meta_loto,
          basin_name = basin_i,
          variable = geo_i
        )

        if (
          !identical(
            attr(
              d_comm_loto,
              "Labels"
            ),
            attr(
              d_geo_loto,
              "Labels"
            )
          )
        ) {
          stop(
            "Geochemical LOTO label mismatch."
          )
        }

        r_loto <- suppressWarnings(
          cor(
            as.vector(
              d_comm_loto
            ),
            as.vector(
              d_geo_loto
            ),
            method = "spearman"
          )
        )

        geo_loto_list[
          [
            counter
          ]
        ] <- data.frame(
          Basin = basin_i,
          Community_metric =
            metric_i,
          Geochemistry =
            geo_i,
          Dropped_Sample =
            drop_sample,
          Dropped_Time_Point =
            drop_time,
          Mantel_r_LOTO =
            r_loto,
          stringsAsFactors = FALSE
        )

        counter <- counter +
          1
      }
    }
  }
}

geo_loto_raw <- bind_rows(
  geo_loto_list
)

geo_loto_summary <- geo_loto_raw |>
  group_by(
    Basin,
    Community_metric,
    Geochemistry
  ) |>
  summarise(
    n_LOTO = n(),
    r_LOTO_min = min(
      Mantel_r_LOTO,
      na.rm = TRUE
    ),
    r_LOTO_median = median(
      Mantel_r_LOTO,
      na.rm = TRUE
    ),
    r_LOTO_max = max(
      Mantel_r_LOTO,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  left_join(
    mantel_geo_results |>
      select(
        Basin,
        Community_metric,
        Geochemistry,
        Mantel_r,
        p_value,
        p_adj_BH_global
      ),
    by = c(
      "Basin",
      "Community_metric",
      "Geochemistry"
    )
  )

geo_sign_stability <- geo_loto_raw |>
  left_join(
    mantel_geo_results |>
      select(
        Basin,
        Community_metric,
        Geochemistry,
        Mantel_r
      ),
    by = c(
      "Basin",
      "Community_metric",
      "Geochemistry"
    )
  ) |>
  group_by(
    Basin,
    Community_metric,
    Geochemistry
  ) |>
  summarise(
    same_sign_proportion = mean(
      sign(
        Mantel_r_LOTO
      ) ==
        sign(
          Mantel_r
        ),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

geo_loto_summary <- geo_loto_summary |>
  left_join(
    geo_sign_stability,
    by = c(
      "Basin",
      "Community_metric",
      "Geochemistry"
    )
  )

write.table(
  geo_loto_summary,
  file.path(
    tables_dir,
    "Mantel_geochemistry_LOTO_sensitivity.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================================ #
# Core robustness | Leave-one-time-point-out
# ============================================================================ #

core_loto_list <- list()
counter <- 1

for (
  basin_i in
    basin_levels
) {

  for (
    metric_i in
      community_metrics
  ) {

    key_i <- paste(
      basin_i,
      metric_i,
      sep = "__"
    )

    d_comm_full <- core_comm_dist[
      [
        key_i
      ]
    ]

    comm_matrix_full <- as.matrix(
      d_comm_full
    )

    for (
      env_i in
        vars_core
    ) {

      d_env_full <- build_temporal_env_dist(
        meta_df = meta_mantel,
        basin_name = basin_i,
        variable = env_i
      )

      env_matrix_full <- as.matrix(
        d_env_full
      )

      if (
        !identical(
          rownames(
            comm_matrix_full
          ),
          rownames(
            env_matrix_full
          )
        )
      ) {
        stop(
          "Core LOTO label mismatch."
        )
      }

      time_labels <- rownames(
        comm_matrix_full
      )

      for (
        drop_i in
          seq_along(
            time_labels
          )
      ) {

        drop_time <- time_labels[
          drop_i
        ]

        keep_times <- time_labels[
          -drop_i
        ]

        d_comm_loto <- as.dist(
          comm_matrix_full[
            keep_times,
            keep_times,
            drop = FALSE
          ]
        )

        d_env_loto <- as.dist(
          env_matrix_full[
            keep_times,
            keep_times,
            drop = FALSE
          ]
        )

        r_loto <- suppressWarnings(
          cor(
            as.vector(
              d_comm_loto
            ),
            as.vector(
              d_env_loto
            ),
            method = "spearman"
          )
        )

        core_loto_list[
          [
            counter
          ]
        ] <- data.frame(
          Basin = basin_i,
          Community_metric =
            metric_i,
          Environment =
            env_i,
          Dropped_Time_Point =
            drop_time,
          Mantel_r_LOTO =
            r_loto,
          stringsAsFactors = FALSE
        )

        counter <- counter +
          1
      }
    }
  }
}

core_loto_raw <- bind_rows(
  core_loto_list
)

core_loto_summary <- core_loto_raw |>
  group_by(
    Basin,
    Community_metric,
    Environment
  ) |>
  summarise(
    n_LOTO = n(),
    r_LOTO_min = min(
      Mantel_r_LOTO,
      na.rm = TRUE
    ),
    r_LOTO_median = median(
      Mantel_r_LOTO,
      na.rm = TRUE
    ),
    r_LOTO_max = max(
      Mantel_r_LOTO,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  left_join(
    mantel_results |>
      select(
        Basin,
        Community_metric,
        Environment,
        Mantel_r,
        p_value,
        p_adj_BH_global
      ),
    by = c(
      "Basin",
      "Community_metric",
      "Environment"
    )
  )

core_sign_stability <- core_loto_raw |>
  left_join(
    mantel_results |>
      select(
        Basin,
        Community_metric,
        Environment,
        Mantel_r
      ),
    by = c(
      "Basin",
      "Community_metric",
      "Environment"
    )
  ) |>
  group_by(
    Basin,
    Community_metric,
    Environment
  ) |>
  summarise(
    same_sign_proportion = mean(
      sign(
        Mantel_r_LOTO
      ) ==
        sign(
          Mantel_r
        ),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

core_loto_summary <- core_loto_summary |>
  left_join(
    core_sign_stability,
    by = c(
      "Basin",
      "Community_metric",
      "Environment"
    )
  )

write.table(
  core_loto_summary,
  file.path(
    tables_dir,
    "Mantel_core_LOTO_sensitivity.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================================ #
# Supplementary Table 5 | Manuscript-style summary
# ============================================================================ #

table_s5_core <- mantel_results |>
  transmute(
    Basin,
    Community_metric,
    Variable_group =
      "Core physicochemical",
    Variable =
      Environment,
    Mantel_r,
    P = p_value,
    BH_FDR_q =
      p_adj_BH_global,
    BH_FDR_lt_0_05 =
      ifelse(
        p_adj_BH_global <
          0.05,
        "Yes",
        "No"
      )
  )

# The current Supplementary Table 5 lists all 24 core tests and the
# FDR-supported geochemical associations.
table_s5_geo <- mantel_geo_results |>
  filter(
    p_adj_BH_global <
      0.05
  ) |>
  transmute(
    Basin,
    Community_metric,
    Variable_group =
      "Geochemistry",
    Variable =
      Geochemistry,
    Mantel_r,
    P = p_value,
    BH_FDR_q =
      p_adj_BH_global,
    BH_FDR_lt_0_05 =
      "Yes"
  )

table_s5 <- bind_rows(
  table_s5_core,
  table_s5_geo
)

openxlsx::write.xlsx(
  list(
    Table_Supplementary_5 =
      table_s5,
    Core_all_24 =
      mantel_results,
    Geochemistry_all_54 =
      mantel_geo_results,
    Core_LOTO =
      core_loto_summary,
    Geochemistry_LOTO =
      geo_loto_summary,
    Core_cyclic_permutations =
      core_series_sensitivity
  ),
  file.path(
    tables_dir,
    "Table_Supplementary5_Mantel.xlsx"
  ),
  overwrite = TRUE
)


# ============================================================================ #
# Validation summary
# ============================================================================ #

validation_summary <- tibble(
  Analysis_block = c(
    "Core physicochemical",
    "Geochemistry"
  ),
  Total_tests = c(
    nrow(
      mantel_results
    ),
    nrow(
      mantel_geo_results
    )
  ),
  BH_significant = c(
    sum(
      mantel_results$p_adj_BH_global <
        0.05,
      na.rm = TRUE
    ),
    sum(
      mantel_geo_results$p_adj_BH_global <
        0.05,
      na.rm = TRUE
    )
  ),
  Minimum_LOTO_same_sign = c(
    min(
      core_loto_summary$same_sign_proportion[
        core_loto_summary$p_adj_BH_global <
          0.05
      ],
      na.rm = TRUE
    ),
    min(
      geo_loto_summary$same_sign_proportion[
        geo_loto_summary$p_adj_BH_global <
          0.05
      ],
      na.rm = TRUE
    )
  )
)

openxlsx::write.xlsx(
  validation_summary,
  file.path(
    tables_dir,
    "Mantel_validation_summary.xlsx"
  ),
  overwrite = TRUE
)


# ============================================================================ #
# Final summary
# ============================================================================ #

cat(
  "\nMantel community-environment analysis completed.\n"
)

print(
  validation_summary
)

cat(
  "\nExpected historical totals:\n",
  "  Core: 24 tests, 18 BH-FDR significant.\n",
  "  Geochemistry: 54 tests, 9 BH-FDR significant.\n",
  sep = ""
)

cat(
  "\nThe core and geochemical Mantel figure components were exported ",
  "separately, as in the master script.\n",
  sep = ""
)
