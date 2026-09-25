# ============================================================================ #
# 15_environmental_beta_diversity.R
# ============================================================================ #
#
# Purpose:
# Reproduce the pairwise community-environment relationships supporting
# Supplementary Figures 5 and 6.
#
# Supplementary Figure 5:
#   Pairwise environmental differences vs weighted UniFrac distance.
#
# Supplementary Figure 6:
#   Pairwise environmental differences vs Bray-Curtis dissimilarity.
#
# Environmental variables:
#   - Water content (Wc)
#   - Salinity
#   - pH
#   - Oxidation-reduction potential (ORP)
#
# Historical analysis:
# - non-rarefied final phyloseq object
# - Bray-Curtis calculated on sample relative abundances
# - weighted UniFrac calculated on the final phyloseq object
# - only within-basin sample pairs retained
# - absolute pairwise environmental differences
# - basin-specific Gaussian GLM:
#
#       Community_distance ~ Environmental_difference
#
# - basin-specific Mantel tests:
#       Spearman correlation
#       9,999 permutations
#
# Curation note:
# The master script contains repeated variable-specific blocks with a
# changeable distance selector. Literal saved selector states are not present
# for every variable × distance combination used in the final two
# supplementary figures. This curated version therefore evaluates BOTH final
# metrics explicitly with the same historical calculations instead of
# requiring manual selector changes and reruns.
#
# Exploratory analyses of within-environmental-bin community heterogeneity
# and Jaccard distance are not included because they are not part of the
# final Supplementary Figures 5-6.
#
# Input:
# - output_phyloseq_tree/ps_mdecon_v2_FTree.rds
#
# Main outputs:
# - Supplementary_Figure5_weighted_UniFrac_environment.*
# - Supplementary_Figure6_Bray_Curtis_environment.*
# - pairwise source-data CSV files
# - basin-specific GLM and Mantel result tables
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
library(broom)
library(ggplot2)
library(cowplot)
library(patchwork)
library(scales)


# ============================================================================ #
# Paths and reproducibility
# ============================================================================ #

set.seed(2026)

phyloseq_path <- file.path(
  "output_phyloseq_tree",
  "ps_mdecon_v2_FTree.rds"
)

output_dir <- "output_environment_beta_diversity"

figure_dir <- file.path(
  output_dir,
  "figures"
)

table_dir <- file.path(
  output_dir,
  "tables"
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
# Load final phyloseq object
# ============================================================================ #

ps <- readRDS(
  phyloseq_path
)

ps_ra <- transform_sample_counts(
  ps,
  function(x) {
    100 *
      x /
      sum(
        x
      )
  }
)

metadata <- data.frame(
  sample_data(
    ps
  ),
  check.names = FALSE
)

metadata$SampleID <- rownames(
  metadata
)


# ============================================================================ #
# Factor order
# ============================================================================ #

basin_levels <- c(
  "Herradura Clay Pan",
  "Ckoirama Halite Field",
  "Yungay Station Basin"
)

metadata$Ephemeral_Basin <- factor(
  as.character(
    metadata$Ephemeral_Basin
  ),
  levels = basin_levels
)


# ============================================================================ #
# Environmental-variable settings
# ============================================================================ #

environment_settings <- list(
  Wc = list(
    column = "Wc",
    label = "Water content difference (%)",
    short = "water_content",
    line_color = "#6baed6",
    jitter_width = 0.5
  ),
  Salinity = list(
    column = "Salinity",
    label = "Salinity difference (PSU)",
    short = "salinity",
    line_color = "red",
    jitter_width = 1.5
  ),
  pH = list(
    column = "pH",
    label = "pH difference",
    short = "pH",
    line_color = "#31a354",
    jitter_width = 0.02
  ),
  ORP = list(
    column = "ORP",
    label = "ORP difference (mV)",
    short = "ORP",
    line_color = "#756bb1",
    jitter_width = 1.5
  )
)


# ============================================================================ #
# Distance helpers
# ============================================================================ #

get_distance_label <- function(
    distance_method
) {

  switch(
    distance_method,
    bray =
      "Bray-Curtis community dissimilarity",
    wunifrac =
      "Weighted UniFrac distance",
    stop(
      "Unsupported distance method."
    )
  )
}


get_distance_tag <- function(
    distance_method
) {

  switch(
    distance_method,
    bray =
      "bray",
    wunifrac =
      "wunifrac",
    stop(
      "Unsupported distance method."
    )
  )
}


compute_community_distance <- function(
    ps_object,
    ps_ra_object,
    distance_method
) {

  switch(
    distance_method,
    bray =
      phyloseq::distance(
        ps_ra_object,
        method = "bray"
      ),
    wunifrac =
      phyloseq::distance(
        ps_object,
        method = "wunifrac"
      ),
    stop(
      "Unsupported distance method."
    )
  )
}


# ============================================================================ #
# Format p values
# ============================================================================ #

format_p_plot <- function(
    p
) {

  dplyr::case_when(
    is.na(
      p
    ) ~
      NA_character_,
    p <
      0.001 ~
      "< 0.001",
    p <
      0.01 ~
      "< 0.01",
    p <
      0.05 ~
      "< 0.05",
    TRUE ~
      paste0(
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
# Build unique within-basin sample pairs
# ============================================================================ #

build_pairwise_environment_data <- function(
    distance_object,
    metadata_table,
    environmental_variable
) {

  distance_matrix <- as.matrix(
    distance_object
  )

  distance_long <- as.data.frame(
    distance_matrix
  ) |>
    rownames_to_column(
      "Sample1"
    ) |>
    pivot_longer(
      cols = -Sample1,
      names_to = "Sample2",
      values_to = "Community_distance"
    ) |>
    filter(
      Sample1 <
        Sample2
    )

  meta_1 <- metadata_table |>
    transmute(
      Sample1 = SampleID,
      Basin1 = Ephemeral_Basin,
      Env1 = as.numeric(
        .data[
          [
            environmental_variable
          ]
        ]
      )
    )

  meta_2 <- metadata_table |>
    transmute(
      Sample2 = SampleID,
      Basin2 = Ephemeral_Basin,
      Env2 = as.numeric(
        .data[
          [
            environmental_variable
          ]
        ]
      )
    )

  distance_long |>
    left_join(
      meta_1,
      by = "Sample1"
    ) |>
    left_join(
      meta_2,
      by = "Sample2"
    ) |>
    filter(
      Basin1 ==
        Basin2,
      !is.na(
        Env1
      ),
      !is.na(
        Env2
      ),
      !is.na(
        Community_distance
      )
    ) |>
    mutate(
      Ephemeral_Basin = factor(
        Basin1,
        levels = basin_levels
      ),
      Environmental_difference = abs(
        Env2 -
          Env1
      )
    ) |>
    select(
      Sample1,
      Sample2,
      Ephemeral_Basin,
      Env1,
      Env2,
      Environmental_difference,
      Community_distance
    )
}


# ============================================================================ #
# Basin-specific Gaussian GLM
# ============================================================================ #

fit_pairwise_glm <- function(
    pairwise_data
) {

  pairwise_data |>
    group_by(
      Ephemeral_Basin
    ) |>
    nest() |>
    mutate(
      model = map(
        data,
        ~ glm(
          Community_distance ~
            Environmental_difference,
          data = .x,
          family = gaussian(
            link = "identity"
          )
        )
      ),
      tidy = map(
        model,
        broom::tidy
      )
    ) |>
    unnest(
      tidy
    ) |>
    filter(
      term ==
        "Environmental_difference"
    ) |>
    transmute(
      Ephemeral_Basin,
      Slope = estimate,
      SE = std.error,
      Statistic = statistic,
      P_value = p.value
    )
}


# ============================================================================ #
# Basin-specific Mantel tests
# ============================================================================ #

run_mantel_by_basin <- function(
    ps_object,
    ps_ra_object,
    metadata_table,
    environmental_variable,
    distance_method,
    permutations = 9999
) {

  metadata_complete <- metadata_table |>
    select(
      SampleID,
      Ephemeral_Basin,
      all_of(
        environmental_variable
      )
    ) |>
    filter(
      !is.na(
        .data[
          [
            environmental_variable
          ]
        ]
      ),
      !is.na(
        Ephemeral_Basin
      )
    )

  results <- list()

  for (
    basin_i in basin_levels
  ) {

    meta_basin <- metadata_complete |>
      filter(
        Ephemeral_Basin ==
          basin_i
      )

    sample_ids <- meta_basin$SampleID

    if (
      length(
        sample_ids
      ) <
        3
    ) {
      next
    }

    ps_basin <- prune_samples(
      sample_ids,
      ps_object
    )

    ps_ra_basin <- prune_samples(
      sample_ids,
      ps_ra_object
    )

    community_distance <- compute_community_distance(
      ps_object = ps_basin,
      ps_ra_object = ps_ra_basin,
      distance_method = distance_method
    )

    community_labels <- attr(
      community_distance,
      "Labels"
    )

    environmental_values <- meta_basin[
      [
        environmental_variable
      ]
    ][
      match(
        community_labels,
        meta_basin$SampleID
      )
    ]

    environmental_distance <- dist(
      as.numeric(
        environmental_values
      ),
      method = "manhattan"
    )

    attr(
      environmental_distance,
      "Labels"
    ) <- community_labels

    set.seed(
      2026
    )

    mantel_result <- vegan::mantel(
      xdis = community_distance,
      ydis = environmental_distance,
      method = "spearman",
      permutations = permutations,
      na.rm = TRUE
    )

    results[
      [
        basin_i
      ]
    ] <- tibble(
      Ephemeral_Basin = basin_i,
      Environmental_variable =
        environmental_variable,
      Distance =
        distance_method,
      Mantel_method =
        "spearman",
      Permutations =
        permutations,
      Mantel_r =
        as.numeric(
          mantel_result$statistic
        ),
      P_value =
        mantel_result$signif
    )
  }

  bind_rows(
    results
  ) |>
    mutate(
      Ephemeral_Basin = factor(
        Ephemeral_Basin,
        levels = basin_levels
      )
    )
}


# ============================================================================ #
# Figure helper
# ============================================================================ #

make_environment_distance_plot <- function(
    pairwise_data,
    glm_results,
    mantel_results,
    settings,
    distance_method,
    show_x_title = TRUE
) {

  distance_label <- get_distance_label(
    distance_method
  )

  strip_labels <- setNames(
    paste0(
      as.character(
        mantel_results$Ephemeral_Basin
      ),
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
    as.character(
      mantel_results$Ephemeral_Basin
    )
  )

  label_data <- glm_results |>
    mutate(
      Label = paste0(
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
      x = Environmental_difference,
      y = Community_distance
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
        width =
          settings$jitter_width,
        height = 0
      )
    ) +
    geom_smooth(
      method = "glm",
      method.args = list(
        family = gaussian(
          link = "identity"
        )
      ),
      se = TRUE,
      color =
        settings$line_color,
      fill = "grey50",
      alpha = 0.4,
      linewidth = 1.15
    ) +
    facet_wrap(
      ~ Ephemeral_Basin,
      nrow = 1,
      labeller = labeller(
        Ephemeral_Basin =
          strip_labels
      ),
      scales = "free_x"
    ) +
    geom_text(
      data = label_data,
      aes(
        x = Inf,
        y = -Inf,
        label = Label
      ),
      hjust = 1.02,
      vjust = -0.2,
      size = 3.0,
      colour = "black",
      lineheight = 1.05,
      inherit.aes = FALSE
    ) +
    labs(
      x = if (
        show_x_title
      ) {
        settings$label
      } else {
        NULL
      },
      y = distance_label
    ) +
    coord_cartesian(
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
        size = 10,
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
        size = 10.5
      ),
      axis.text = element_text(
        size = 8.5
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
# Run one environmental variable × distance metric
# ============================================================================ #

run_environment_distance_analysis <- function(
    environmental_name,
    distance_method
) {

  settings <- environment_settings[
    [
      environmental_name
    ]
  ]

  environmental_variable <-
    settings$column

  community_distance <- compute_community_distance(
    ps_object = ps,
    ps_ra_object = ps_ra,
    distance_method = distance_method
  )

  pairwise_data <- build_pairwise_environment_data(
    distance_object = community_distance,
    metadata_table = metadata,
    environmental_variable =
      environmental_variable
  )

  glm_results <- fit_pairwise_glm(
    pairwise_data
  )

  mantel_results <- run_mantel_by_basin(
    ps_object = ps,
    ps_ra_object = ps_ra,
    metadata_table = metadata,
    environmental_variable =
      environmental_variable,
    distance_method =
      distance_method,
    permutations = 9999
  )

  figure <- make_environment_distance_plot(
    pairwise_data =
      pairwise_data,
    glm_results =
      glm_results,
    mantel_results =
      mantel_results,
    settings =
      settings,
    distance_method =
      distance_method
  )

  distance_tag <- get_distance_tag(
    distance_method
  )

  write.csv(
    pairwise_data,
    file.path(
      table_dir,
      paste0(
        "pairwise_",
        settings$short,
        "_",
        distance_tag,
        ".csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    glm_results,
    file.path(
      table_dir,
      paste0(
        "GLM_",
        settings$short,
        "_",
        distance_tag,
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
        "Mantel_",
        settings$short,
        "_",
        distance_tag,
        "_spearman_9999.csv"
      )
    ),
    row.names = FALSE
  )

  list(
    plot = figure,
    pairwise = pairwise_data,
    glm = glm_results,
    mantel = mantel_results
  )
}


# ============================================================================ #
# Run Supplementary Figure 5 | weighted UniFrac
# ============================================================================ #

wunifrac_results <- list()

for (
  environmental_name in
    names(
      environment_settings
    )
) {

  wunifrac_results[
    [
      environmental_name
    ]
  ] <- run_environment_distance_analysis(
    environmental_name =
      environmental_name,
    distance_method =
      "wunifrac"
  )
}

supplementary_figure5 <- (
  wunifrac_results$Wc$plot /
    wunifrac_results$Salinity$plot /
    wunifrac_results$pH$plot /
    wunifrac_results$ORP$plot
) +
  plot_annotation(
    tag_levels = "A"
  )

ggsave(
  file.path(
    figure_dir,
    "Supplementary_Figure5_weighted_UniFrac_environment.tiff"
  ),
  plot =
    supplementary_figure5,
  width = 12,
  height = 15,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Supplementary_Figure5_weighted_UniFrac_environment.svg"
  ),
  plot =
    supplementary_figure5,
  width = 12,
  height = 15,
  units = "in"
)


# ============================================================================ #
# Run Supplementary Figure 6 | Bray-Curtis
# ============================================================================ #

bray_results <- list()

for (
  environmental_name in
    names(
      environment_settings
    )
) {

  bray_results[
    [
      environmental_name
    ]
  ] <- run_environment_distance_analysis(
    environmental_name =
      environmental_name,
    distance_method =
      "bray"
  )
}

supplementary_figure6 <- (
  bray_results$Wc$plot /
    bray_results$Salinity$plot /
    bray_results$pH$plot /
    bray_results$ORP$plot
) +
  plot_annotation(
    tag_levels = "A"
  )

ggsave(
  file.path(
    figure_dir,
    "Supplementary_Figure6_Bray_Curtis_environment.tiff"
  ),
  plot =
    supplementary_figure6,
  width = 12,
  height = 15,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Supplementary_Figure6_Bray_Curtis_environment.svg"
  ),
  plot =
    supplementary_figure6,
  width = 12,
  height = 15,
  units = "in"
)


# ============================================================================ #
# Combined Mantel summary
# ============================================================================ #

mantel_summary <- bind_rows(
  lapply(
    names(
      wunifrac_results
    ),
    function(v) {
      wunifrac_results[
        [
          v
        ]
      ]$mantel |>
        mutate(
          Variable = v,
          .before = 1
        )
    }
  ),
  lapply(
    names(
      bray_results
    ),
    function(v) {
      bray_results[
        [
          v
        ]
      ]$mantel |>
        mutate(
          Variable = v,
          .before = 1
        )
    }
  )
) |>
  mutate(
    P_adjust_BH_within_distance =
      ave(
        P_value,
        Distance,
        FUN = function(x) {
          p.adjust(
            x,
            method = "BH"
          )
        }
      )
  )

write.csv(
  mantel_summary,
  file.path(
    table_dir,
    "Supplementary_Figures5-6_Mantel_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================================ #
# Combined GLM summary
# ============================================================================ #

glm_summary <- bind_rows(
  lapply(
    names(
      wunifrac_results
    ),
    function(v) {
      wunifrac_results[
        [
          v
        ]
      ]$glm |>
        mutate(
          Variable = v,
          Distance = "wunifrac",
          .before = 1
        )
    }
  ),
  lapply(
    names(
      bray_results
    ),
    function(v) {
      bray_results[
        [
          v
        ]
      ]$glm |>
        mutate(
          Variable = v,
          Distance = "bray",
          .before = 1
        )
    }
  )
)

write.csv(
  glm_summary,
  file.path(
    table_dir,
    "Supplementary_Figures5-6_GLM_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================================ #
# Final summary
# ============================================================================ #

cat(
  "\nEnvironmental beta-diversity analysis completed.\n"
)

cat(
  "\nSupplementary Figure 5: weighted UniFrac.\n"
)

cat(
  "Supplementary Figure 6: Bray-Curtis.\n"
)

cat(
  "\nEach figure contains Wc, salinity, pH, and ORP pairwise differences.\n"
)

cat(
  "Mantel tests: Spearman, 9,999 permutations, basin-specific.\n"
)

cat(
  "Pairwise regression: Gaussian GLM within each basin.\n"
)
