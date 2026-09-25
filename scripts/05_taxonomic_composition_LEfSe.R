# ============================================================================ #
# 05_taxonomic_composition_LEfSe.R
# ============================================================================ #
#
# Purpose:
# Reproduce the taxonomic-composition and LEfSe analyses used in Figure 4:
#
#   - Figure 4A: order-level relative abundance through time
#   - Figure 4B: family-level LEfSe biomarkers distinguishing basins
#                within the Wet and Dry phases
#   - Figure 4C: family-level LEfSe biomarkers distinguishing Wet and Dry
#                phases within each basin
#
# Notes:
# - Taxonomic composition is calculated from non-rarefied counts.
# - Relative abundance is calculated independently within each sample.
# - Selected orders are displayed individually; remaining orders are grouped
#   as "Other", while missing order-level assignments are "Unclassified".
# - Samples absent from the final sequencing dataset but needed to preserve
#   the complete time-series layout are added as zero-abundance placeholders.
#   They do not contribute sequence abundance.
# - LEfSe is performed at Family level with CPM normalization, P < 0.05 for
#   Kruskal-Wallis and Wilcoxon tests, and an LDA cutoff of 4.
# - ANCOM-BC2, cladograms, and other exploratory alternatives present in the
#   master script are intentionally omitted.
#
# Input:
# - output_phyloseq_tree/ps_mdecon_v2_FTree.rds
#
# Main outputs:
# - Figure4A_order_relative_abundance.*
# - Figure4B_LEfSe_basins_by_phase.*
# - Figure4C_LEfSe_wet_dry_by_basin.*
# - LEfSe marker tables
#
# ============================================================================ #


# ============================================================================ #
# Packages
# ============================================================================ #

library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(openxlsx)
library(microbiomeMarker)
library(svglite)


# ============================================================================ #
# Reproducibility
# ============================================================================ #

set.seed(2026)

Sys.setlocale(
  "LC_TIME",
  "C"
)


# ============================================================================ #
# User paths
# ============================================================================ #

phyloseq_path <- file.path(
  "output_phyloseq_tree",
  "ps_mdecon_v2_FTree.rds"
)

output_dir <- "output_taxonomic_composition"
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

basin_colors <- c(
  "Herradura Clay Pan" = "#1b9e77",
  "Ckoirama Halite Field" = "#d95f02",
  "Yungay Station Basin" = "#7570b3"
)

phase_colors <- c(
  "Wet" = "#80A6CC",
  "Dry" = "#D9B35F"
)


# ============================================================================ #
# Load phyloseq object
# ============================================================================ #

ps <- readRDS(
  phyloseq_path
)

sample_data(ps)$Ephemeral_Basin <- factor(
  sample_data(ps)$Ephemeral_Basin,
  levels = basin_levels
)

sample_data(ps)$Phase <- factor(
  sample_data(ps)$Phase,
  levels = c(
    "Wet",
    "Dry"
  )
)


# ============================================================================ #
# Figure 4A | Order-level relative abundance
# ============================================================================ #

# -------------------------------------------------------------------------- #
# Transform counts to sample-wise relative abundance (%)
# -------------------------------------------------------------------------- #

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

otu_table(ps_ra)[
  is.na(
    otu_table(ps_ra)
  )
] <- 0


# -------------------------------------------------------------------------- #
# Agglomerate at Order level
# -------------------------------------------------------------------------- #

ps_order <- tax_glom(
  ps_ra,
  taxrank = "Order",
  NArm = FALSE,
  bad_empty = c(
    "",
    " ",
    "\t"
  )
)

cols_keep <- c(
  "OTU",
  "Sample",
  "Abundance",
  "Sample_Time",
  "Time",
  "Ephemeral_Basin",
  "Point",
  "Sampling_Point",
  "Basin_Time",
  "Matrix",
  "Date",
  "Kingdom",
  "Phylum",
  "Class",
  "Order"
)

ra_order <- psmelt(
  ps_order
) |>
  dplyr::mutate(
    dplyr::across(
      c(
        Kingdom,
        Phylum,
        Class,
        Order
      ),
      function(x) {

        x <- as.character(
          x
        )

        ifelse(
          is.na(x) |
            x == "",
          "Unclassified",
          x
        )
      }
    ),
    Date = as.Date(
      as.character(
        Date
      )
    )
  ) |>
  dplyr::select(
    dplyr::all_of(
      cols_keep
    )
  )


# -------------------------------------------------------------------------- #
# Add zero-abundance placeholders for missing time-series samples
# -------------------------------------------------------------------------- #

# These samples were absent from the final sequence table but were inserted
# as zero-abundance placeholders in the original plotting workflow so that
# missing dates remained visible in the longitudinal layout.

missing_samples <- c(
  "CHF3_T16",
  "CHF1_T18",
  "YSB1_T1",
  "YSB2_T1",
  "YSB3_T1",
  "YSB3_T20",
  "YSB2_T18",
  "YSB1_T9",
  "YSB1_T3"
)

missing_data <- data.frame(
  OTU = "ASV_placeholder",
  Sample = missing_samples,
  Abundance = 0,
  Sample_Time = missing_samples,
  Point = c(
    "Point 3",
    "Point 1",
    "Point 1",
    "Point 2",
    "Point 3",
    "Point 3",
    "Point 2",
    "Point 1",
    "Point 1"
  ),
  Time = c(
    "T16",
    "T18",
    "T1",
    "T1",
    "T1",
    "T20",
    "T18",
    "T9",
    "T3"
  ),
  Ephemeral_Basin = c(
    rep(
      "Ckoirama Halite Field",
      2
    ),
    rep(
      "Yungay Station Basin",
      7
    )
  ),
  Sampling_Point = c(
    "CHF3",
    "CHF1",
    "YSB1",
    "YSB2",
    "YSB3",
    "YSB3",
    "YSB2",
    "YSB1",
    "YSB1"
  ),
  Basin_Time = c(
    "CHF_T16",
    "CHF_T18",
    "YSB_T1",
    "YSB_T1",
    "YSB_T1",
    "YSB_T20",
    "YSB_T18",
    "YSB_T9",
    "YSB_T3"
  ),
  Matrix = c(
    "Sediment",
    "Sediment",
    "Water",
    "Water",
    "Water",
    "Sediment",
    "Sediment",
    "Water",
    "Water"
  ),
  Date = as.Date(
    c(
      "2018-03-28",
      "2018-07-03",
      "2017-06-16",
      "2017-06-16",
      "2017-06-16",
      "2018-09-14",
      "2018-07-03",
      "2017-07-14",
      "2017-06-24"
    )
  ),
  Kingdom = "Archaea",
  Phylum = "Halobacteriota",
  Class = "Halobacteria",
  Order = "Halobacterales",
  stringsAsFactors = FALSE
)

ra_order <- dplyr::bind_rows(
  ra_order,
  missing_data
)


# -------------------------------------------------------------------------- #
# Orders displayed individually in Figure 4A
# -------------------------------------------------------------------------- #

order_levels <- c(
  "Azospirillales",
  "Bacillales",
  "Bacteroidales",
  "Bradymonadales",
  "Burkholderiales",
  "Cardiobacteriales",
  "Caulobacterales",
  "Cytophagales",
  "Deinococcales",
  "Enterobacterales",
  "Flavobacteriales",
  "Halobacterales",
  "Hyphomicrobiales",
  "Lachnospirales",
  "Longimicrobiales",
  "Lysobacterales",
  "Micrococcales",
  "Mycobacteriales",
  "Nanosalinales",
  "Oscillospirales",
  "Paenibacillales",
  "Propionibacteriales",
  "Pseudomonadales",
  "Rhodobacterales",
  "Sphingobacteriales",
  "Sphingomonadales",
  "Other",
  "Unclassified"
)

selected_orders <- setdiff(
  order_levels,
  c(
    "Other",
    "Unclassified"
  )
)

ra_order <- ra_order |>
  dplyr::mutate(
    Order = dplyr::case_when(
      is.na(Order) |
        Order == "" ~ "Unclassified",
      Order %in% c(
        selected_orders,
        "Unclassified"
      ) ~ Order,
      TRUE ~ "Other"
    )
  )


# -------------------------------------------------------------------------- #
# Sum relative abundances after grouping non-displayed orders as Other
# -------------------------------------------------------------------------- #

group_cols <- c(
  "Sample_Time",
  "Sample",
  "Time",
  "Date",
  "Ephemeral_Basin",
  "Sampling_Point",
  "Point",
  "Basin_Time",
  "Matrix",
  "Order"
)

ra_order <- ra_order |>
  dplyr::group_by(
    dplyr::across(
      dplyr::all_of(
        group_cols
      )
    )
  ) |>
  dplyr::summarise(
    Kingdom = ifelse(
      dplyr::n_distinct(
        Kingdom
      ) == 1,
      dplyr::first(
        Kingdom
      ),
      "Mixed"
    ),
    Phylum = ifelse(
      dplyr::n_distinct(
        Phylum
      ) == 1,
      dplyr::first(
        Phylum
      ),
      "Mixed"
    ),
    Class = ifelse(
      dplyr::n_distinct(
        Class
      ) == 1,
      dplyr::first(
        Class
      ),
      "Mixed"
    ),
    Abundance = sum(
      Abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# -------------------------------------------------------------------------- #
# Factor ordering
# -------------------------------------------------------------------------- #

samples_by_time <- unlist(
  lapply(
    1:20,
    function(t) {
      paste0(
        rep(
          c(
            "HCP",
            "CHF",
            "YSB"
          ),
          each = 3
        ),
        rep(
          1:3,
          times = 3
        ),
        "_T",
        t
      )
    }
  )
)

sample_levels <- c(
  paste0(
    "YB",
    1:5,
    "_T0"
  ),
  samples_by_time
)

order_levels_rank <- c(
  "Other",
  "Unclassified",
  "Nanosalinales",
  "Longimicrobiales",
  "Bradymonadales",
  "Oscillospirales",
  "Lachnospirales",
  "Deinococcales",
  "Paenibacillales",
  "Bacillales",
  "Mycobacteriales",
  "Propionibacteriales",
  "Micrococcales",
  "Bacteroidales",
  "Cytophagales",
  "Sphingobacteriales",
  "Flavobacteriales",
  "Azospirillales",
  "Rhodobacterales",
  "Hyphomicrobiales",
  "Caulobacterales",
  "Sphingomonadales",
  "Cardiobacteriales",
  "Enterobacterales",
  "Lysobacterales",
  "Burkholderiales",
  "Pseudomonadales",
  "Halobacterales"
)

ra_order <- ra_order |>
  dplyr::mutate(
    Sample_Time = factor(
      Sample_Time,
      levels = sample_levels
    ),
    Ephemeral_Basin = factor(
      Ephemeral_Basin,
      levels = basin_levels
    ),
    Time = factor(
      Time,
      levels = paste0(
        "T",
        1:20
      )
    ),
    Order = factor(
      Order,
      levels = order_levels
    ),
    Date_plot = factor(
      Date,
      levels = sort(
        unique(
          Date
        )
      )
    )
  ) |>
  dplyr::arrange(
    Ephemeral_Basin,
    Point,
    Date,
    Order
  )

ra_order_rank <- ra_order |>
  dplyr::mutate(
    Order = factor(
      as.character(
        Order
      ),
      levels = order_levels_rank
    )
  )

if (
  any(
    is.na(
      ra_order_rank$Order
    )
  )
) {
  stop(
    "Some Order values became NA after factor ordering."
  )
}


# -------------------------------------------------------------------------- #
# Order colors used in the original figure
# -------------------------------------------------------------------------- #

order_colors <- c(
  "Other" = "#B8B8B8",
  "Unclassified" = "#495057",
  "Nanosalinales" = "gold",
  "Longimicrobiales" = "chartreuse4",
  "Bradymonadales" = "lightpink",
  "Oscillospirales" = "peachpuff4",
  "Lachnospirales" = "lightgoldenrod",
  "Deinococcales" = "#ae017e",
  "Paenibacillales" = "#F4A261",
  "Bacillales" = "#e6550d",
  "Mycobacteriales" = "darkred",
  "Propionibacteriales" = "#fc9272",
  "Micrococcales" = "firebrick2",
  "Bacteroidales" = "olivedrab1",
  "Cytophagales" = "limegreen",
  "Sphingobacteriales" = "#66C2A5",
  "Flavobacteriales" = "#d9f0a3",
  "Azospirillales" = "#c7e9b4",
  "Rhodobacterales" = "#2c7fb8",
  "Hyphomicrobiales" = "#41b6c4",
  "Caulobacterales" = "#a1dab4",
  "Sphingomonadales" = "#253494",
  "Cardiobacteriales" = "#AAB0E6",
  "Enterobacterales" = "#D6BCEB",
  "Lysobacterales" = "#810f7c",
  "Burkholderiales" = "#dadaeb",
  "Pseudomonadales" = "#7A1FA2",
  "Halobacterales" = "#48D1CC"
)


# -------------------------------------------------------------------------- #
# Plot Figure 4A
# -------------------------------------------------------------------------- #

date_labels <- function(x) {
  tools::toTitleCase(
    format(
      as.Date(
        as.character(
          x
        )
      ),
      "%d %b"
    )
  )
}

figure4a <- ggplot(
  ra_order_rank,
  aes(
    Date_plot,
    Abundance,
    fill = Order
  )
) +
  geom_col(
    width = 0.95
  ) +
  facet_grid(
    rows = vars(Point),
    cols = vars(Ephemeral_Basin),
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = order_colors,
    na.value = "grey80"
  ) +
  scale_x_discrete(
    labels = date_labels,
    expand = expansion(
      add = c(
        0,
        0
      )
    )
  ) +
  scale_y_continuous(
    limits = c(
      0,
      100
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  guides(
    fill = guide_legend(
      ncol = 1
    )
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "Relative abundance (%)"
  ) +
  theme_cowplot() +
  theme(
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.5
    ),
    panel.spacing.x = unit(
      1.2,
      "lines"
    ),
    panel.spacing.y = unit(
      0.75,
      "lines"
    ),
    strip.background = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "right"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure4A_order_relative_abundance.tiff"
  ),
  plot = figure4a,
  width = 18,
  height = 10,
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure4A_order_relative_abundance.svg"
  ),
  plot = figure4a,
  width = 18,
  height = 10,
  device = "svg"
)


# -------------------------------------------------------------------------- #
# Export Figure 4A plotting data
# -------------------------------------------------------------------------- #

ra_order_plot_export <- ra_order_rank |>
  dplyr::mutate(
    dplyr::across(
      where(
        is.factor
      ),
      as.character
    )
  ) |>
  dplyr::arrange(
    Ephemeral_Basin,
    Sampling_Point,
    Date,
    Order
  )

ra_order_sample_totals <- ra_order_plot_export |>
  dplyr::group_by(
    Sample,
    Sample_Time,
    Date,
    Ephemeral_Basin,
    Sampling_Point,
    Point,
    Matrix
  ) |>
  dplyr::summarise(
    Total_abundance = sum(
      Abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

ra_order_legend <- data.frame(
  Order = order_levels_rank,
  Color = unname(
    order_colors[
      order_levels_rank
    ]
  ),
  stringsAsFactors = FALSE
)

openxlsx::write.xlsx(
  list(
    plot_data = ra_order_plot_export,
    sample_totals = ra_order_sample_totals,
    order_legend = ra_order_legend
  ),
  file = file.path(
    table_dir,
    "Figure4A_order_relative_abundance_data.xlsx"
  )
)


# ============================================================================ #
# LEfSe helper functions
# ============================================================================ #

taxa_rank_lefse <- "Family"
lda_cutoff_lefse <- 4

run_lefse_safe <- function(
    ps_object,
    group_var
) {

  ps_object <- prune_taxa(
    taxa_sums(
      ps_object
    ) > 0,
    ps_object
  )

  microbiomeMarker::run_lefse(
    ps_object,
    group = group_var,
    taxa_rank = taxa_rank_lefse,
    transform = "identity",
    norm = "CPM",
    kw_cutoff = 0.05,
    wilcoxon_cutoff = 0.05,
    lda_cutoff = lda_cutoff_lefse
  )
}


clean_family_name <- function(x) {

  x <- gsub(
    "^.*_f__",
    "",
    x
  )

  x <- gsub(
    "^f__",
    "",
    x
  )

  ifelse(
    is.na(x) |
      x == "",
    "Unclassified family",
    x
  )
}


# ============================================================================ #
# Figure 4C | Wet vs Dry biomarkers within each basin
# ============================================================================ #

run_phase_lefse_basin <- function(
    ps_object,
    basin_name
) {

  ps_sub <- subset_samples(
    ps_object,
    Ephemeral_Basin == basin_name
  )

  lefse_res <- run_lefse_safe(
    ps_sub,
    "Phase"
  )

  marker_df <- microbiomeMarker::marker_table(
    lefse_res
  ) |>
    data.frame() |>
    dplyr::mutate(
      ef_lda = ifelse(
        enrich_group == "Wet",
        -abs(
          ef_lda
        ),
        abs(
          ef_lda
        )
      ),
      tax_short = clean_family_name(
        feature
      ),
      enrich_group = factor(
        enrich_group,
        levels = c(
          "Wet",
          "Dry"
        )
      ),
      Basin = basin_name
    )

  marker_df$tax_short <- factor(
    marker_df$tax_short,
    levels = marker_df$tax_short[
      order(
        marker_df$ef_lda,
        decreasing = TRUE
      )
    ]
  )

  max_abs <- ceiling(
    max(
      abs(
        marker_df$ef_lda
      ),
      na.rm = TRUE
    )
  )

  p <- ggplot(
    marker_df,
    aes(
      x = tax_short,
      y = ef_lda,
      fill = enrich_group
    )
  ) +
    geom_col(
      color = "black",
      alpha = 0.9
    ) +
    coord_flip() +
    theme_cowplot() +
    scale_fill_manual(
      values = phase_colors
    ) +
    scale_y_continuous(
      limits = c(
        -max_abs,
        max_abs
      )
    ) +
    labs(
      x = NULL,
      y = "LDA Score (log10)",
      fill = NULL,
      title = basin_name
    ) +
    theme(
      axis.text.y = element_text(
        size = 10
      ),
      panel.grid.major.x = element_line(
        color = "grey",
        linetype = "dotted"
      ),
      panel.grid.major.y = element_line(
        color = "grey",
        linetype = "dotted"
      )
    )

  list(
    result = lefse_res,
    table = marker_df,
    plot = p
  )
}

lefse_HCP <- run_phase_lefse_basin(
  ps,
  "Herradura Clay Pan"
)

lefse_CHF <- run_phase_lefse_basin(
  ps,
  "Ckoirama Halite Field"
)

lefse_YSB <- run_phase_lefse_basin(
  ps,
  "Yungay Station Basin"
)

p_HCP <- lefse_HCP$plot +
  theme(
    legend.position = "none"
  )

p_CHF <- lefse_CHF$plot +
  theme(
    legend.position = "none"
  )

p_YSB <- lefse_YSB$plot +
  theme(
    legend.position = "right"
  )

figure4c <- cowplot::plot_grid(
  p_HCP,
  p_CHF,
  p_YSB,
  ncol = 3,
  labels = c(
    "A",
    "B",
    "C"
  ),
  label_size = 15,
  rel_widths = c(
    1,
    1,
    1.25
  )
)

ggsave(
  file.path(
    figure_dir,
    "Figure4C_LEfSe_wet_dry_by_basin.tiff"
  ),
  plot = figure4c,
  width = 14,
  height = 7,
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure4C_LEfSe_wet_dry_by_basin.svg"
  ),
  plot = figure4c,
  width = 14,
  height = 7,
  device = "svg"
)

openxlsx::write.xlsx(
  list(
    HCP = lefse_HCP$table,
    CHF = lefse_CHF$table,
    YSB = lefse_YSB$table
  ),
  file.path(
    table_dir,
    "Figure4C_LEfSe_wet_dry_by_basin.xlsx"
  )
)


# ============================================================================ #
# Figure 4B | Basin biomarkers within Wet and Dry phases
# ============================================================================ #

run_basin_lefse_phase <- function(
    ps_object,
    phase_name
) {

  ps_sub <- subset_samples(
    ps_object,
    Phase == phase_name
  )

  lefse_res <- run_lefse_safe(
    ps_sub,
    "Ephemeral_Basin"
  )

  marker_df <- microbiomeMarker::marker_table(
    lefse_res
  ) |>
    data.frame() |>
    dplyr::mutate(
      tax_short = clean_family_name(
        feature
      ),
      enrich_group = factor(
        enrich_group,
        levels = basin_levels
      ),
      Phase = phase_name
    ) |>
    dplyr::arrange(
      enrich_group,
      ef_lda
    )

  marker_df$tax_short <- factor(
    marker_df$tax_short,
    levels = marker_df$tax_short
  )

  max_lda <- ceiling(
    max(
      marker_df$ef_lda,
      na.rm = TRUE
    )
  )

  p <- ggplot(
    marker_df,
    aes(
      x = tax_short,
      y = ef_lda,
      fill = enrich_group
    )
  ) +
    geom_col(
      color = "black",
      alpha = 0.9
    ) +
    coord_flip() +
    theme_cowplot() +
    scale_fill_manual(
      values = basin_colors,
      drop = FALSE
    ) +
    scale_y_continuous(
      limits = c(
        0,
        max_lda
      ),
      breaks = seq(
        0,
        max_lda,
        by = 1
      ),
      expand = expansion(
        mult = c(
          0,
          0.05
        )
      )
    ) +
    labs(
      x = NULL,
      y = "LDA Score (log10)",
      fill = NULL,
      title = paste(
        "Phase",
        phase_name
      )
    ) +
    theme(
      axis.text.y = element_text(
        size = 10
      ),
      panel.grid.major.x = element_line(
        color = "grey",
        linetype = "dotted"
      ),
      panel.grid.major.y = element_line(
        color = "grey",
        linetype = "dotted"
      )
    )

  list(
    result = lefse_res,
    table = marker_df,
    plot = p
  )
}

lefse_wet <- run_basin_lefse_phase(
  ps,
  "Wet"
)

lefse_dry <- run_basin_lefse_phase(
  ps,
  "Dry"
)

p_wet <- lefse_wet$plot +
  theme(
    legend.position = "none"
  )

p_dry <- lefse_dry$plot +
  theme(
    legend.position = "right"
  )

figure4b <- cowplot::plot_grid(
  p_wet,
  p_dry,
  ncol = 2,
  labels = c(
    "A",
    "B"
  ),
  label_size = 15,
  rel_widths = c(
    1,
    1.5
  )
)

ggsave(
  file.path(
    figure_dir,
    "Figure4B_LEfSe_basins_by_phase.tiff"
  ),
  plot = figure4b,
  width = 12,
  height = 7,
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure4B_LEfSe_basins_by_phase.svg"
  ),
  plot = figure4b,
  width = 12,
  height = 7,
  device = "svg"
)

openxlsx::write.xlsx(
  list(
    Wet = lefse_wet$table,
    Dry = lefse_dry$table
  ),
  file.path(
    table_dir,
    "Figure4B_LEfSe_basins_by_phase.xlsx"
  )
)


# ============================================================================ #
# Final summary
# ============================================================================ #

cat(
  "\nTaxonomic-composition analysis completed.\n"
)

cat(
  "Environmental samples: ",
  nsamples(ps),
  "\n",
  sep = ""
)

cat(
  "ASVs: ",
  ntaxa(ps),
  "\n",
  sep = ""
)

cat(
  "Displayed order categories (including Other/Unclassified): ",
  length(order_levels_rank),
  "\n",
  sep = ""
)

cat(
  "LEfSe rank: Family\n"
)

cat(
  "LEfSe LDA cutoff: ",
  lda_cutoff_lefse,
  "\n",
  sep = ""
)
