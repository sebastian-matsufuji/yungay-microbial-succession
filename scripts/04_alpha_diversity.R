# ============================================================================ #
# 04_alpha_diversity.R
# ============================================================================ #
#
# Purpose:
# Reproduce the alpha-diversity analyses used in the Yungay manuscript:
#
#   - repeated rarefaction with RTK to 359 reads
#   - observed ASV richness
#   - Shannon diversity
#   - Faith's phylogenetic diversity (PD)
#   - Figure 3A: temporal trajectories
#   - Figure 3B: basin comparisons
#   - Figure 3C: Wet vs Dry comparisons
#   - Supplementary Figure 3: alpha diversity vs physicochemical variables
#   - Supplementary Figure 4: within-basin differences among sampling points
#
# Notes:
# - Five samples below 359 reads are excluded before RTK rarefaction.
# - RTK uses 1,000 iterations and seed 2026, matching the original script.
# - Basin comparisons use Kruskal-Wallis + BH-adjusted Dunn tests, matching
#   the final Results/Figure 3 code.
# - Wet vs Dry comparisons use Wilcoxon rank-sum tests within each basin.
# - Supplementary Figure 3 is exported as four separate 3-panel figures
#   (Wc, pH, Salinity, ORP). These were combined externally in Illustrator
#   in the original figure assembly.
# - The final salinity panels retained GAM fits for all basin-metric
#   combinations, reproducing the final compiled figure/statistics.
#
# Input:
# - output_phyloseq_tree/ps_mdecon_v2_FTree.rds
#
# Main outputs:
# - ps_mdecon_v2_FTree_rtk359.rds
# - alpha_meta_rtk359.xlsx
# - Figure 3 component figures and statistics
# - Supplementary Figure 3 component figures and model statistics
# - Supplementary Figure 4 figure and statistics
#
# ============================================================================ #


# ============================================================================ #
# Packages
# ============================================================================ #

library(phyloseq)
library(rtk)
library(vegan)
library(picante)
library(dplyr)
library(tidyr)
library(tibble)
library(openxlsx)
library(ggplot2)
library(cowplot)
library(rstatix)
library(ggpubr)
library(ggforce)
library(ggh4x)
library(scales)
library(mgcv)
library(broom)
library(purrr)


# ============================================================================ #
# User paths
# ============================================================================ #

phyloseq_path <- file.path(
  "output_phyloseq_tree",
  "ps_mdecon_v2_FTree.rds"
)

output_dir <- "output_alpha_diversity"
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
object_dir <- file.path(output_dir, "objects")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(object_dir, showWarnings = FALSE, recursive = TRUE)


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
  "Wet" = "steelblue4",
  "Dry" = "#d8b365"
)

metric_labels <- c(
  "Sp_richness" = "Species richness",
  "Shannon" = "Shannon",
  "PD" = "Faith's PD"
)


# ============================================================================ #
# Load final phyloseq object
# ============================================================================ #

ps <- readRDS(phyloseq_path)

sample_data(ps)$Ephemeral_Basin <- factor(
  sample_data(ps)$Ephemeral_Basin,
  levels = basin_levels
)


# ============================================================================ #
# RTK repeated rarefaction to 359 reads
# ============================================================================ #

raref_depth <- 359
raref_repeats <- 1000
raref_seed <- 2026

# Exclude samples below the target depth.
ps_filt <- prune_samples(
  sample_sums(ps) >= raref_depth,
  ps
)

otu <- as(
  otu_table(ps_filt),
  "matrix"
)

if (!taxa_are_rows(ps_filt)) {
  otu <- t(otu)
}

rtk_res <- rtk(
  input = otu,
  depth = raref_depth,
  repeats = raref_repeats,
  ReturnMatrix = 1,
  margin = 2,
  seed = raref_seed,
  verbose = FALSE
)

rarefied_mat <- rtk_res$raremat

if (is.list(rarefied_mat)) {
  rarefied_mat <- rarefied_mat[[1]]
}

kept_samples <- colnames(rarefied_mat)
kept_taxa <- rownames(rarefied_mat)

ps_rtk <- phyloseq(
  otu_table(rarefied_mat, taxa_are_rows = TRUE),
  sample_data(
    sample_data(ps_filt)[kept_samples, , drop = FALSE]
  ),
  tax_table(
    tax_table(ps_filt)[kept_taxa, , drop = FALSE]
  ),
  phy_tree(
    prune_taxa(kept_taxa, ps_filt)
  )
)

# Remove taxa with zero total abundance after rarefaction.
ps_rtk <- prune_taxa(
  taxa_sums(ps_rtk) > 0,
  ps_rtk
)

saveRDS(
  ps_rtk,
  file.path(
    object_dir,
    "ps_mdecon_v2_FTree_rtk359.rds"
  )
)


# ============================================================================ #
# Calculate alpha-diversity metrics
# ============================================================================ #

comm <- as(
  otu_table(ps_rtk),
  "matrix"
)

if (taxa_are_rows(ps_rtk)) {
  comm <- t(comm)
}

# Samples in rows, ASVs in columns.
Sp_richness <- vegan::specnumber(comm)

Shannon <- vegan::diversity(
  comm,
  index = "shannon"
)

PD_index <- picante::pd(
  comm,
  phy_tree(ps_rtk),
  include.root = FALSE
)

alpha <- data.frame(
  Sp_richness = Sp_richness,
  Shannon = Shannon,
  PD = PD_index$PD,
  row.names = rownames(comm),
  check.names = FALSE
)

metadata <- data.frame(
  sample_data(ps_rtk),
  check.names = FALSE
)

metadata <- metadata[
  rownames(alpha),
  ,
  drop = FALSE
]

alpha_meta <- cbind(
  alpha,
  metadata
)

openxlsx::write.xlsx(
  alpha_meta,
  file.path(
    table_dir,
    "alpha_meta_rtk359.xlsx"
  ),
  rowNames = TRUE
)

saveRDS(
  alpha_meta,
  file.path(
    object_dir,
    "alpha_meta_rtk359.rds"
  )
)


# ============================================================================ #
# Date parser
# ============================================================================ #

parse_date_column <- function(x) {

  if (inherits(x, "Date")) {
    return(as.Date(x))
  }

  if (inherits(x, "POSIXt")) {
    return(as.Date(x))
  }

  x_chr <- as.character(x)
  x_num <- suppressWarnings(as.numeric(x_chr))

  n_valid <- sum(
    !is.na(x_chr) &
      x_chr != ""
  )

  n_num <- sum(
    !is.na(x_num)
  )

  if (
    n_valid > 0 &&
      n_num / n_valid >= 0.8 &&
      max(x_num, na.rm = TRUE) > 1000
  ) {
    return(
      as.Date(
        x_num,
        origin = "1899-12-30"
      )
    )
  }

  as.Date(x_chr)
}

alpha_meta$Date <- parse_date_column(
  alpha_meta$Date
)

alpha_meta$Ephemeral_Basin <- factor(
  alpha_meta$Ephemeral_Basin,
  levels = basin_levels
)

alpha_meta$Phase <- factor(
  alpha_meta$Phase,
  levels = c("Wet", "Dry")
)


# ============================================================================ #
# Figure 3A | Temporal alpha-diversity trajectories
# ============================================================================ #

alpha_meta_long <- alpha_meta |>
  tibble::rownames_to_column("Sample") |>
  tidyr::pivot_longer(
    cols = c(
      Sp_richness,
      Shannon,
      PD
    ),
    names_to = "Diversity_Index",
    values_to = "Value"
  ) |>
  dplyr::mutate(
    Diversity_Index = factor(
      Diversity_Index,
      levels = c(
        "Sp_richness",
        "Shannon",
        "PD"
      )
    ),
    Ephemeral_Basin = factor(
      Ephemeral_Basin,
      levels = basin_levels
    )
  ) |>
  dplyr::arrange(Date)

date_breaks <- as.Date(
  c(
    "2017-06-16",
    "2017-06-24",
    "2017-06-30",
    "2017-07-07",
    "2017-07-14",
    "2017-08-25",
    "2017-09-26",
    "2017-11-06",
    "2017-12-15",
    "2018-01-22",
    "2018-03-02",
    "2018-03-28",
    "2018-05-25",
    "2018-07-03",
    "2018-07-27",
    "2018-09-14"
  )
)

alpha_date_loess <- ggplot(
  alpha_meta_long,
  aes(
    x = Date,
    y = Value,
    colour = Ephemeral_Basin
  )
) +
  geom_smooth(
    aes(
      group = Ephemeral_Basin,
      fill = Ephemeral_Basin,
      color = Ephemeral_Basin
    ),
    method = "loess",
    se = TRUE,
    span = 0.5,
    alpha = 0.15,
    linewidth = 1
  ) +
  geom_point(
    aes(fill = Ephemeral_Basin),
    shape = 21,
    size = 2.2,
    color = "black"
  ) +
  facet_grid(
    rows = vars(Diversity_Index),
    scales = "free_y",
    space = "free_x",
    labeller = as_labeller(metric_labels)
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      Diversity_Index == "Sp_richness" ~ scale_y_continuous(
        breaks = c(0, 50, 100),
        limits = c(0, 100),
        expand = expansion(
          mult = c(0, 0.02)
        )
      ),
      Diversity_Index == "Shannon" ~ scale_y_continuous(
        breaks = c(0, 2, 4),
        limits = c(0, 5),
        expand = expansion(
          mult = c(0, 0.05)
        )
      ),
      Diversity_Index == "PD" ~ scale_y_continuous(
        breaks = c(0, 5, 10),
        limits = c(0, 12),
        expand = expansion(
          mult = c(0, 0.02)
        )
      )
    )
  ) +
  scale_color_manual(
    values = basin_colors
  ) +
  scale_fill_manual(
    values = basin_colors
  ) +
  scale_x_date(
    breaks = date_breaks,
    labels = scales::date_format("%d-%b"),
    expand = expansion(
      mult = 0.01
    )
  ) +
  theme_cowplot() +
  labs(
    x = NULL,
    y = "Diversity index value",
    title = NULL
  ) +
  theme(
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.5
    ),
    strip.background = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    ),
    legend.position = "bottom",
    legend.justification = "center",
    panel.grid.major.x = element_line(
      color = "grey15",
      linetype = "dotted"
    ),
    plot.margin = margin(
      t = 5.5,
      r = 5.5,
      b = 25,
      l = 5.5
    )
  ) +
  coord_cartesian(
    clip = "on"
  )

ggsave(
  plot = alpha_date_loess,
  filename = file.path(
    figure_dir,
    "Figure3A_alpha_date_loess.tiff"
  ),
  width = 10,
  height = 6,
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  plot = alpha_date_loess,
  filename = file.path(
    figure_dir,
    "Figure3A_alpha_date_loess.svg"
  ),
  width = 10,
  height = 6,
  device = "svg"
)


# ============================================================================ #
# Figure 3B | Alpha diversity among basins
# ============================================================================ #

run_kw_dunn <- function(
    df,
    response
) {

  df2 <- df |>
    dplyr::select(
      Ephemeral_Basin,
      dplyr::all_of(response)
    ) |>
    dplyr::rename(
      Value = !!response
    ) |>
    dplyr::filter(
      !is.na(Value),
      !is.na(Ephemeral_Basin)
    )

  kw <- rstatix::kruskal_test(
    df2,
    Value ~ Ephemeral_Basin
  )

  dunn <- rstatix::dunn_test(
    df2,
    Value ~ Ephemeral_Basin,
    p.adjust.method = "BH"
  )

  list(
    kw = kw,
    dunn = dunn
  )
}


make_basin_plot <- function(
    df,
    response,
    ylab
) {

  stats <- run_kw_dunn(
    df,
    response
  )

  dunn_plot <- stats$dunn |>
    rstatix::add_xy_position(
      x = "Ephemeral_Basin"
    ) |>
    dplyr::filter(
      p.adj <= 0.05
    )

  p <- ggplot(
    df,
    aes(
      x = Ephemeral_Basin,
      y = .data[[response]],
      fill = Ephemeral_Basin
    )
  ) +
    geom_violin(
      width = 0.92,
      alpha = 0.50,
      trim = FALSE,
      colour = "black",
      linewidth = 0.45
    ) +
    ggforce::geom_sina(
      shape = 21,
      size = 1.70,
      stroke = 0.25,
      colour = "black",
      alpha = 0.85,
      maxwidth = 0.34
    ) +
    geom_boxplot(
      width = 0.07,
      alpha = 1,
      colour = "black",
      linewidth = 0.45,
      outlier.shape = NA
    ) +
    scale_fill_manual(
      values = basin_colors
    ) +
    labs(
      x = NULL,
      y = ylab,
      subtitle = paste0(
        "Kruskal-Wallis, p = ",
        signif(
          stats$kw$p,
          3
        )
      )
    ) +
    theme_cowplot() +
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.30
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      plot.subtitle = element_text(
        size = 8
      )
    )

  if (nrow(dunn_plot) > 0) {
    p <- p +
      ggpubr::stat_pvalue_manual(
        dunn_plot,
        label = "p.adj.signif",
        hide.ns = TRUE,
        tip.length = 0.01,
        size = 4.2,
        bracket.size = 0.35
      )
  }

  list(
    plot = p,
    kw = stats$kw,
    dunn = stats$dunn
  )
}

basin_rich <- make_basin_plot(
  alpha_meta,
  "Sp_richness",
  "Species richness"
)

basin_shannon <- make_basin_plot(
  alpha_meta,
  "Shannon",
  "Shannon"
)

basin_pd <- make_basin_plot(
  alpha_meta,
  "PD",
  "Faith's PD"
)

legend_basin <- ggplot(
  alpha_meta,
  aes(
    x = Ephemeral_Basin,
    y = Sp_richness,
    fill = Ephemeral_Basin
  )
) +
  geom_violin() +
  scale_fill_manual(
    values = basin_colors,
    name = "Basin"
  ) +
  theme_void() +
  theme(
    legend.position = "bottom"
  )

legend_basin <- cowplot::get_legend(
  legend_basin
)

figure3b_grid <- cowplot::plot_grid(
  basin_rich$plot,
  basin_shannon$plot,
  basin_pd$plot,
  ncol = 3,
  align = "hv",
  axis = "tblr"
)

figure3b <- cowplot::plot_grid(
  figure3b_grid,
  legend_basin,
  ncol = 1,
  rel_heights = c(
    1,
    0.10
  )
)

ggsave(
  plot = figure3b,
  filename = file.path(
    figure_dir,
    "Figure3B_alpha_basins_kruskal_dunn.tiff"
  ),
  width = 6,
  height = 3,
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  plot = figure3b,
  filename = file.path(
    figure_dir,
    "Figure3B_alpha_basins_kruskal_dunn.svg"
  ),
  width = 6,
  height = 3,
  device = "svg"
)

kw_global_stats <- dplyr::bind_rows(
  basin_rich$kw |>
    mutate(Metric = "Sp_richness"),
  basin_shannon$kw |>
    mutate(Metric = "Shannon"),
  basin_pd$kw |>
    mutate(Metric = "PD")
)

kw_pairwise_stats <- dplyr::bind_rows(
  basin_rich$dunn |>
    mutate(Metric = "Sp_richness"),
  basin_shannon$dunn |>
    mutate(Metric = "Shannon"),
  basin_pd$dunn |>
    mutate(Metric = "PD")
)

openxlsx::write.xlsx(
  kw_global_stats,
  file.path(
    table_dir,
    "Figure3B_alpha_basins_kruskal_global.xlsx"
  )
)

openxlsx::write.xlsx(
  kw_pairwise_stats,
  file.path(
    table_dir,
    "Figure3B_alpha_basins_dunn_BH_pairwise.xlsx"
  )
)


# ============================================================================ #
# Figure 3C | Wet vs Dry alpha-diversity comparisons
# ============================================================================ #

alpha_phase_long <- alpha_meta |>
  tibble::rownames_to_column("Sample") |>
  tidyr::pivot_longer(
    cols = c(
      Sp_richness,
      Shannon,
      PD
    ),
    names_to = "Diversity_Index",
    values_to = "Value"
  ) |>
  dplyr::mutate(
    Phase = factor(
      Phase,
      levels = c(
        "Wet",
        "Dry"
      )
    ),
    Diversity_Index = factor(
      Diversity_Index,
      levels = c(
        "Sp_richness",
        "Shannon",
        "PD"
      )
    ),
    Ephemeral_Basin = factor(
      Ephemeral_Basin,
      levels = basin_levels
    )
  )

valid_phase_facets <- alpha_phase_long |>
  dplyr::filter(
    !is.na(Value),
    !is.na(Phase),
    !is.na(Ephemeral_Basin)
  ) |>
  dplyr::count(
    Diversity_Index,
    Ephemeral_Basin,
    Phase,
    name = "n"
  ) |>
  dplyr::group_by(
    Diversity_Index,
    Ephemeral_Basin
  ) |>
  dplyr::summarise(
    n_phases = n_distinct(Phase),
    min_n = min(n),
    .groups = "drop"
  ) |>
  dplyr::filter(
    n_phases == 2,
    min_n >= 2
  )

alpha_phase_stats_data <- alpha_phase_long |>
  dplyr::semi_join(
    valid_phase_facets,
    by = c(
      "Diversity_Index",
      "Ephemeral_Basin"
    )
  )

phase_stats <- alpha_phase_stats_data |>
  dplyr::group_by(
    Diversity_Index,
    Ephemeral_Basin
  ) |>
  rstatix::wilcox_test(
    Value ~ Phase
  ) |>
  rstatix::add_significance("p") |>
  dplyr::ungroup()

phase_stats_export <- phase_stats |>
  dplyr::mutate(
    dplyr::across(
      where(is.list),
      ~ vapply(
        .x,
        toString,
        character(1)
      )
    )
  ) |>
  as.data.frame()

write.table(
  phase_stats_export,
  file = file.path(
    table_dir,
    "Figure3C_alpha_phase_wilcox_results.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

phase_sig_labels <- phase_stats |>
  dplyr::left_join(
    alpha_phase_long |>
      dplyr::group_by(
        Diversity_Index,
        Ephemeral_Basin
      ) |>
      dplyr::summarise(
        y = max(
          Value,
          na.rm = TRUE
        ) +
          0.40 *
          diff(
            range(
              Value,
              na.rm = TRUE
            )
          ),
        .groups = "drop"
      ),
    by = c(
      "Diversity_Index",
      "Ephemeral_Basin"
    )
  ) |>
  dplyr::filter(
    p <= 0.05
  ) |>
  dplyr::mutate(
    Basin_x = as.numeric(
      Ephemeral_Basin
    )
  )

half_violin_max_width <- 0.36
split_box_offset <- 0.075

half_violin_data <- alpha_phase_long |>
  dplyr::filter(
    !is.na(Value),
    !is.na(Phase),
    !is.na(Ephemeral_Basin)
  ) |>
  dplyr::group_by(
    Diversity_Index,
    Ephemeral_Basin,
    Phase
  ) |>
  dplyr::group_modify(
    ~ {
      vals <- .x$Value
      center <- as.numeric(
        .y$Ephemeral_Basin
      )
      side <- if (
        as.character(
          .y$Phase
        ) == "Wet"
      ) {
        -1
      } else {
        1
      }

      if (
        length(vals) < 2 ||
          length(
            unique(vals)
          ) < 2
      ) {
        return(
          tibble(
            x = numeric(0),
            y = numeric(0)
          )
        )
      }

      dens <- density(
        vals,
        na.rm = TRUE,
        n = 256
      )

      dens_width <- dens$y /
        max(dens$y) *
        half_violin_max_width

      tibble(
        x = c(
          center,
          center +
            side *
            dens_width,
          center
        ),
        y = c(
          min(dens$x),
          dens$x,
          max(dens$x)
        )
      )
    }
  ) |>
  dplyr::ungroup()

split_box_data <- alpha_phase_long |>
  dplyr::filter(
    !is.na(Value),
    !is.na(Phase),
    !is.na(Ephemeral_Basin)
  ) |>
  dplyr::mutate(
    Basin_x = as.numeric(
      Ephemeral_Basin
    ),
    Box_x = Basin_x +
      ifelse(
        Phase == "Wet",
        -split_box_offset,
        split_box_offset
      )
  )

figure3c <- ggplot() +
  geom_polygon(
    data = half_violin_data,
    aes(
      x = x,
      y = y,
      group = interaction(
        Diversity_Index,
        Ephemeral_Basin,
        Phase
      ),
      fill = Phase
    ),
    alpha = 0.50,
    colour = "black",
    linewidth = 0.45
  ) +
  geom_boxplot(
    data = split_box_data,
    aes(
      x = Box_x,
      y = Value,
      group = interaction(
        Ephemeral_Basin,
        Phase
      ),
      fill = Phase
    ),
    width = 0.10,
    outlier.shape = NA,
    alpha = 1,
    colour = "black",
    linewidth = 0.45
  ) +
  geom_text(
    data = phase_sig_labels,
    aes(
      x = Basin_x,
      y = y,
      label = p.signif
    ),
    inherit.aes = FALSE,
    size = 4.2
  ) +
  facet_grid(
    rows = vars(Diversity_Index),
    scales = "free_y",
    labeller = as_labeller(
      metric_labels
    )
  ) +
  scale_x_continuous(
    breaks = 1:3,
    labels = basin_levels
  ) +
  scale_fill_manual(
    values = phase_colors,
    name = "Phase"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.02,
        0.18
      )
    )
  ) +
  labs(
    x = NULL,
    y = "Diversity index value"
  ) +
  theme_cowplot() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.30
    ),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    ),
    strip.background = element_blank(),
    legend.position = "right"
  ) +
  coord_cartesian(
    clip = "off"
  )

ggsave(
  plot = figure3c,
  filename = file.path(
    figure_dir,
    "Figure3C_alpha_phase_split_violin_wilcox.tiff"
  ),
  width = 5,
  height = 4,
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  plot = figure3c,
  filename = file.path(
    figure_dir,
    "Figure3C_alpha_phase_split_violin_wilcox.svg"
  ),
  width = 5,
  height = 4,
  device = "svg"
)


# ============================================================================ #
# Supplementary Figure 4 | Sampling-point heterogeneity
# ============================================================================ #

alpha_site_long <- alpha_meta |>
  tibble::rownames_to_column("Sample") |>
  tidyr::pivot_longer(
    cols = c(
      Sp_richness,
      Shannon,
      PD
    ),
    names_to = "Diversity_Index",
    values_to = "Value"
  ) |>
  dplyr::mutate(
    Diversity_Index = factor(
      Diversity_Index,
      levels = c(
        "Sp_richness",
        "Shannon",
        "PD"
      )
    ),
    Ephemeral_Basin = factor(
      Ephemeral_Basin,
      levels = basin_levels
    ),
    Sampling_Point = as.factor(
      Sampling_Point
    ),
    Basin_Site_fill = interaction(
      Ephemeral_Basin,
      Sampling_Point,
      sep = "_",
      drop = TRUE
    )
  )

facet_ranges <- alpha_site_long |>
  dplyr::group_by(
    Diversity_Index,
    Ephemeral_Basin
  ) |>
  dplyr::summarise(
    ymin = min(
      Value,
      na.rm = TRUE
    ),
    ymax = max(
      Value,
      na.rm = TRUE
    ),
    yrange = ymax - ymin,
    .groups = "drop"
  ) |>
  dplyr::mutate(
    yrange = ifelse(
      yrange == 0,
      1,
      yrange
    )
  )

site_kw <- alpha_site_long |>
  dplyr::group_by(
    Diversity_Index,
    Ephemeral_Basin
  ) |>
  rstatix::kruskal_test(
    Value ~ Sampling_Point
  ) |>
  dplyr::mutate(
    test = "Kruskal-Wallis",
    p_label = dplyr::case_when(
      p < 0.0001 ~ "KW, p < 0.0001",
      p < 0.001 ~ "KW, p < 0.001",
      TRUE ~ paste0(
        "KW, p = ",
        signif(
          p,
          3
        )
      )
    )
  )

site_dunn <- alpha_site_long |>
  dplyr::group_by(
    Diversity_Index,
    Ephemeral_Basin
  ) |>
  rstatix::dunn_test(
    Value ~ Sampling_Point,
    p.adjust.method = "BH"
  ) |>
  dplyr::left_join(
    facet_ranges,
    by = c(
      "Diversity_Index",
      "Ephemeral_Basin"
    )
  ) |>
  dplyr::group_by(
    Diversity_Index,
    Ephemeral_Basin
  ) |>
  dplyr::arrange(
    p.adj,
    .by_group = TRUE
  ) |>
  dplyr::mutate(
    y.position = ymax +
      yrange *
      (
        0.10 +
          0.09 *
          row_number()
      )
  ) |>
  dplyr::ungroup()

site_dunn_plot <- site_dunn |>
  dplyr::filter(
    !is.na(p.adj),
    p.adj <= 0.05
  )

site_kw_positions <- facet_ranges |>
  dplyr::mutate(
    x_kw = Inf,
    y_kw = Inf,
    hjust_kw = 1.05,
    vjust_kw = 1.2
  )

site_kw <- site_kw |>
  dplyr::left_join(
    site_kw_positions,
    by = c(
      "Diversity_Index",
      "Ephemeral_Basin"
    )
  )

color_levels <- alpha_site_long |>
  dplyr::distinct(
    Ephemeral_Basin,
    Sampling_Point,
    Basin_Site_fill
  ) |>
  dplyr::arrange(
    Ephemeral_Basin,
    Sampling_Point
  ) |>
  dplyr::pull(
    Basin_Site_fill
  )

site_colors <- c(
  "#A3E4D7",
  "#10A97B",
  "#00796B",
  "#FFD0A5",
  "#D94E00",
  "#9E3800",
  "#D8BFEA",
  "#7B4CB6",
  "#542F80"
)

site_colors <- setNames(
  site_colors,
  color_levels
)

supp4 <- ggplot(
  alpha_site_long,
  aes(
    x = Sampling_Point,
    y = Value
  )
) +
  geom_boxplot(
    aes(
      fill = Basin_Site_fill
    ),
    linewidth = 0.8,
    outlier.shape = NA,
    staplewidth = 0.55
  ) +
  facet_grid(
    rows = vars(Diversity_Index),
    cols = vars(Ephemeral_Basin),
    scales = "free",
    space = "free_x",
    labeller = as_labeller(
      c(
        metric_labels,
        setNames(
          basin_levels,
          basin_levels
        )
      )
    )
  ) +
  scale_fill_manual(
    values = site_colors
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.05,
        0.30
      )
    )
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL
  ) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.5
    ),
    legend.position = "none",
    strip.background = element_blank(),
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  ) +
  geom_text(
    data = site_kw,
    aes(
      x = x_kw,
      y = y_kw,
      label = p_label,
      hjust = hjust_kw,
      vjust = vjust_kw
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  ggpubr::stat_pvalue_manual(
    site_dunn_plot,
    hide.ns = TRUE,
    label = "p.adj.signif",
    size = 5,
    bracket.size = 0.35,
    tip.length = 0.01,
    step.increase = 0
  ) +
  coord_cartesian(
    clip = "off"
  )

ggsave(
  plot = supp4,
  filename = file.path(
    figure_dir,
    "Supplementary_Figure4_alpha_sampling_points.tiff"
  ),
  width = 8,
  height = 8,
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  plot = supp4,
  filename = file.path(
    figure_dir,
    "Supplementary_Figure4_alpha_sampling_points.svg"
  ),
  width = 8,
  height = 8,
  device = "svg"
)

openxlsx::write.xlsx(
  site_kw,
  file.path(
    table_dir,
    "Supplementary_Figure4_Kruskal_Wallis.xlsx"
  )
)

openxlsx::write.xlsx(
  site_dunn,
  file.path(
    table_dir,
    "Supplementary_Figure4_Dunn_BH.xlsx"
  )
)


# ============================================================================ #
# Supplementary Figure 3 | Model helper functions
# ============================================================================ #

fit_glm <- function(
    data,
    response,
    predictor
) {

  f <- reformulate(
    predictor,
    response = response
  )

  glm(
    f,
    data = data,
    family = gaussian(
      link = "identity"
    )
  )
}


safe_gam_k <- function(
    data,
    predictor,
    k = 4
) {

  k_safe <- min(
    k,
    dplyr::n_distinct(
      data[[predictor]]
    ) - 1,
    nrow(data) - 1
  )

  if (
    is.na(k_safe) ||
      k_safe < 3
  ) {
    return(
      NA_integer_
    )
  }

  as.integer(
    k_safe
  )
}


fit_gam_safe <- function(
    data,
    response,
    predictor,
    k = 4
) {

  k_safe <- safe_gam_k(
    data,
    predictor,
    k = k
  )

  if (
    is.na(k_safe)
  ) {
    stop(
      "Not enough unique predictor values to fit GAM safely."
    )
  }

  f <- stats::as.formula(
    paste0(
      response,
      " ~ s(",
      predictor,
      ", k = ",
      k_safe,
      ")"
    )
  )

  list(
    model = mgcv::gam(
      f,
      data = data,
      method = "REML"
    ),
    k_used = k_safe
  )
}


calc_glm_r2 <- function(
    model
) {

  if (
    is.na(
      model$null.deviance
    ) ||
      model$null.deviance == 0
  ) {
    return(
      NA_real_
    )
  }

  1 -
    (
      model$deviance /
        model$null.deviance
    )
}


add_delta_aic <- function(
    tab
) {

  if (
    all(
      is.na(
        tab$AIC
      )
    )
  ) {
    tab$delta_AIC <- NA_real_
  } else {
    tab$delta_AIC <- tab$AIC -
      min(
        tab$AIC,
        na.rm = TRUE
      )
  }

  tab
}


subset_model_data <- function(
    data,
    basin_name,
    response,
    predictor
) {

  data |>
    dplyr::select(
      Basin,
      all_of(
        c(
          response,
          predictor
        )
      )
    ) |>
    tidyr::drop_na() |>
    dplyr::filter(
      Basin == basin_name,
      is.finite(
        .data[[response]]
      ),
      is.finite(
        .data[[predictor]]
      )
    )
}


compare_glm_gam <- function(
    data,
    basin_name,
    response,
    predictor,
    k = 4
) {

  dat <- subset_model_data(
    data,
    basin_name,
    response,
    predictor
  )

  m_glm <- fit_glm(
    dat,
    response,
    predictor
  )

  gam_fit <- fit_gam_safe(
    dat,
    response,
    predictor,
    k = k
  )

  tibble(
    Basin = basin_name,
    Response = response,
    Predictor = predictor,
    Model = c(
      "GLM",
      "GAM"
    ),
    AIC = c(
      AIC(m_glm),
      AIC(gam_fit$model)
    ),
    k_used = c(
      NA_integer_,
      gam_fit$k_used
    )
  ) |>
    add_delta_aic()
}


extract_glm_gam_stats <- function(
    data,
    basin_name,
    response,
    predictor,
    k = 4
) {

  dat <- subset_model_data(
    data,
    basin_name,
    response,
    predictor
  )

  m_glm <- fit_glm(
    dat,
    response,
    predictor
  )

  glm_tidy <- broom::tidy(
    m_glm
  )

  glm_stats <- tibble(
    Basin = basin_name,
    Response = response,
    Predictor = predictor,
    Model = "GLM",
    N = nrow(dat),
    Slope = glm_tidy$estimate[
      glm_tidy$term == predictor
    ][1],
    EDF = NA_real_,
    R2 = calc_glm_r2(
      m_glm
    ),
    Adj_R2 = NA_real_,
    P_value = glm_tidy$p.value[
      glm_tidy$term == predictor
    ][1],
    AIC = AIC(
      m_glm
    ),
    k_used = NA_integer_,
    Note = NA_character_
  )

  gam_fit <- fit_gam_safe(
    dat,
    response,
    predictor,
    k = k
  )

  gam_summary <- summary(
    gam_fit$model
  )

  gam_stats <- tibble(
    Basin = basin_name,
    Response = response,
    Predictor = predictor,
    Model = "GAM",
    N = nrow(dat),
    Slope = NA_real_,
    EDF = gam_summary$s.table[
      1,
      "edf"
    ],
    R2 = gam_summary$r.sq,
    Adj_R2 = NA_real_,
    P_value = gam_summary$s.table[
      1,
      "p-value"
    ],
    AIC = AIC(
      gam_fit$model
    ),
    k_used = gam_fit$k_used,
    Note = NA_character_
  )

  dplyr::bind_rows(
    glm_stats,
    gam_stats
  ) |>
    add_delta_aic()
}


make_model_predictions <- function(
    data,
    basin_name,
    response,
    predictor,
    model_type,
    k = 4,
    n_grid = 100
) {

  dat <- subset_model_data(
    data,
    basin_name,
    response,
    predictor
  )

  model_fit <- if (
    model_type == "GLM"
  ) {
    fit_glm(
      dat,
      response,
      predictor
    )
  } else {
    fit_gam_safe(
      dat,
      response,
      predictor,
      k = k
    )$model
  }

  x_grid <- seq(
    min(
      dat[[predictor]],
      na.rm = TRUE
    ),
    max(
      dat[[predictor]],
      na.rm = TRUE
    ),
    length.out = n_grid
  )

  newdat <- tibble(
    !!predictor := x_grid
  )

  pred <- predict(
    model_fit,
    newdata = newdat,
    se.fit = TRUE
  )

  newdat |>
    dplyr::mutate(
      Basin = basin_name,
      Response = response,
      Predictor = predictor,
      Final_Model = model_type,
      fit = as.numeric(
        pred$fit
      ),
      se = as.numeric(
        pred$se.fit
      ),
      ymin = fit -
        1.96 *
        se,
      ymax = fit +
        1.96 *
        se
    )
}


format_p_value <- function(
    p
) {

  ifelse(
    is.na(p),
    "NA",
    format.pval(
      p,
      digits = 2,
      eps = 1e-3
    )
  )
}


make_stat_label <- function(
    model,
    slope,
    edf,
    r2,
    p_value
) {

  ifelse(
    model == "GLM",
    paste0(
      "GLM | slope=",
      formatC(
        slope,
        format = "e",
        digits = 2
      ),
      "  R²=",
      sprintf(
        "%.2f",
        r2
      ),
      "  p=",
      format_p_value(
        p_value
      )
    ),
    paste0(
      "GAM | edf=",
      sprintf(
        "%.2f",
        edf
      ),
      "  R²=",
      sprintf(
        "%.2f",
        r2
      ),
      "  p=",
      format_p_value(
        p_value
      )
    )
  )
}


# ============================================================================ #
# Supplementary Figure 3 | Final model settings
# ============================================================================ #

# These settings reproduce the four component figures used in the final
# Illustrator compilation.
#
# Wc, pH and ORP:
#   GLM vs GAM chosen by minimum AIC within each basin × metric.
#
# Salinity:
#   GAM retained for all basin × metric panels in the final compilation.

supp3_settings <- tibble::tribble(
  ~Predictor, ~Predictor_Label, ~Final_Model_Mode,
  "Wc", "Water content (%)", "best_aic",
  "pH", "pH", "best_aic",
  "Salinity", "Salinity (PSU)", "all_gam",
  "ORP", "ORP (mV)", "best_aic"
)

metric_table <- tibble::tribble(
  ~Response, ~Y_Label,
  "Rich_log1p", "Richness (log1p)",
  "Shannon", "Shannon",
  "PD", "Faith's PD"
)


# ============================================================================ #
# Supplementary Figure 3 | Analysis function
# ============================================================================ #

run_supp3_predictor <- function(
    predictor_var,
    predictor_label,
    final_model_mode
) {

  alpha_model <- alpha_meta |>
    tibble::rownames_to_column(
      "Sample"
    ) |>
    dplyr::mutate(
      Basin = factor(
        Ephemeral_Basin,
        levels = basin_levels
      ),
      Rich_log1p = log1p(
        Sp_richness
      ),
      !!predictor_var := suppressWarnings(
        as.numeric(
          .data[[predictor_var]]
        )
      )
    ) |>
    dplyr::filter(
      !is.na(Basin),
      !is.na(
        .data[[predictor_var]]
      ),
      is.finite(
        .data[[predictor_var]]
      )
    )

  layout <- tidyr::expand_grid(
    Basin = factor(
      basin_levels,
      levels = basin_levels
    ),
    Response = metric_table$Response
  )

  comparison_aic <- purrr::pmap_dfr(
    list(
      layout$Basin,
      layout$Response
    ),
    function(
      Basin,
      Response
    ) {
      compare_glm_gam(
        data = alpha_model,
        basin_name = Basin,
        response = Response,
        predictor = predictor_var,
        k = 4
      )
    }
  )

  comparison_stats <- purrr::pmap_dfr(
    list(
      layout$Basin,
      layout$Response
    ),
    function(
      Basin,
      Response
    ) {
      extract_glm_gam_stats(
        data = alpha_model,
        basin_name = Basin,
        response = Response,
        predictor = predictor_var,
        k = 4
      )
    }
  )

  if (
    final_model_mode == "best_aic"
  ) {

    final_model_table <- comparison_aic |>
      dplyr::group_by(
        Basin,
        Response
      ) |>
      dplyr::slice_min(
        AIC,
        n = 1,
        with_ties = FALSE
      ) |>
      dplyr::ungroup() |>
      dplyr::transmute(
        Basin,
        Response,
        Final_Model = Model
      )

  } else if (
    final_model_mode == "all_gam"
  ) {

    final_model_table <- layout |>
      dplyr::mutate(
        Final_Model = "GAM"
      )

  } else {

    final_model_table <- layout |>
      dplyr::mutate(
        Final_Model = "GLM"
      )
  }

  final_stats <- comparison_stats |>
    dplyr::left_join(
      final_model_table,
      by = c(
        "Basin",
        "Response"
      )
    ) |>
    dplyr::filter(
      Model == Final_Model
    ) |>
    dplyr::arrange(
      Basin,
      Response
    )

  prediction_data <- purrr::pmap_dfr(
    list(
      final_model_table$Basin,
      final_model_table$Response,
      final_model_table$Final_Model
    ),
    function(
      Basin,
      Response,
      Final_Model
    ) {
      make_model_predictions(
        data = alpha_model,
        basin_name = Basin,
        response = Response,
        predictor = predictor_var,
        model_type = Final_Model,
        k = 4,
        n_grid = 100
      )
    }
  )

  plot_combined_metric <- function(
      response,
      ylab
  ) {

    dat <- alpha_model |>
      dplyr::filter(
        !is.na(
          .data[[response]]
        ),
        !is.na(
          .data[[predictor_var]]
        )
      )

    pred_dat <- prediction_data |>
      dplyr::filter(
        Response == response
      )

    label_dat <- final_stats |>
      dplyr::filter(
        Response == response
      ) |>
      dplyr::mutate(
        Label = make_stat_label(
          Model,
          Slope,
          EDF,
          R2,
          P_value
        ),
        label_vjust = 1.25 +
          1.35 *
          (
            row_number() -
              1
          )
      )

    ggplot() +
      geom_point(
        data = dat,
        aes(
          x = .data[[predictor_var]],
          y = .data[[response]],
          color = Basin
        ),
        size = 1.5,
        alpha = 0.75
      ) +
      geom_ribbon(
        data = pred_dat,
        aes(
          x = .data[[predictor_var]],
          ymin = ymin,
          ymax = ymax,
          group = Basin,
          fill = Basin
        ),
        alpha = 0.13,
        colour = NA
      ) +
      geom_line(
        data = pred_dat,
        aes(
          x = .data[[predictor_var]],
          y = fit,
          color = Basin,
          group = Basin
        ),
        linewidth = 0.8
      ) +
      geom_text(
        data = label_dat,
        aes(
          x = Inf,
          y = Inf,
          label = Label,
          color = Basin,
          vjust = label_vjust
        ),
        hjust = 1.02,
        size = 2.7,
        show.legend = FALSE,
        inherit.aes = FALSE
      ) +
      scale_color_manual(
        values = basin_colors,
        drop = FALSE,
        name = "Basin"
      ) +
      scale_fill_manual(
        values = basin_colors,
        drop = FALSE,
        guide = "none"
      ) +
      labs(
        x = predictor_label,
        y = ylab
      ) +
      theme_cowplot() +
      theme(
        panel.border = element_rect(
          fill = NA,
          colour = "black",
          linewidth = 0.7
        ),
        legend.position = "right",
        legend.box = "horizontal",
        legend.title = element_text(
          size = 10
        ),
        legend.text = element_text(
          size = 9
        ),
        legend.key.width = grid::unit(
          1.7,
          "cm"
        )
      )
  }

  plot_list <- purrr::pmap(
    list(
      metric_table$Response,
      metric_table$Y_Label
    ),
    function(
      Response,
      Y_Label
    ) {
      plot_combined_metric(
        response = Response,
        ylab = Y_Label
      )
    }
  )

  shared_legend <- cowplot::get_legend(
    plot_list[[1]] +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center"
      )
  )

  plot_list_nolegend <- purrr::map(
    plot_list,
    ~ .x +
      theme(
        legend.position = "none"
      )
  )

  panels <- cowplot::plot_grid(
    plotlist = plot_list_nolegend,
    ncol = 1,
    labels = c(
      "A",
      "B",
      "C"
    ),
    label_size = 12
  )

  final_figure <- cowplot::plot_grid(
    panels,
    shared_legend,
    ncol = 1,
    rel_heights = c(
      1,
      0.07
    )
  )

  openxlsx::write.xlsx(
    comparison_aic,
    file.path(
      table_dir,
      paste0(
        "compare_",
        predictor_var,
        "_combined_single_glm_vs_gam_aic.xlsx"
      )
    )
  )

  openxlsx::write.xlsx(
    comparison_stats,
    file.path(
      table_dir,
      paste0(
        "compare_",
        predictor_var,
        "_combined_single_glm_vs_gam_stats.xlsx"
      )
    )
  )

  openxlsx::write.xlsx(
    final_stats,
    file.path(
      table_dir,
      paste0(
        "final_",
        predictor_var,
        "_combined_single_chosen_model_stats.xlsx"
      )
    )
  )

  ggsave(
    plot = final_figure,
    filename = file.path(
      figure_dir,
      paste0(
        "final_",
        predictor_var,
        "_combined_single_3panels.tiff"
      )
    ),
    width = 10,
    height = 4,
    dpi = 300,
    device = "tiff",
    compression = "lzw"
  )

  ggsave(
    plot = final_figure,
    filename = file.path(
      figure_dir,
      paste0(
        "final_",
        predictor_var,
        "_combined_single_3panels.svg"
      )
    ),
    width = 10,
    height = 4,
    device = "svg"
  )

  final_stats
}


# ============================================================================ #
# Supplementary Figure 3 | Run four predictors
# ============================================================================ #

supp3_final_stats <- purrr::pmap_dfr(
  supp3_settings,
  function(
    Predictor,
    Predictor_Label,
    Final_Model_Mode
  ) {
    run_supp3_predictor(
      predictor_var = Predictor,
      predictor_label = Predictor_Label,
      final_model_mode = Final_Model_Mode
    )
  }
)

openxlsx::write.xlsx(
  supp3_final_stats,
  file.path(
    table_dir,
    "compilation_alpha_models.xlsx"
  )
)


# ============================================================================ #
# Final summary
# ============================================================================ #

cat(
  "\nAlpha-diversity analysis completed.\n"
)

cat(
  "Samples before RTK filtering: ",
  nsamples(ps),
  "\n",
  sep = ""
)

cat(
  "Samples retained at >=359 reads: ",
  nsamples(ps_rtk),
  "\n",
  sep = ""
)

cat(
  "ASVs retained after RTK rarefaction: ",
  ntaxa(ps_rtk),
  "\n",
  sep = ""
)

cat(
  "Supplementary Figure 3 component figures written for Wc, pH, Salinity, and ORP.\n"
)
