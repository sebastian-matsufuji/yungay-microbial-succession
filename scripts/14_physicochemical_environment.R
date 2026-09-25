# ============================================================================ #
# 14_physicochemical_environment.R
# ============================================================================ #
#
# Purpose:
# Reproduce the environmental analyses supporting Figure 2 and
# Supplementary Figure 2.
#
# Main components:
#
#   Figure 2A-D
#     Temporal trajectories of water content, salinity, ORP, and pH.
#
#   Figure 2E-F
#     Wet- and dry-phase basin comparisons using Kruskal-Wallis tests
#     followed by Holm-adjusted Dunn tests.
#
#   Figure 2G
#     Temporal total organic carbon (TOC) in water and sediment.
#
#   Supplementary Figure 2
#     Pearson correlation matrices for wet geochemistry, dry geochemistry,
#     and the core physicochemical variables.
#
# Important:
# - The curated script uses metadata stored in the final phyloseq object.
#   Historical code read the same final metadata from metadata.xlsx.
# - Controls are already absent from the final phyloseq object.
# - Only the four physicochemical variables retained in the manuscript
#   Figure 2 are included in the phase-comparison panels. Historical
#   exploratory Wet-phase plots also included TDS and temperature.
# - The Wet/Dry basin comparisons are descriptive basin-level comparisons
#   across each complete phase, not date-matched tests.
#
# Input:
# - output_phyloseq_tree/ps_mdecon_v2_FTree.rds
#
# Main outputs:
# - Figure2_A-D_physicochemical_time.*
# - Figure2_E_wet_basin_comparisons.*
# - Figure2_F_dry_basin_comparisons.*
# - Figure2_G_TOC_water_sediment.*
# - Supplementary_Figure2_environmental_correlations.*
# - corresponding statistical tables
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
library(patchwork)
library(scales)
library(rstatix)
library(ggpubr)
library(ggforce)
library(caret)
library(usdm)
library(openxlsx)


# ============================================================================ #
# Paths and reproducibility
# ============================================================================ #

set.seed(2026)

phyloseq_path <- file.path(
  "output_phyloseq_tree",
  "ps_mdecon_v2_FTree.rds"
)

output_dir <- "output_environment"

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
# Load final metadata
# ============================================================================ #

ps <- readRDS(
  phyloseq_path
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
# Date parser
# ============================================================================ #

parse_date_column <- function(
    x
) {

  if (
    inherits(
      x,
      "Date"
    )
  ) {
    return(
      as.Date(
        x
      )
    )
  }

  if (
    inherits(
      x,
      "POSIXt"
    )
  ) {
    return(
      as.Date(
        x
      )
    )
  }

  x_chr <- as.character(
    x
  )

  x_num <- suppressWarnings(
    as.numeric(
      x_chr
    )
  )

  valid <- !is.na(
    x_chr
  ) &
    x_chr != ""

  if (
    any(
      valid
    ) &&
      sum(
        !is.na(
          x_num
        )
      ) /
        sum(
          valid
        ) >=
        0.8 &&
      max(
        x_num,
        na.rm = TRUE
      ) >
        1000
  ) {

    return(
      as.Date(
        x_num,
        origin = "1899-12-30"
      )
    )
  }

  as.Date(
    x_chr
  )
}

metadata$Date <- parse_date_column(
  metadata$Date
)


# ============================================================================ #
# Factor order
# ============================================================================ #

basin_levels <- c(
  "Herradura Clay Pan",
  "Ckoirama Halite Field",
  "Yungay Station Basin"
)

basin_codes <- c(
  "HCP",
  "CHF",
  "YSB"
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

metadata$Ephemeral_Basin <- factor(
  as.character(
    metadata$Ephemeral_Basin
  ),
  levels = basin_levels
)

metadata$Basin <- factor(
  as.character(
    metadata$Basin
  ),
  levels = basin_codes
)

metadata$Sampling_Point <- factor(
  as.character(
    metadata$Sampling_Point
  ),
  levels = sampling_point_levels
)

metadata$Time <- factor(
  as.character(
    metadata$Time
  ),
  levels = paste0(
    "T",
    1:20
  )
)

metadata$Phase <- factor(
  as.character(
    metadata$Phase
  ),
  levels = c(
    "Wet",
    "Dry"
  )
)

metadata <- metadata |>
  arrange(
    Date
  )


# ============================================================================ #
# Validate required variables
# ============================================================================ #

required_columns <- c(
  "Date",
  "Time",
  "Phase",
  "Basin",
  "Ephemeral_Basin",
  "Sampling_Point",
  "Matrix",
  "Wc",
  "Salinity",
  "ORP",
  "pH",
  "TOC"
)

missing_required <- setdiff(
  required_columns,
  colnames(
    metadata
  )
)

if (
  length(
    missing_required
  ) >
    0
) {
  stop(
    "Missing required metadata columns: ",
    paste(
      missing_required,
      collapse = ", "
    )
  )
}


# ============================================================================ #
# Shared colors
# ============================================================================ #

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
  "HCP" = "#1b9e77",
  "CHF" = "#d95f02",
  "YSB" = "#7570b3"
)


# ============================================================================ #
# Figure 2A-D | Temporal physicochemical trajectories
# ============================================================================ #

loess_span <- 0.50
line_width <- 0.75
line_alpha <- 0.95

point_size <- 2.3
point_alpha <- 0.95
point_stroke <- 0.25
point_jitter_days <- 0.60

shade_fill <- "grey75"
shade_alpha <- 0.25

panel_border_width <- 0.55
grid_color <- "grey75"
grid_linewidth <- 0.40
grid_linetype <- "dotted"

facet_title_size <- 9.5
axis_title_size <- 10.5
axis_text_size <- 8.5
legend_title_size <- 9
legend_text_size <- 8.5

surface_water_periods <- data.frame(
  Ephemeral_Basin = factor(
    basin_levels,
    levels = basin_levels
  ),
  xmin = as.Date(
    c(
      "2017-06-16",
      "2017-06-16",
      "2017-06-16"
    )
  ),
  xmax = as.Date(
    c(
      "2017-06-20",
      "2017-07-14",
      "2017-07-14"
    )
  )
)

date_breaks <- as.Date(
  c(
    "2017-06-16",
    "2017-06-30",
    "2017-07-14",
    "2017-08-25",
    "2017-11-06",
    "2018-01-22",
    "2018-03-28",
    "2018-05-25",
    "2018-07-27",
    "2018-09-14"
  )
)


# ---------------------------------------------------------------------------- #
# Helper | Temporal parameter plot
# ---------------------------------------------------------------------------- #

make_parameter_plot <- function(
    data,
    variable,
    y_label,
    show_x_text = TRUE,
    show_x_title = FALSE,
    show_strip_text = TRUE
) {

  plot_data <- data |>
    filter(
      !is.na(
        .data[
          [
            variable
          ]
        ]
      ),
      !is.na(
        Date
      ),
      !is.na(
        Ephemeral_Basin
      ),
      !is.na(
        Sampling_Point
      )
    )

  strip_text_color <- if (
    show_strip_text
  ) {
    "black"
  } else {
    "transparent"
  }

  p <- ggplot(
    plot_data,
    aes(
      x = Date
    )
  ) +
    geom_rect(
      data = surface_water_periods,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = shade_fill,
      alpha = shade_alpha
    ) +
    geom_smooth(
      aes(
        y = .data[
          [
            variable
          ]
        ],
        colour = Sampling_Point,
        group = Sampling_Point
      ),
      method = "loess",
      se = FALSE,
      span = loess_span,
      linewidth = line_width,
      alpha = line_alpha
    ) +
    geom_point(
      aes(
        y = .data[
          [
            variable
          ]
        ],
        fill = Sampling_Point
      ),
      shape = 21,
      size = point_size,
      stroke = point_stroke,
      alpha = point_alpha,
      colour = "black",
      position = position_jitter(
        width = point_jitter_days,
        height = 0
      )
    ) +
    facet_wrap(
      ~ Ephemeral_Basin,
      ncol = 1,
      scales = "fixed"
    ) +
    scale_color_manual(
      values = sampling_point_colors,
      guide = "none"
    ) +
    scale_fill_manual(
      values = sampling_point_colors,
      name = "Sampling point",
      drop = FALSE
    ) +
    scale_x_date(
      breaks = date_breaks,
      labels = date_format(
        "%d-%b"
      ),
      expand = expansion(
        mult = c(
          0.01,
          0.01
        )
      )
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0.05,
          0.08
        )
      )
    ) +
    labs(
      x = if (
        show_x_title
      ) {
        "Date"
      } else {
        NULL
      },
      y = y_label
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    theme_cowplot() +
    theme(
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = panel_border_width
      ),
      panel.grid.major.x = element_line(
        colour = grid_color,
        linetype = grid_linetype,
        linewidth = grid_linewidth
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(
        size = facet_title_size,
        face = "bold",
        colour = strip_text_color,
        margin = margin(
          t = 1,
          r = 0,
          b = 2,
          l = 0
        )
      ),
      panel.spacing.y = grid::unit(
        0.08,
        "cm"
      ),
      axis.title.x = element_text(
        size = axis_title_size,
        margin = margin(
          t = 3
        )
      ),
      axis.title.y = element_text(
        size = axis_title_size,
        margin = margin(
          r = 3
        )
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = axis_text_size
      ),
      axis.text.y = element_text(
        size = axis_text_size
      ),
      axis.ticks = element_line(
        linewidth = 0.35
      ),
      legend.position = "right",
      legend.title = element_text(
        size = legend_title_size
      ),
      legend.text = element_text(
        size = legend_text_size
      ),
      legend.key.width = grid::unit(
        0.55,
        "cm"
      ),
      legend.key.height = grid::unit(
        0.35,
        "cm"
      ),
      plot.margin = margin(
        t = 2,
        r = 3,
        b = 2,
        l = 3
      )
    )

  if (
    !show_x_text
  ) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  p
}


# ---------------------------------------------------------------------------- #
# Build Figure 2A-D
# ---------------------------------------------------------------------------- #

p_Wc <- make_parameter_plot(
  metadata,
  "Wc",
  "Water content (%)",
  show_x_text = FALSE,
  show_x_title = FALSE,
  show_strip_text = TRUE
)

p_Salinity <- make_parameter_plot(
  metadata,
  "Salinity",
  "Salinity (PSU)",
  show_x_text = FALSE,
  show_x_title = FALSE,
  show_strip_text = FALSE
)

p_ORP <- make_parameter_plot(
  metadata,
  "ORP",
  "Redox potential (mV)",
  show_x_text = TRUE,
  show_x_title = FALSE,
  show_strip_text = TRUE
)

p_pH <- make_parameter_plot(
  metadata,
  "pH",
  "pH",
  show_x_text = TRUE,
  show_x_title = TRUE,
  show_strip_text = FALSE
)

figure2_ad <- (
  (
    p_Wc |
      p_Salinity
  ) /
    (
      p_ORP |
        p_pH
    )
) +
  plot_layout(
    guides = "collect",
    widths = c(
      1,
      1
    ),
    heights = c(
      1,
      1
    )
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "right",
    plot.tag = element_text(
      size = 11,
      face = "bold"
    )
  )

ggsave(
  file.path(
    figure_dir,
    "Figure2_A-D_physicochemical_time.tiff"
  ),
  plot = figure2_ad,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure2_A-D_physicochemical_time.svg"
  ),
  plot = figure2_ad,
  width = 12,
  height = 8,
  units = "in",
  device = "svg"
)


# ============================================================================ #
# Figure 2E-F | Wet/Dry basin comparisons
# ============================================================================ #

core_variables <- tribble(
  ~Variable, ~Y_label,
  "Wc", "Water content (%)",
  "pH", "pH",
  "ORP", "ORP (mV)",
  "Salinity", "Salinity (PSU)"
)

wet_manual_positions <- list(
  Wc = c(
    100,
    104,
    95
  ),
  pH = c(
    8.6,
    8.8
  ),
  ORP = c(
    255,
    260
  ),
  Salinity = c(
    410,
    450,
    430
  )
)

dry_manual_positions <- list(
  Wc = c(
    23,
    25
  ),
  pH = c(
    8.75,
    8.8
  ),
  ORP = NULL,
  Salinity = c(
    105,
    110
  )
)


# ---------------------------------------------------------------------------- #
# Helper | Phase comparison
# ---------------------------------------------------------------------------- #

run_phase_comparison <- function(
    data,
    phase_name,
    variable,
    y_label,
    manual_positions = NULL
) {

  phase_data <- data |>
    filter(
      Phase ==
        phase_name,
      !is.na(
        .data[
          [
            variable
          ]
        ]
      )
    ) |>
    mutate(
      Basin = factor(
        Basin,
        levels = basin_codes
      )
    )

  formula <- as.formula(
    paste(
      variable,
      "~ Basin"
    )
  )

  kw <- phase_data |>
    kruskal_test(
      formula
    )

  dunn <- phase_data |>
    dunn_test(
      formula,
      p.adjust.method = "holm"
    )

  dunn_sig <- dunn |>
    filter(
      p.adj.signif !=
        "ns"
    )

  if (
    nrow(
      dunn_sig
    ) >
      0
  ) {

    if (
      !is.null(
        manual_positions
      )
    ) {

      if (
        length(
          manual_positions
        ) !=
          nrow(
            dunn_sig
          )
      ) {
        stop(
          "Manual annotation positions do not match significant Dunn ",
          "comparisons for ",
          phase_name,
          " / ",
          variable,
          "."
        )
      }

      dunn_sig$y.position <-
        manual_positions

    } else {

      y_max <- max(
        phase_data[
          [
            variable
          ]
        ],
        na.rm = TRUE
      )

      step <- max(
        abs(
          y_max
        ) *
          0.06,
        diff(
          range(
            phase_data[
              [
                variable
              ]
            ],
            na.rm = TRUE
          )
        ) *
          0.06
      )

      dunn_sig$y.position <- seq(
        y_max +
          step,
        by = step,
        length.out =
          nrow(
            dunn_sig
          )
      )
    }
  }

  kw_label <- paste0(
    "\u03c7\u00b2(",
    kw$df,
    ") = ",
    sprintf(
      "%.2f",
      kw$statistic
    ),
    ", ",
    ifelse(
      kw$p <
        0.0001,
      "p < 0.0001",
      paste0(
        "p = ",
        formatC(
          kw$p,
          format = "f",
          digits = 4
        )
      )
    )
  )

  p <- ggplot(
    phase_data,
    aes(
      x = Basin,
      y = .data[
        [
          variable
        ]
      ],
      fill = Basin
    )
  ) +
    geom_boxplot(
      linewidth = 0.45,
      outlier.shape = NA
    ) +
    geom_jitter(
      shape = 21,
      width = 0.16,
      height = 0,
      size = 1.55,
      stroke = 0.25,
      alpha = 0.85,
      colour = "black"
    ) +
    scale_fill_manual(
      values = basin_colors
    ) +
    labs(
      x = NULL,
      y = y_label,
      subtitle = kw_label
    ) +
    theme_cowplot() +
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.55
      ),
      axis.title = element_text(
        size = 10.5
      ),
      axis.text = element_text(
        size = 9.5
      ),
      plot.subtitle = element_text(
        size = 9.7,
        margin = margin(
          b = 3
        )
      ),
      legend.position = "none",
      plot.margin = margin(
        4,
        4,
        4,
        4
      )
    )

  if (
    nrow(
      dunn_sig
    ) >
      0
  ) {
    p <- p +
      stat_pvalue_manual(
        dunn_sig,
        hide.ns = TRUE,
        label = "p.adj.signif",
        size = 4.2,
        bracket.size = 0.35
      )
  }

  list(
    plot = p,
    kruskal = kw,
    dunn = dunn
  )
}


# ---------------------------------------------------------------------------- #
# Wet phase
# ---------------------------------------------------------------------------- #

wet_results <- list()

for (
  i in seq_len(
    nrow(
      core_variables
    )
  )
) {

  variable_i <- core_variables$Variable[
    i
  ]

  wet_results[
    [
      variable_i
    ]
  ] <- run_phase_comparison(
    data = metadata,
    phase_name = "Wet",
    variable = variable_i,
    y_label = core_variables$Y_label[
      i
    ],
    manual_positions =
      wet_manual_positions[
        [
          variable_i
        ]
      ]
  )
}

figure2_e <- (
  wet_results$Wc$plot |
    wet_results$pH$plot
) /
  (
    wet_results$ORP$plot |
      wet_results$Salinity$plot
  ) +
  plot_annotation(
    title = "Wet phase"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure2_E_wet_basin_comparisons.tiff"
  ),
  plot = figure2_e,
  width = 8,
  height = 7,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure2_E_wet_basin_comparisons.svg"
  ),
  plot = figure2_e,
  width = 8,
  height = 7,
  units = "in"
)


# ---------------------------------------------------------------------------- #
# Dry phase
# ---------------------------------------------------------------------------- #

dry_results <- list()

for (
  i in seq_len(
    nrow(
      core_variables
    )
  )
) {

  variable_i <- core_variables$Variable[
    i
  ]

  dry_results[
    [
      variable_i
    ]
  ] <- run_phase_comparison(
    data = metadata,
    phase_name = "Dry",
    variable = variable_i,
    y_label = core_variables$Y_label[
      i
    ],
    manual_positions =
      dry_manual_positions[
        [
          variable_i
        ]
      ]
  )
}

figure2_f <- (
  dry_results$Wc$plot |
    dry_results$pH$plot |
    dry_results$ORP$plot |
    dry_results$Salinity$plot
) +
  plot_layout(
    ncol = 4
  ) +
  plot_annotation(
    title = "Dry phase"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure2_F_dry_basin_comparisons.tiff"
  ),
  plot = figure2_f,
  width = 12,
  height = 4,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure2_F_dry_basin_comparisons.svg"
  ),
  plot = figure2_f,
  width = 12,
  height = 4,
  units = "in"
)


# ---------------------------------------------------------------------------- #
# Export Wet/Dry statistics
# ---------------------------------------------------------------------------- #

wet_kw <- bind_rows(
  lapply(
    names(
      wet_results
    ),
    function(v) {
      wet_results[
        [
          v
        ]
      ]$kruskal |>
        mutate(
          Variable = v,
          .before = 1
        )
    }
  )
)

wet_dunn <- bind_rows(
  lapply(
    names(
      wet_results
    ),
    function(v) {
      wet_results[
        [
          v
        ]
      ]$dunn |>
        mutate(
          Variable = v,
          .before = 1
        )
    }
  )
)

dry_kw <- bind_rows(
  lapply(
    names(
      dry_results
    ),
    function(v) {
      dry_results[
        [
          v
        ]
      ]$kruskal |>
        mutate(
          Variable = v,
          .before = 1
        )
    }
  )
)

dry_dunn <- bind_rows(
  lapply(
    names(
      dry_results
    ),
    function(v) {
      dry_results[
        [
          v
        ]
      ]$dunn |>
        mutate(
          Variable = v,
          .before = 1
        )
    }
  )
)

openxlsx::write.xlsx(
  list(
    Wet_Kruskal_Wallis =
      wet_kw,
    Wet_Dunn_Holm =
      wet_dunn,
    Dry_Kruskal_Wallis =
      dry_kw,
    Dry_Dunn_Holm =
      dry_dunn
  ),
  file.path(
    table_dir,
    "Figure2_E-F_environmental_basin_statistics.xlsx"
  ),
  overwrite = TRUE
)


# ============================================================================ #
# Figure 2G | TOC in water and sediment
# ============================================================================ #

toc_points <- c(
  "HCP1",
  "CHF2",
  "YSB3"
)

toc_colors <- c(
  "HCP1" = "#1b9e77",
  "CHF2" = "#d95f02",
  "YSB3" = "#7570b3"
)

toc_metadata <- metadata |>
  filter(
    Sampling_Point %in%
      toc_points
  ) |>
  mutate(
    Sampling_Point = factor(
      as.character(
        Sampling_Point
      ),
      levels = toc_points
    )
  )

water_toc <- toc_metadata |>
  filter(
    Matrix ==
      "Water",
    !is.na(
      Date
    ),
    !is.na(
      TOC
    )
  ) |>
  mutate(
    Matrix = "Water",
    TOC_gL =
      TOC /
      1000
  )

sediment_toc <- toc_metadata |>
  filter(
    Matrix ==
      "Sediment",
    !is.na(
      Date
    ),
    !is.na(
      TOC
    )
  ) |>
  mutate(
    Matrix = "Sediment"
  )

if (
  nrow(
    water_toc
  ) ==
    0 ||
    nrow(
      sediment_toc
    ) ==
    0
) {
  stop(
    "Water and/or sediment TOC data were not found."
  )
}

toc_scale_factor <- max(
  water_toc$TOC_gL,
  na.rm = TRUE
) /
  max(
    sediment_toc$TOC,
    na.rm = TRUE
  )

toc_date_breaks <- as.Date(
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

toc_colors_sediment <- scales::alpha(
  toc_colors,
  0.45
)

names(
  toc_colors_sediment
) <- paste0(
  names(
    toc_colors_sediment
  ),
  "_Sediment"
)

figure2_g <- ggplot() +
  geom_area(
    data = sediment_toc,
    aes(
      x = Date,
      y = TOC *
        toc_scale_factor,
      fill = Sampling_Point,
      group = Sampling_Point
    ),
    alpha = 0.4,
    position = "identity"
  ) +
  geom_line(
    data = sediment_toc,
    aes(
      x = Date,
      y = TOC *
        toc_scale_factor,
      color = paste0(
        Sampling_Point,
        "_Sediment"
      ),
      group = Sampling_Point
    ),
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_area(
    data = water_toc,
    aes(
      x = Date,
      y = TOC_gL,
      fill = Sampling_Point,
      group = Sampling_Point
    ),
    alpha = 0.4,
    position = "identity"
  ) +
  geom_line(
    data = water_toc,
    aes(
      x = Date,
      y = TOC_gL,
      color = Sampling_Point,
      group = Sampling_Point
    ),
    linewidth = 1
  ) +
  geom_point(
    data = sediment_toc,
    aes(
      x = Date,
      y = TOC *
        toc_scale_factor,
      fill = Sampling_Point,
      shape = Matrix
    ),
    size = 2.5,
    stroke = 0.9,
    color = "black"
  ) +
  geom_point(
    data = water_toc,
    aes(
      x = Date,
      y = TOC_gL,
      fill = Sampling_Point,
      shape = Matrix
    ),
    size = 2.5,
    stroke = 0.9,
    color = "black"
  ) +
  scale_y_continuous(
    name =
      "Total organic carbon (water, g/L)",
    sec.axis = sec_axis(
      ~ . /
        toc_scale_factor,
      name =
        "Total organic carbon (sediment, %)"
    )
  ) +
  scale_x_date(
    breaks = toc_date_breaks,
    labels = scales::label_date(
      "%d %b"
    ),
    expand = expansion(
      mult = 0.01
    )
  ) +
  scale_color_manual(
    values = c(
      toc_colors,
      toc_colors_sediment
    ),
    breaks = toc_points
  ) +
  scale_fill_manual(
    values = toc_colors,
    breaks = toc_points
  ) +
  scale_shape_manual(
    values = c(
      "Water" = 21,
      "Sediment" = 24
    ),
    breaks = c(
      "Water",
      "Sediment"
    )
  ) +
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(
        fill = c(
          "black",
          "grey35"
        ),
        color = c(
          "black",
          "grey35"
        )
      )
    ),
    color = guide_legend(
      order = 2,
      override.aes = list(
        shape = 21,
        fill = toc_colors
      )
    ),
    fill = "none"
  ) +
  theme_cowplot() +
  theme(
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    ),
    legend.position = "right",
    panel.grid.major.x = element_line(
      color = "grey15",
      linetype = "dotted"
    ),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.box.background = element_blank(),
    axis.title.y.right = element_text(
      color = "grey30"
    ),
    axis.text.y.right = element_text(
      color = "grey30"
    )
  ) +
  labs(
    x = NULL,
    title = NULL,
    color = "Sampling point",
    shape = "Matrix"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure2_G_TOC_water_sediment.tiff"
  ),
  plot = figure2_g,
  width = 12,
  height = 4,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure2_G_TOC_water_sediment.svg"
  ),
  plot = figure2_g,
  width = 12,
  height = 4,
  units = "in"
)


# ============================================================================ #
# Supplementary Figure 2 | Collinearity and correlation matrices
# ============================================================================ #

vars_phys <- c(
  "pH",
  "Salinity",
  "ORP",
  "Wc",
  "Temperature",
  "T_sup",
  "T_5cm"
)

vars_geochem <- c(
  "TOC",
  "Cl",
  "SO4",
  "NO3",
  "Al",
  "As",
  "S",
  "Ba",
  "B",
  "Cd",
  "Ca",
  "Cr",
  "Co",
  "Cu",
  "Fe",
  "K",
  "Mg",
  "Mn",
  "Ni",
  "P",
  "Mo",
  "Na",
  "Li",
  "Pb",
  "Ag",
  "Se",
  "Sr",
  "Si",
  "Te",
  "Tl",
  "Ti",
  "V",
  "Zn",
  "Rb"
)


# ---------------------------------------------------------------------------- #
# Helper | Collinearity screen
# ---------------------------------------------------------------------------- #

run_collinearity_screen <- function(
    data,
    vars,
    label,
    vif_threshold = 5,
    cor_threshold = 0.80,
    max_na_prop = 0.20
) {

  vars_present <- intersect(
    vars,
    names(
      data
    )
  )

  vars_missing <- setdiff(
    vars,
    names(
      data
    )
  )

  if (
    length(
      vars_present
    ) <
      2
  ) {
    stop(
      "Fewer than two requested variables were found for ",
      label,
      "."
    )
  }

  x <- data[
    ,
    vars_present,
    drop = FALSE
  ]

  x[] <- lapply(
    x,
    function(z) {

      z <- suppressWarnings(
        as.numeric(
          as.character(
            z
          )
        )
      )

      z[
        !is.finite(
          z
        )
      ] <- NA_real_

      z
    }
  )

  na_prop <- vapply(
    x,
    function(z) {
      mean(
        is.na(
          z
        )
      )
    },
    numeric(
      1
    )
  )

  removed_missing <- names(
    na_prop[
      na_prop >
        max_na_prop
    ]
  )

  x <- x[
    ,
    !names(
      x
    ) %in%
      removed_missing,
    drop = FALSE
  ]

  x <- x[
    complete.cases(
      x
    ),
    ,
    drop = FALSE
  ]

  if (
    nrow(
      x
    ) <
      5
  ) {
    stop(
      "Fewer than five complete samples remain for ",
      label,
      "."
    )
  }

  nzv_idx <- caret::nearZeroVar(
    x
  )

  if (
    length(
      nzv_idx
    ) >
      0
  ) {
    removed_nzv <- names(
      x
    )[
      nzv_idx
    ]

    x <- x[
      ,
      -nzv_idx,
      drop = FALSE
    ]
  } else {
    removed_nzv <- character(
      0
    )
  }

  cor_matrix <- cor(
    x,
    method = "pearson"
  )

  idx_high_cor <- which(
    abs(
      cor_matrix
    ) >=
      cor_threshold &
      upper.tri(
        cor_matrix
      ),
    arr.ind = TRUE
  )

  if (
    nrow(
      idx_high_cor
    ) >
      0
  ) {

    high_cor_pairs <- data.frame(
      Variable_1 =
        rownames(
          cor_matrix
        )[
          idx_high_cor[
            ,
            1
          ]
        ],
      Variable_2 =
        colnames(
          cor_matrix
        )[
          idx_high_cor[
            ,
            2
          ]
        ],
      Pearson_r =
        cor_matrix[
          idx_high_cor
        ],
      row.names = NULL
    ) |>
      arrange(
        desc(
          abs(
            Pearson_r
          )
        )
      )

  } else {

    high_cor_pairs <- data.frame()
  }

  linear_combos <- caret::findLinearCombos(
    x
  )

  if (
    !is.null(
      linear_combos$remove
    ) &&
      length(
        linear_combos$remove
      ) >
        0
  ) {

    removed_linear <- names(
      x
    )[
      linear_combos$remove
    ]

    x <- x[
      ,
      -linear_combos$remove,
      drop = FALSE
    ]

  } else {

    removed_linear <- character(
      0
    )
  }

  removed_rank_guard <- character(
    0
  )

  while (
    ncol(
      x
    ) >=
      (
        nrow(
          x
        ) -
          1
      ) &&
      ncol(
        x
      ) >
        2
  ) {

    cor_tmp <- abs(
      cor(
        x,
        method = "pearson"
      )
    )

    diag(
      cor_tmp
    ) <- NA

    mean_abs_cor <- colMeans(
      cor_tmp,
      na.rm = TRUE
    )

    drop_var <- names(
      which.max(
        mean_abs_cor
      )
    )

    removed_rank_guard <- c(
      removed_rank_guard,
      drop_var
    )

    x <- x[
      ,
      !names(
        x
      ) %in%
        drop_var,
      drop = FALSE
    ]
  }

  if (
    ncol(
      x
    ) >=
      2
  ) {

    vif_object <- usdm::vifstep(
      x,
      th = vif_threshold,
      method = "pearson"
    )

    vif_table <- as.data.frame(
      vif_object@results
    )

    if (
      "Variables" %in%
        names(
          vif_table
        )
    ) {
      retained_vars <- as.character(
        vif_table$Variables
      )
    } else {
      retained_vars <- as.character(
        vif_table[
          [
            1
          ]
        ]
      )
    }

    removed_vif <- setdiff(
      names(
        x
      ),
      retained_vars
    )

  } else {

    retained_vars <- names(
      x
    )

    removed_vif <- character(
      0
    )

    vif_table <- data.frame(
      Variables = retained_vars,
      VIF = 1
    )
  }

  list(
    label = label,
    n_samples = nrow(
      x
    ),
    requested_vars = vars,
    missing_vars = vars_missing,
    removed_missing = removed_missing,
    removed_nzv = removed_nzv,
    high_cor_pairs = high_cor_pairs,
    cor_matrix = cor_matrix,
    removed_linear = removed_linear,
    removed_rank_guard = removed_rank_guard,
    removed_vif = removed_vif,
    retained = retained_vars,
    vif_table = vif_table
  )
}


# ---------------------------------------------------------------------------- #
# Run physicochemical and phase-specific geochemical screens
# ---------------------------------------------------------------------------- #

col_phys <- run_collinearity_screen(
  data = metadata,
  vars = vars_phys,
  label = "Physicochemical_All",
  vif_threshold = 5,
  cor_threshold = 0.80,
  max_na_prop = 0.20
)

meta_wet_geochem <- metadata |>
  filter(
    Phase ==
      "Wet",
    !is.na(
      Al
    )
  )

meta_dry_geochem <- metadata |>
  filter(
    Phase ==
      "Dry",
    !is.na(
      Al
    )
  )

col_geochem_wet <- run_collinearity_screen(
  data = meta_wet_geochem,
  vars = vars_geochem,
  label = "Geochemistry_Wet",
  vif_threshold = 5,
  cor_threshold = 0.80,
  max_na_prop = 0.20
)

col_geochem_dry <- run_collinearity_screen(
  data = meta_dry_geochem,
  vars = vars_geochem,
  label = "Geochemistry_Dry",
  vif_threshold = 5,
  cor_threshold = 0.80,
  max_na_prop = 0.20
)


# ---------------------------------------------------------------------------- #
# Export VIF/collinearity diagnostics
# ---------------------------------------------------------------------------- #

openxlsx::write.xlsx(
  list(
    Physicochemical_VIF =
      col_phys$vif_table,
    Physicochemical_high_cor =
      col_phys$high_cor_pairs,
    Wet_geochemistry_VIF =
      col_geochem_wet$vif_table,
    Wet_geochemistry_high_cor =
      col_geochem_wet$high_cor_pairs,
    Dry_geochemistry_VIF =
      col_geochem_dry$vif_table,
    Dry_geochemistry_high_cor =
      col_geochem_dry$high_cor_pairs
  ),
  file.path(
    table_dir,
    "environmental_collinearity_VIF.xlsx"
  ),
  overwrite = TRUE
)


# ---------------------------------------------------------------------------- #
# Helper | Publication-style correlation matrix
# ---------------------------------------------------------------------------- #

make_corr_heatmap <- function(
    cor_mat,
    var_order,
    title = NULL,
    axis_size = 8,
    show_legend = FALSE,
    show_values = FALSE,
    value_size = 3
) {

  var_order <- intersect(
    var_order,
    rownames(
      cor_mat
    )
  )

  cor_mat <- cor_mat[
    var_order,
    var_order,
    drop = FALSE
  ]

  corr_df <- expand.grid(
    row_var = var_order,
    col_var = var_order,
    stringsAsFactors = FALSE
  )

  corr_df$row_id <- match(
    corr_df$row_var,
    var_order
  )

  corr_df$col_id <- match(
    corr_df$col_var,
    var_order
  )

  corr_df$r <- mapply(
    function(
        i,
        j
    ) {
      cor_mat[
        i,
        j
      ]
    },
    corr_df$row_id,
    corr_df$col_id
  )

  corr_df <- corr_df |>
    filter(
      row_id <
        col_id
    ) |>
    mutate(
      abs_r = abs(
        r
      ),
      square_side =
        0.08 +
        0.84 *
        sqrt(
          abs_r
        ),
      text_colour =
        ifelse(
          abs_r >=
            0.50,
          "white",
          "black"
        ),
      row_var = factor(
        row_var,
        levels = rev(
          var_order
        )
      ),
      col_var = factor(
        col_var,
        levels = var_order
      )
    )

  p <- ggplot(
    corr_df,
    aes(
      x = col_var,
      y = row_var
    )
  ) +
    geom_tile(
      width = 0.96,
      height = 0.96,
      fill = NA,
      color = "grey80",
      linewidth = 0.30
    ) +
    geom_tile(
      aes(
        width = square_side,
        height = square_side,
        fill = r
      ),
      color = "grey65",
      linewidth = 0.15
    ) +
    geom_point(
      aes(
        size = abs_r
      ),
      shape = 22,
      alpha = 0,
      fill = "grey60",
      color = "grey40",
      show.legend = show_legend
    ) +
    scale_fill_gradient2(
      low = "#A6611A",
      mid = "white",
      high = "#018571",
      midpoint = 0,
      limits = c(
        -1,
        1
      ),
      breaks = c(
        -1,
        -0.5,
        0,
        0.5,
        1
      ),
      name = "Pearson r"
    ) +
    scale_size_continuous(
      limits = c(
        0,
        1
      ),
      breaks = c(
        0.25,
        0.50,
        0.75,
        1.00
      ),
      range = c(
        2,
        7
      ),
      name = "|Pearson r|"
    ) +
    scale_x_discrete(
      position = "top",
      expand = c(
        0,
        0
      )
    ) +
    scale_y_discrete(
      position = "right",
      expand = c(
        0,
        0
      )
    ) +
    coord_fixed(
      clip = "off"
    ) +
    labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    guides(
      fill = guide_colorbar(
        order = 1,
        title.position = "top",
        title.hjust = 0.5,
        barheight = grid::unit(
          4.5,
          "cm"
        ),
        barwidth = grid::unit(
          0.45,
          "cm"
        )
      ),
      size = guide_legend(
        order = 2,
        title.position = "top",
        override.aes = list(
          alpha = 1,
          shape = 22,
          fill = "grey55",
          color = "grey40"
        )
      )
    ) +
    theme_classic(
      base_size = 10
    ) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 0,
        vjust = 0.5,
        size = axis_size,
        color = "black"
      ),
      axis.text.y = element_text(
        hjust = 0,
        size = axis_size,
        color = "black",
        margin = margin(
          l = 4
        )
      ),
      plot.title = element_text(
        size = 11,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 8
        )
      ),
      plot.margin = margin(
        t = 8,
        r = 8,
        b = 5,
        l = 5
      ),
      legend.position =
        if (
          show_legend
        ) {
          "right"
        } else {
          "none"
        },
      legend.box = "vertical",
      legend.title = element_text(
        size = 9.5,
        face = "bold"
      ),
      legend.text = element_text(
        size = 8.5
      )
    )

  if (
    show_values
  ) {

    p <- p +
      geom_text(
        aes(
          label = sprintf(
            "%.2f",
            r
          ),
          color = text_colour
        ),
        size = value_size,
        fontface = "bold",
        show.legend = FALSE
      ) +
      scale_color_identity()
  }

  p
}


# ---------------------------------------------------------------------------- #
# Shared Wet/Dry geochemical ordering
# ---------------------------------------------------------------------------- #

geo_common <- intersect(
  rownames(
    col_geochem_wet$cor_matrix
  ),
  rownames(
    col_geochem_dry$cor_matrix
  )
)

cor_wet_common <- col_geochem_wet$cor_matrix[
  geo_common,
  geo_common,
  drop = FALSE
]

cor_dry_common <- col_geochem_dry$cor_matrix[
  geo_common,
  geo_common,
  drop = FALSE
]

cor_reference <- (
  abs(
    cor_wet_common
  ) +
    abs(
      cor_dry_common
    )
) /
  2

diag(
  cor_reference
) <- 1

geo_hclust <- hclust(
  as.dist(
    1 -
      cor_reference
  ),
  method = "complete"
)

geo_order <- geo_common[
  geo_hclust$order
]

phys_order <- c(
  "Salinity",
  "Wc",
  "pH",
  "ORP"
)


# ---------------------------------------------------------------------------- #
# Build Supplementary Figure 2
# ---------------------------------------------------------------------------- #

p_corr_wet <- make_corr_heatmap(
  cor_mat = cor_wet_common,
  var_order = geo_order,
  title = "Wet geochemistry",
  axis_size = 7.3,
  show_legend = FALSE,
  show_values = FALSE
)

p_corr_dry <- make_corr_heatmap(
  cor_mat = cor_dry_common,
  var_order = geo_order,
  title = "Dry geochemistry",
  axis_size = 7.3,
  show_legend = FALSE,
  show_values = FALSE
)

p_corr_phys <- make_corr_heatmap(
  cor_mat = col_phys$cor_matrix,
  var_order = phys_order,
  title = "Physicochemical variables",
  axis_size = 9.5,
  show_legend = FALSE,
  show_values = TRUE,
  value_size = 3.2
)

p_legend <- make_corr_heatmap(
  cor_mat = cor_wet_common,
  var_order = geo_order,
  title = NULL,
  axis_size = 7,
  show_legend = TRUE,
  show_values = FALSE
) +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )

shared_legend <- cowplot::get_legend(
  p_legend
)

p_corr_wet_labeled <- plot_grid(
  p_corr_wet,
  labels = "A",
  label_fontface = "bold",
  label_size = 13,
  label_x = 0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1
)

p_corr_dry_labeled <- plot_grid(
  p_corr_dry,
  labels = "B",
  label_fontface = "bold",
  label_size = 13,
  label_x = 0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1
)

p_corr_phys_labeled <- plot_grid(
  p_corr_phys,
  labels = "C",
  label_fontface = "bold",
  label_size = 13,
  label_x = 0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1
)

top_row <- plot_grid(
  p_corr_wet_labeled,
  p_corr_dry_labeled,
  ncol = 2,
  rel_widths = c(
    1,
    1
  ),
  align = "h"
)

bottom_row <- plot_grid(
  NULL,
  p_corr_phys_labeled,
  NULL,
  ncol = 3,
  rel_widths = c(
    0.34,
    0.32,
    0.34
  )
)

combined_corr <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(
    1,
    0.40
  )
)

supplementary_figure2 <- plot_grid(
  combined_corr,
  shared_legend,
  ncol = 2,
  rel_widths = c(
    1,
    0.075
  )
)

ggsave(
  file.path(
    figure_dir,
    "Supplementary_Figure2_environmental_correlations.tiff"
  ),
  plot = supplementary_figure2,
  width = 15,
  height = 9.5,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Supplementary_Figure2_environmental_correlations.svg"
  ),
  plot = supplementary_figure2,
  width = 15,
  height = 9.5,
  units = "in"
)


# ============================================================================ #
# Export correlation matrices
# ============================================================================ #

openxlsx::write.xlsx(
  list(
    Wet_geochemistry =
      as.data.frame(
        cor_wet_common
      ),
    Dry_geochemistry =
      as.data.frame(
        cor_dry_common
      ),
    Physicochemical =
      as.data.frame(
        col_phys$cor_matrix
      )
  ),
  file.path(
    table_dir,
    "Supplementary_Figure2_correlation_matrices.xlsx"
  ),
  rowNames = TRUE,
  overwrite = TRUE
)


# ============================================================================ #
# Final summary
# ============================================================================ #

cat(
  "\nEnvironmental analysis completed.\n"
)

cat(
  "\nFigure 2 components exported:\n",
  "  A-D | Wc, salinity, ORP, pH trajectories\n",
  "  E   | Wet basin comparisons\n",
  "  F   | Dry basin comparisons\n",
  "  G   | TOC water/sediment\n",
  sep = ""
)

cat(
  "\nSupplementary Figure 2 exported.\n"
)

cat(
  "\nPhysicochemical variables retained after VIF screening:\n"
)

print(
  col_phys$retained
)

cat(
  "\nWet geochemical variables retained after VIF screening:\n"
)

print(
  col_geochem_wet$retained
)

cat(
  "\nDry geochemical variables retained after VIF screening:\n"
)

print(
  col_geochem_dry$retained
)
