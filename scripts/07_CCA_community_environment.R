# ============================================================================ #
# 07_CCA_community_environment.R
# ============================================================================ #
#
# Purpose:
# Reproduce the constrained-ordination analyses used in Figure 6A-C:
#
#   - Detrended correspondence analysis (DCA) used to justify CCA
#   - Figure 6A: full-dataset CCA using pH, water content, ORP, and salinity
#   - Figure 6B: exploratory dry-phase geochemical CCA
#   - Figure 6C: exploratory wet-phase geochemical CCA
#   - envfit overlays for the ten most abundant families
#
# Notes:
# - Community analyses use non-rarefied data.
# - Relative abundance is used for constrained ordination.
# - Rare ASVs are retained if they reach >=0.1% relative abundance in >=2%
#   of samples OR >=1% relative abundance in any sample.
# - In the full CCA, Wc and Salinity are log1p-transformed; pH and ORP are
#   retained on their original scales.
# - Phase-specific CCAs use the same filtered community table and complete
#   geochemical subsets. All selected geochemical variables are log1p-
#   transformed, matching the original workflow.
# - Global models, constrained axes, sequential terms, and marginal effects
#   are evaluated with 999 permutations.
# - Family envfit P values are corrected with Benjamini-Hochberg across the
#   ten displayed families within each ordination.
#
# Input:
# - output_phyloseq_tree/ps_mdecon_v2_FTree.rds
#
# Main outputs:
# - Figure6A_CCA_full.*
# - Figure6B_CCA_dry.*
# - Figure6C_CCA_wet.*
# - CCA permutation-test tables
# - family envfit tables
# - DCA gradient-length summary
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
library(ggvegan)
library(ggforce)
library(microViz)
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

output_dir <- "output_CCA"
figure_dir <- file.path(
  output_dir,
  "figures"
)
table_dir <- file.path(
  output_dir,
  "tables"
)
text_dir <- file.path(
  output_dir,
  "summaries"
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

dir.create(
  text_dir,
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

shape_values <- c(
  22,
  21,
  24
)


# ============================================================================ #
# Helper functions
# ============================================================================ #

psotu2veg <- function(
    physeq
) {

  OTU <- otu_table(
    physeq
  )

  if (
    taxa_are_rows(
      OTU
    )
  ) {
    OTU <- t(
      OTU
    )
  }

  as(
    OTU,
    "matrix"
  )
}


pssd2veg <- function(
    physeq
) {

  as(
    sample_data(
      physeq
    ),
    "data.frame"
  )
}


anova_to_df <- function(
    x,
    test_name
) {

  out <- as.data.frame(
    x
  )

  out <- tibble::rownames_to_column(
    out,
    "Term"
  )

  out$Test <- test_name

  out
}


get_constrained_axis_variance <- function(
    model
) {

  prop_var <- summary(
    model
  )$concont$importance[
    "Proportion Explained",
    ,
    drop = TRUE
  ]

  c(
    CCA1 = round(
      prop_var[1] * 100,
      1
    ),
    CCA2 = round(
      prop_var[2] * 100,
      1
    )
  )
}


get_constrained_inertia_percent <- function(
    model
) {

  100 *
    model$CCA$tot.chi /
    model$tot.chi
}


make_family_envfit <- function(
    ps_object,
    cca_model,
    env_rows,
    n_families = 10,
    permutations = 999
) {

  ps_family <- tax_glom(
    ps_object,
    taxrank = "Family"
  )

  otu_family <- psotu2veg(
    ps_family
  )

  tax_df <- as.data.frame(
    tax_table(
      ps_family
    )
  )

  family_names <- as.character(
    tax_df$Family
  )

  family_names[
    is.na(
      family_names
    ) |
      family_names == ""
  ] <- "Unclassified family"

  colnames(
    otu_family
  ) <- make.unique(
    family_names
  )

  otu_family <- otu_family[
    env_rows,
    ,
    drop = FALSE
  ]

  top_taxa <- sort(
    colSums(
      otu_family
    ),
    decreasing = TRUE
  )

  top_names <- names(
    top_taxa
  )[
    seq_len(
      min(
        n_families,
        length(
          top_taxa
        )
      )
    )
  ]

  otu_top <- otu_family[
    ,
    top_names,
    drop = FALSE
  ]

  fit_taxa <- envfit(
    cca_model,
    otu_top,
    permutations = permutations
  )

  scores_taxa <- vegan::scores(
    fit_taxa,
    display = "vectors"
  )

  pvals <- fit_taxa$vectors$pvals
  padj <- p.adjust(
    pvals,
    method = "BH"
  )

  out <- data.frame(
    Family = names(
      pvals
    ),
    CCA1 = scores_taxa[
      names(
        pvals
      ),
      1
    ],
    CCA2 = scores_taxa[
      names(
        pvals
      ),
      2
    ],
    r2 = fit_taxa$vectors$r[
      names(
        pvals
      )
    ],
    P_value = pvals,
    P_adjust_BH = padj,
    Significant_BH = ifelse(
      padj <= 0.05,
      "Significant",
      "Not significant"
    ),
    row.names = NULL
  )

  list(
    fit = fit_taxa,
    table = out
  )
}


prepare_plot_scores <- function(
    model,
    metadata,
    family_table
) {

  sites <- ggvegan::fortify(
    model,
    layers = "wa",
    scaling = 3
  )

  bp <- ggvegan::fortify(
    model,
    layers = "bp",
    scaling = 3
  )

  if (
    "cca1" %in%
      colnames(
        sites
      )
  ) {
    sites <- sites |>
      dplyr::rename(
        CCA1 = cca1,
        CCA2 = cca2
      )
  }

  if (
    "cca1" %in%
      colnames(
        bp
      )
  ) {
    bp <- bp |>
      dplyr::rename(
        CCA1 = cca1,
        CCA2 = cca2
      )
  }

  sites <- sites |>
    dplyr::mutate(
      Phase = metadata$Phase,
      Time_Interval = metadata$Time_Interval,
      Ephemeral_Basin = metadata$Ephemeral_Basin
    )

  taxa_df <- family_table |>
    dplyr::transmute(
      label = Family,
      CCA1 = CCA1 * 1.5,
      CCA2 = CCA2 * 1.5,
      pval = P_value,
      p_adj_BH = P_adjust_BH,
      signif = Significant_BH
    )

  list(
    sites = sites,
    bp = bp,
    taxa = taxa_df
  )
}


run_cca_tests <- function(
    model,
    permutations = 999
) {

  list(
    global = anova(
      model,
      permutations = permutations
    ),
    axis = anova(
      model,
      by = "axis",
      permutations = permutations
    ),
    terms = anova(
      model,
      by = "terms",
      permutations = permutations
    ),
    margin = anova(
      model,
      by = "margin",
      permutations = permutations
    )
  )
}


# ============================================================================ #
# Load and prepare final phyloseq object
# ============================================================================ #

ps <- readRDS(
  phyloseq_path
)

sample_data(
  ps
)$Ephemeral_Basin <- factor(
  sample_data(
    ps
  )$Ephemeral_Basin,
  levels = basin_levels
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
      x /
        sum(
          x,
          na.rm = TRUE
        )
    }
  }
)

# These validation steps were present in the original CCA workflow.
ps_ra <- microViz::tax_fix(
  ps_ra
)

ps_ra <- microViz::phyloseq_validate(
  ps_ra,
  remove_undetected = TRUE,
  verbose = TRUE
)


# ============================================================================ #
# Rare-ASV filter used for constrained ordination
# ============================================================================ #

min_prev_samples <- ceiling(
  0.02 *
    nsamples(
      ps_ra
    )
)

ps_ra_filt <- filter_taxa(
  ps_ra,
  function(x) {
    sum(
      x >= 0.001
    ) >=
      min_prev_samples ||
      max(
        x
      ) >=
      0.01
  },
  prune = TRUE
)

cat(
  "ASVs before CCA filter: ",
  ntaxa(
    ps_ra
  ),
  "\n",
  sep = ""
)

cat(
  "ASVs after CCA filter: ",
  ntaxa(
    ps_ra_filt
  ),
  "\n",
  sep = ""
)


# ============================================================================ #
# DCA gradient-length diagnostic
# ============================================================================ #

# The original DCA diagnostic was run on the final non-rarefied community
# matrix and returned a first-axis length of approximately 6.04 SD units.

otu_dca <- psotu2veg(
  ps
)

meta_dca <- pssd2veg(
  ps
)

dca_env_raw <- meta_dca |>
  dplyr::select(
    pH,
    Salinity,
    ORP,
    Wc,
    CE,
    RES
  )

dca_env <- dca_env_raw

dca_env[
  c(
    "Salinity",
    "Wc",
    "CE",
    "RES"
  )
] <- lapply(
  dca_env[
    c(
      "Salinity",
      "Wc",
      "CE",
      "RES"
    )
  ],
  log1p
)

keep_dca <- complete.cases(
  dca_env
)

otu_dca <- otu_dca[
  keep_dca,
  ,
  drop = FALSE
]

decorana_res <- vegan::decorana(
  otu_dca
)

writeLines(
  capture.output(
    decorana_res
  ),
  file.path(
    text_dir,
    "DCA_gradient_length.txt"
  )
)


# ============================================================================ #
# Figure 6A | Full-dataset CCA
# ============================================================================ #

otu_full <- psotu2veg(
  ps_ra_filt
)

meta_full <- pssd2veg(
  ps_ra_filt
)

meta_full$pH <- as.numeric(
  as.character(
    meta_full$pH
  )
)

meta_full$Wc <- log1p(
  as.numeric(
    as.character(
      meta_full$Wc
    )
  )
)

meta_full$ORP <- as.numeric(
  as.character(
    meta_full$ORP
  )
)

meta_full$Salinity <- log1p(
  as.numeric(
    as.character(
      meta_full$Salinity
    )
  )
)

env_full <- meta_full |>
  dplyr::select(
    pH,
    Wc,
    ORP,
    Salinity
  )

keep_full <- complete.cases(
  env_full
)

env_full <- env_full[
  keep_full,
  ,
  drop = FALSE
]

otu_full <- otu_full[
  rownames(
    env_full
  ),
  ,
  drop = FALSE
]

meta_full <- meta_full[
  rownames(
    env_full
  ),
  ,
  drop = FALSE
]

cca_full <- vegan::cca(
  otu_full ~
    pH +
    Wc +
    ORP +
    Salinity,
  data = env_full
)

vif_full <- vegan::vif.cca(
  cca_full
)

fit_env_full <- vegan::envfit(
  cca_full,
  env_full,
  permutations = 999
)

tests_full <- run_cca_tests(
  cca_full,
  permutations = 999
)

family_fit_full <- make_family_envfit(
  ps_object = ps_ra_filt,
  cca_model = cca_full,
  env_rows = rownames(
    env_full
  ),
  n_families = 10,
  permutations = 999
)

axis_var_full <- get_constrained_axis_variance(
  cca_full
)

inertia_full <- get_constrained_inertia_percent(
  cca_full
)

scores_full <- prepare_plot_scores(
  model = cca_full,
  metadata = meta_full,
  family_table = family_fit_full$table
)

figure6a <- ggplot() +
  ggforce::geom_mark_hull(
    data = scores_full$sites,
    aes(
      x = CCA1,
      y = CCA2,
      group = Phase,
      label = Phase
    ),
    concavity = 3,
    expand = grid::unit(
      1,
      "mm"
    ),
    alpha = 0.15,
    fill = "gray90",
    color = "gray40",
    linewidth = 0.5,
    label.fontsize = 10,
    label.lineheight = 0.75,
    label.fontface = c(
      "bold",
      "plain"
    ),
    label.fill = NA,
    label.colour = "gray40",
    con.colour = "gray40",
    con.size = 0.5,
    con.type = "straight"
  ) +
  geom_point(
    data = scores_full$sites,
    aes(
      x = CCA1,
      y = CCA2,
      fill = Time_Interval,
      shape = Ephemeral_Basin
    ),
    size = 4,
    alpha = 0.9,
    color = "black",
    stroke = 0.8
  ) +
  scale_shape_manual(
    values = shape_values
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    name = "Time interval (days)",
    trans = "sqrt",
    direction = 1,
    begin = 0.1,
    end = 1,
    breaks = c(
      0,
      28,
      66,
      100,
      200,
      300,
      455
    )
  ) +
  geom_segment(
    data = scores_full$bp,
    aes(
      x = 0,
      y = 0,
      xend = CCA1 * 2,
      yend = CCA2 * 2
    ),
    arrow = arrow(
      length = grid::unit(
        0.3,
        "cm"
      )
    ),
    color = "black",
    linewidth = 1
  ) +
  geom_text(
    data = scores_full$bp,
    aes(
      x = CCA1 * 2,
      y = CCA2 * 2,
      label = label
    ),
    color = "black",
    fontface = "bold",
    size = 5,
    vjust = -1
  ) +
  geom_point(
    data = scores_full$taxa,
    aes(
      x = CCA1,
      y = CCA2,
      color = signif
    ),
    size = 3
  ) +
  geom_text(
    data = scores_full$taxa,
    aes(
      x = CCA1,
      y = CCA2,
      label = label,
      color = signif
    ),
    fontface = "italic",
    size = 4,
    vjust = -0.7,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name = "Top 10 families",
    values = c(
      "Significant" = "darkred",
      "Not significant" = "gray55"
    )
  ) +
  theme_cowplot() +
  cowplot::background_grid(
    major = "xy",
    minor = "none"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(
      fill = NA,
      colour = "black"
    )
  ) +
  labs(
    title = NULL,
    x = paste0(
      "CCA1 (",
      axis_var_full[
        "CCA1"
      ],
      "%)"
    ),
    y = paste0(
      "CCA2 (",
      axis_var_full[
        "CCA2"
      ],
      "%)"
    )
  )

ggsave(
  file.path(
    figure_dir,
    "Figure6A_CCA_full.tiff"
  ),
  plot = figure6a,
  width = 11,
  height = 7,
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure6A_CCA_full.svg"
  ),
  plot = figure6a,
  width = 11,
  height = 7,
  device = "svg"
)


# ============================================================================ #
# Export full-dataset CCA statistics
# ============================================================================ #

full_summary <- data.frame(
  Analysis = "Full CCA",
  N_samples = nrow(
    env_full
  ),
  N_ASVs = ncol(
    otu_full
  ),
  Constrained_inertia_percent = inertia_full,
  CCA1_percent_of_constrained = axis_var_full[
    "CCA1"
  ],
  CCA2_percent_of_constrained = axis_var_full[
    "CCA2"
  ]
)

openxlsx::write.xlsx(
  list(
    Summary = full_summary,
    VIF = data.frame(
      Variable = names(
        vif_full
      ),
      VIF = as.numeric(
        vif_full
      )
    ),
    ANOVA_global = anova_to_df(
      tests_full$global,
      "global"
    ),
    ANOVA_axis = anova_to_df(
      tests_full$axis,
      "axis"
    ),
    ANOVA_terms = anova_to_df(
      tests_full$terms,
      "terms"
    ),
    ANOVA_margin = anova_to_df(
      tests_full$margin,
      "margin"
    ),
    Family_envfit = family_fit_full$table
  ),
  file = file.path(
    table_dir,
    "Figure6A_CCA_full_statistics.xlsx"
  )
)

writeLines(
  c(
    "CCA MODEL",
    capture.output(
      cca_full
    ),
    "",
    "CCA SUMMARY",
    capture.output(
      summary(
        cca_full
      )
    ),
    "",
    "VIF",
    capture.output(
      vif_full
    ),
    "",
    "ENVFIT ENVIRONMENT",
    capture.output(
      fit_env_full
    )
  ),
  file.path(
    text_dir,
    "Figure6A_CCA_full_summary.txt"
  )
)


# ============================================================================ #
# Phase-specific CCA helper
# ============================================================================ #

run_phase_cca <- function(
    ps_ra_filtered,
    phase_name,
    variables,
    viridis_option,
    time_breaks,
    figure_name
) {

  ps_phase <- subset_samples(
    ps_ra_filtered,
    Phase == phase_name
  )

  # Original workflow used Al availability as the indicator of samples with
  # geochemical measurements.
  ps_phase <- subset_samples(
    ps_phase,
    !is.na(
      Al
    )
  )

  meta_phase <- data.frame(
    sample_data(
      ps_phase
    ),
    check.names = FALSE
  )

  # Hg and Rb were removed from the phase-specific metadata in the original
  # workflow before model construction.
  meta_phase <- meta_phase[
    ,
    !names(
      meta_phase
    ) %in%
      c(
        "Hg",
        "Rb"
      ),
    drop = FALSE
  ]

  for (
    variable in variables
  ) {
    meta_phase[
      [
        variable
      ]
    ] <- log1p(
      as.numeric(
        as.character(
          meta_phase[
            [
              variable
            ]
          ]
        )
      )
    )
  }

  keep_samples <- complete.cases(
    meta_phase[
      ,
      variables,
      drop = FALSE
    ]
  )

  keep_ids <- rownames(
    meta_phase
  )[
    keep_samples
  ]

  ps_phase <- prune_samples(
    keep_ids,
    ps_phase
  )

  meta_phase <- data.frame(
    sample_data(
      ps_phase
    ),
    check.names = FALSE
  )

  for (
    variable in variables
  ) {
    meta_phase[
      [
        variable
      ]
    ] <- log1p(
      as.numeric(
        as.character(
          meta_phase[
            [
              variable
            ]
          ]
        )
      )
    )
  }

  sample_data(
    ps_phase
  ) <- sample_data(
    meta_phase
  )

  env_phase <- meta_phase[
    ,
    variables,
    drop = FALSE
  ]

  otu_phase <- psotu2veg(
    ps_phase
  )

  env_phase <- env_phase[
    complete.cases(
      env_phase
    ),
    ,
    drop = FALSE
  ]

  otu_phase <- otu_phase[
    rownames(
      env_phase
    ),
    ,
    drop = FALSE
  ]

  meta_phase <- meta_phase[
    rownames(
      env_phase
    ),
    ,
    drop = FALSE
  ]

  cca_phase <- vegan::cca(
    otu_phase ~ .,
    data = env_phase
  )

  tests_phase <- run_cca_tests(
    cca_phase,
    permutations = 999
  )

  fit_env_phase <- vegan::envfit(
    cca_phase,
    env_phase,
    permutations = 999
  )

  family_fit_phase <- make_family_envfit(
    ps_object = ps_phase,
    cca_model = cca_phase,
    env_rows = rownames(
      env_phase
    ),
    n_families = 10,
    permutations = 999
  )

  axis_var_phase <- get_constrained_axis_variance(
    cca_phase
  )

  inertia_phase <- get_constrained_inertia_percent(
    cca_phase
  )

  scores_phase <- prepare_plot_scores(
    model = cca_phase,
    metadata = meta_phase,
    family_table = family_fit_phase$table
  )

  p <- ggplot() +
    geom_point(
      data = scores_phase$sites,
      aes(
        x = CCA1,
        y = CCA2,
        fill = Time_Interval,
        shape = Ephemeral_Basin
      ),
      size = 4,
      alpha = 0.9,
      color = "black",
      stroke = 0.8
    ) +
    scale_shape_manual(
      values = shape_values
    ) +
    scale_fill_viridis_c(
      option = viridis_option,
      name = "Time interval (days)",
      trans = "sqrt",
      direction = 1,
      begin = 0.1,
      end = 1,
      breaks = time_breaks
    ) +
    geom_segment(
      data = scores_phase$bp,
      aes(
        x = 0,
        y = 0,
        xend = CCA1 * 2,
        yend = CCA2 * 2
      ),
      arrow = arrow(
        length = grid::unit(
          0.3,
          "cm"
        )
      ),
      color = "black",
      linewidth = 1
    ) +
    geom_text(
      data = scores_phase$bp,
      aes(
        x = CCA1,
        y = CCA2,
        label = label
      ),
      color = "black",
      fontface = "bold",
      size = 5,
      vjust = -1
    ) +
    geom_point(
      data = scores_phase$taxa,
      aes(
        x = CCA1,
        y = CCA2,
        color = signif
      ),
      size = 2.2
    ) +
    geom_text(
      data = scores_phase$taxa,
      aes(
        x = CCA1,
        y = CCA2,
        label = label,
        color = signif
      ),
      fontface = "italic",
      size = 4,
      vjust = -0.7,
      show.legend = FALSE
    ) +
    scale_color_manual(
      name = "Top 10 families",
      values = c(
        "Significant" = "darkred",
        "Not significant" = "gray55"
      )
    ) +
    theme_cowplot() +
    cowplot::background_grid(
      major = "xy",
      minor = "none"
    ) +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(
        fill = NA,
        colour = "black"
      )
    ) +
    labs(
      title = NULL,
      x = paste0(
        "CCA1 (",
        axis_var_phase[
          "CCA1"
        ],
        "%)"
      ),
      y = paste0(
        "CCA2 (",
        axis_var_phase[
          "CCA2"
        ],
        "%)"
      )
    )

  ggsave(
    file.path(
      figure_dir,
      paste0(
        figure_name,
        ".tiff"
      )
    ),
    plot = p,
    width = 10,
    height = 6,
    device = "tiff",
    dpi = 300,
    compression = "lzw"
  )

  ggsave(
    file.path(
      figure_dir,
      paste0(
        figure_name,
        ".svg"
      )
    ),
    plot = p,
    width = 10,
    height = 6,
    device = "svg"
  )

  summary_df <- data.frame(
    Phase = phase_name,
    N_samples = nrow(
      env_phase
    ),
    N_ASVs = ncol(
      otu_phase
    ),
    Predictors = paste(
      variables,
      collapse = ", "
    ),
    Constrained_inertia_percent = inertia_phase,
    CCA1_percent_of_constrained = axis_var_phase[
      "CCA1"
    ],
    CCA2_percent_of_constrained = axis_var_phase[
      "CCA2"
    ]
  )

  list(
    model = cca_phase,
    plot = p,
    summary = summary_df,
    tests = tests_phase,
    envfit = fit_env_phase,
    family_envfit = family_fit_phase$table
  )
}


# ============================================================================ #
# Figure 6B | Dry-phase geochemical CCA
# ============================================================================ #

dry_variables <- c(
  "TOC",
  "Na",
  "Mo",
  "Se",
  "Sr",
  "Pb"
)

dry_cca <- run_phase_cca(
  ps_ra_filtered = ps_ra_filt,
  phase_name = "Dry",
  variables = dry_variables,
  viridis_option = "A",
  time_breaks = c(
    66,
    100,
    200,
    300,
    455
  ),
  figure_name = "Figure6B_CCA_dry"
)

openxlsx::write.xlsx(
  list(
    Summary = dry_cca$summary,
    ANOVA_global = anova_to_df(
      dry_cca$tests$global,
      "global"
    ),
    ANOVA_axis = anova_to_df(
      dry_cca$tests$axis,
      "axis"
    ),
    ANOVA_terms = anova_to_df(
      dry_cca$tests$terms,
      "terms"
    ),
    ANOVA_margin = anova_to_df(
      dry_cca$tests$margin,
      "margin"
    ),
    Family_envfit = dry_cca$family_envfit
  ),
  file = file.path(
    table_dir,
    "Figure6B_CCA_dry_statistics.xlsx"
  )
)

writeLines(
  c(
    "DRY CCA MODEL",
    capture.output(
      dry_cca$model
    ),
    "",
    "DRY CCA SUMMARY",
    capture.output(
      summary(
        dry_cca$model
      )
    ),
    "",
    "DRY ENVFIT",
    capture.output(
      dry_cca$envfit
    )
  ),
  file.path(
    text_dir,
    "Figure6B_CCA_dry_summary.txt"
  )
)


# ============================================================================ #
# Figure 6C | Wet-phase geochemical CCA
# ============================================================================ #

wet_variables <- c(
  "TOC",
  "Cl",
  "NO3",
  "Fe",
  "Li",
  "Pb"
)

wet_cca <- run_phase_cca(
  ps_ra_filtered = ps_ra_filt,
  phase_name = "Wet",
  variables = wet_variables,
  viridis_option = "viridis",
  time_breaks = c(
    0,
    7,
    14,
    21,
    28
  ),
  figure_name = "Figure6C_CCA_wet"
)

openxlsx::write.xlsx(
  list(
    Summary = wet_cca$summary,
    ANOVA_global = anova_to_df(
      wet_cca$tests$global,
      "global"
    ),
    ANOVA_axis = anova_to_df(
      wet_cca$tests$axis,
      "axis"
    ),
    ANOVA_terms = anova_to_df(
      wet_cca$tests$terms,
      "terms"
    ),
    ANOVA_margin = anova_to_df(
      wet_cca$tests$margin,
      "margin"
    ),
    Family_envfit = wet_cca$family_envfit
  ),
  file = file.path(
    table_dir,
    "Figure6C_CCA_wet_statistics.xlsx"
  )
)

writeLines(
  c(
    "WET CCA MODEL",
    capture.output(
      wet_cca$model
    ),
    "",
    "WET CCA SUMMARY",
    capture.output(
      summary(
        wet_cca$model
      )
    ),
    "",
    "WET ENVFIT",
    capture.output(
      wet_cca$envfit
    )
  ),
  file.path(
    text_dir,
    "Figure6C_CCA_wet_summary.txt"
  )
)


# ============================================================================ #
# Combined CCA summary
# ============================================================================ #

cca_summary_all <- dplyr::bind_rows(
  full_summary,
  dry_cca$summary,
  wet_cca$summary
)

openxlsx::write.xlsx(
  cca_summary_all,
  file.path(
    table_dir,
    "Figure6A-C_CCA_summary.xlsx"
  )
)


# ============================================================================ #
# Final summary
# ============================================================================ #

cat(
  "\nCCA analyses completed.\n"
)

cat(
  "DCA first-axis length should reproduce approximately 6.04 SD units.\n"
)

cat(
  "Filtered ASVs: ",
  ntaxa(
    ps_ra_filt
  ),
  "\n",
  sep = ""
)

cat(
  "Full CCA samples: ",
  nrow(
    env_full
  ),
  "\n",
  sep = ""
)

cat(
  "Full CCA constrained inertia (%): ",
  round(
    inertia_full,
    2
  ),
  "\n",
  sep = ""
)

cat(
  "Dry CCA samples: ",
  dry_cca$summary$N_samples,
  "\n",
  sep = ""
)

cat(
  "Wet CCA samples: ",
  wet_cca$summary$N_samples,
  "\n",
  sep = ""
)
