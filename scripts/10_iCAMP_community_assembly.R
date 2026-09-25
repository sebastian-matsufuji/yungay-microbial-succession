# ============================================================================ #
# 10_iCAMP_community_assembly.R
# ============================================================================ #
#
# Purpose:
# Reproduce the final iCAMP community-assembly analysis supporting Figure 8
# and the iCAMP wet-dry statistics reported in Table Supplementary 9.
#
# This curated script keeps the final analysis only. Pilot runs and bin-size
# sensitivity tests from the historical master script are intentionally omitted.
#
# Final settings:
# - non-rarefied counts
# - analyses run independently by basin
# - phylogenetic-distance cutoff ds = 0.2
# - minimum bin size = 12 ASVs for HCP, 24 for CHF and YSB
# - 999 randomizations
# - abundance-weighted bMPD
# - Confidence significance index
# - Bray-Curtis taxonomic turnover
# - within-bin phylogenetic randomization
# - basin-wide taxonomic randomization
# - seeds: HCP 20260901, CHF 20260902, YSB 20260903
#
# Temporal summaries use consecutive observed transitions within the same
# permanent sampling point. Wet-Dry boundary-crossing transitions are excluded
# from phase summaries.
#
# Figure 8 was assembled from two R-generated components:
# - overall basin process composition (circular diagrams)
# - temporal process trajectories
#
# ============================================================================ #


# ============================================================================ #
# Packages
# ============================================================================ #

library(phyloseq)
library(ape)
library(iCAMP)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(openxlsx)


# ============================================================================ #
# Paths
# ============================================================================ #

phyloseq_path <- file.path(
  "output_phyloseq_tree",
  "ps_mdecon_v2_FTree.rds"
)

output_dir <- "output_iCAMP"
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")

for (d in c(output_dir, figure_dir, table_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}


# ============================================================================ #
# Final settings
# ============================================================================ #

process_cols <- c(
  "Heterogeneous.Selection",
  "Homogeneous.Selection",
  "Dispersal.Limitation",
  "Homogenizing.Dispersal",
  "Drift.and.Others"
)

final_settings <- data.frame(
  Basin = c("HCP", "CHF", "YSB"),
  BinSize = c(12, 24, 24),
  Seed = c(20260901, 20260902, 20260903),
  stringsAsFactors = FALSE
)

basin_labels <- c(
  "HCP" = "Herradura Clay Pan",
  "CHF" = "Ckoirama Halite Field",
  "YSB" = "Yungay Station Basin"
)

process_colors <- c(
  "Heterogeneous.Selection" = "#8c2d04",
  "Homogeneous.Selection" = "#f4a582",
  "Dispersal.Limitation" = "#2166ac",
  "Homogenizing.Dispersal" = "#92c5de",
  "Drift.and.Others" = "#4D4D4D"
)

process_labels <- c(
  "Heterogeneous.Selection" = "Heterogeneous selection",
  "Homogeneous.Selection" = "Homogeneous selection",
  "Dispersal.Limitation" = "Dispersal limitation",
  "Homogenizing.Dispersal" = "Homogenizing dispersal",
  "Drift.and.Others" = "Drift and Others"
)

rain_date <- as.Date("2017-06-07")


# ============================================================================ #
# Load final non-rarefied phyloseq object
# ============================================================================ #

ps <- readRDS(phyloseq_path)


# ============================================================================ #
# Helper | Prepare one basin for iCAMP
# ============================================================================ #

prepare_basin_icamp <- function(ps_object, basin_name) {

  basin_dir <- file.path(output_dir, basin_name)
  input_dir <- file.path(basin_dir, "input")
  final_dir <- file.path(basin_dir, "final")

  dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

  ps_basin <- subset_samples(
    ps_object,
    Basin == basin_name
  )

  ps_basin <- prune_taxa(
    taxa_sums(ps_basin) > 0,
    ps_basin
  )

  ps_basin <- prune_samples(
    sample_sums(ps_basin) > 0,
    ps_basin
  )

  comm <- as(
    otu_table(ps_basin),
    "matrix"
  )

  if (taxa_are_rows(ps_basin)) {
    comm <- t(comm)
  }

  storage.mode(comm) <- "numeric"

  tree <- phy_tree(ps_basin)

  common_taxa <- intersect(
    colnames(comm),
    tree$tip.label
  )

  comm <- comm[, common_taxa, drop = FALSE]
  tree <- ape::keep.tip(tree, common_taxa)

  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  tree <- ape::keep.tip(tree, colnames(comm))

  metadata <- data.frame(
    sample_data(ps_basin),
    check.names = FALSE
  )

  metadata <- metadata[
    rownames(comm),
    ,
    drop = FALSE
  ]

  if (!all(comm >= 0)) {
    stop("Negative community counts detected in ", basin_name)
  }

  if (anyNA(comm)) {
    stop("NA community counts detected in ", basin_name)
  }

  if (!setequal(colnames(comm), tree$tip.label)) {
    stop("Community/tree taxa mismatch in ", basin_name)
  }

  saveRDS(
    ps_basin,
    file.path(input_dir, paste0("ps_", basin_name, "_icamp.rds"))
  )

  saveRDS(
    comm,
    file.path(input_dir, paste0("comm_", basin_name, "_icamp.rds"))
  )

  saveRDS(
    tree,
    file.path(input_dir, paste0("tree_", basin_name, "_icamp.rds"))
  )

  saveRDS(
    metadata,
    file.path(input_dir, paste0("metadata_", basin_name, "_icamp.rds"))
  )

  pd_dir <- file.path(input_dir, "phylogenetic_distance")
  dir.create(pd_dir, showWarnings = FALSE, recursive = TRUE)

  pd_big <- iCAMP::pdist.big(
    tree = tree,
    wd = pd_dir,
    nworker = 4,
    treepath.file = "tree.path.rda",
    pd.spname.file = "pd.taxon.name.csv",
    pd.backingfile = "pd.bin",
    pd.desc.file = "pd.desc"
  )

  saveRDS(
    pd_big,
    file.path(input_dir, paste0("pd_big_", basin_name, ".rds"))
  )

  tree_root_result <- iCAMP::midpoint.root.big(
    tree = tree,
    pd.desc = pd_big$pd.file,
    pd.spname = pd_big$tip.label,
    pd.wd = pd_big$pd.wd,
    nworker = 4
  )

  tree_rooted <- tree_root_result$tree

  saveRDS(
    tree_rooted,
    file.path(
      input_dir,
      paste0("tree_", basin_name, "_icamp_midpoint_rooted.rds")
    )
  )

  ape::write.tree(
    tree_rooted,
    file = file.path(
      input_dir,
      paste0("tree_", basin_name, "_icamp_midpoint_rooted.nwk")
    )
  )

  list(
    Basin = basin_name,
    Samples = nrow(comm),
    ASVs = ncol(comm),
    Input_dir = input_dir,
    Final_dir = final_dir
  )
}


# ============================================================================ #
# Prepare HCP, CHF, and YSB
# ============================================================================ #

preparation_summary <- purrr::map_dfr(
  final_settings$Basin,
  ~ as.data.frame(
    prepare_basin_icamp(ps, .x),
    stringsAsFactors = FALSE
  )
)

openxlsx::write.xlsx(
  preparation_summary,
  file.path(table_dir, "iCAMP_input_summary.xlsx"),
  overwrite = TRUE
)


# ============================================================================ #
# Helper | Run final iCAMP analysis for one basin
# ============================================================================ #

run_icamp_final <- function(basin_name, bin_size, seed_value) {

  basin_dir <- file.path(output_dir, basin_name)
  input_dir <- file.path(basin_dir, "input")
  final_dir <- file.path(basin_dir, "final")

  comm <- readRDS(
    file.path(input_dir, paste0("comm_", basin_name, "_icamp.rds"))
  )

  tree <- readRDS(
    file.path(
      input_dir,
      paste0("tree_", basin_name, "_icamp_midpoint_rooted.rds")
    )
  )

  pd_big <- readRDS(
    file.path(input_dir, paste0("pd_big_", basin_name, ".rds"))
  )

  prefix <- paste0(
    "iCAMP_",
    basin_name,
    "_bMPD_Confidence_rand999_binsize",
    bin_size
  )

  complete_rds <- file.path(
    final_dir,
    paste0(prefix, "_complete_object.rds")
  )

  if (file.exists(complete_rds)) {

    icamp_result <- readRDS(complete_rds)
    run_status <- "REUSED"
    run_time <- NA

  } else {

    set.seed(seed_value)

    run_time <- system.time({
      icamp_result <- iCAMP::icamp.big(
        comm = comm,
        tree = tree,
        pd.desc = pd_big$pd.file,
        pd.spname = pd_big$tip.label,
        pd.wd = pd_big$pd.wd,
        rand = 999,
        prefix = prefix,
        output.wd = final_dir,
        ds = 0.2,
        bin.size.limit = bin_size,
        phylo.rand.scale = "within.bin",
        taxa.rand.scale = "across.all",
        phylo.metric = "bMPD",
        sig.index = "Confidence",
        taxo.metric = "bray",
        unit.sum = rowSums(comm),
        transform.method = NULL,
        dirichlet = FALSE,
        correct.special = TRUE,
        special.method = "depend",
        omit.option = "no",
        sp.check = TRUE,
        ignore.zero = TRUE,
        nworker = 4,
        memory.G = 16,
        rtree.save = FALSE,
        detail.save = TRUE,
        qp.save = TRUE,
        detail.null = FALSE
      )
    })

    saveRDS(icamp_result, complete_rds)

    saveRDS(
      list(
        Basin = basin_name,
        Randomizations = 999,
        BinSize = bin_size,
        ds = 0.2,
        phylo.metric = "bMPD",
        sig.index = "Confidence",
        phylo.rand.scale = "within.bin",
        taxa.rand.scale = "across.all",
        taxo.metric = "bray",
        Seed = seed_value,
        Runtime = run_time,
        iCAMP_version = as.character(packageVersion("iCAMP"))
      ),
      file.path(final_dir, paste0(prefix, "_run_info.rds"))
    )

    run_status <- "CALCULATED"
  }

  process_table <- icamp_result$CbMPDiCBraya

  write.table(
    process_table,
    file = file.path(final_dir, paste0(prefix, "_community_processes.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  process_sum <- rowSums(
    process_table[, process_cols, drop = FALSE]
  )

  data.frame(
    Basin = basin_name,
    BinSize = bin_size,
    Turnovers = nrow(process_table),
    Expected = choose(nrow(comm), 2),
    NA_values = sum(
      is.na(process_table[, process_cols, drop = FALSE])
    ),
    SumMin = min(process_sum),
    SumMax = max(process_sum),
    Status = run_status,
    stringsAsFactors = FALSE
  )
}


# ============================================================================ #
# Run final iCAMP analyses
# ============================================================================ #

final_summary <- purrr::pmap_dfr(
  final_settings,
  function(Basin, BinSize, Seed) {
    run_icamp_final(
      basin_name = Basin,
      bin_size = BinSize,
      seed_value = Seed
    )
  }
)

openxlsx::write.xlsx(
  final_summary,
  file.path(table_dir, "iCAMP_final_run_summary.xlsx"),
  overwrite = TRUE
)


# ============================================================================ #
# Temporal processing | consecutive observations within sampling point
# ============================================================================ #

make_pair_key <- function(a, b) {
  ifelse(
    a <= b,
    paste(a, b, sep = "||"),
    paste(b, a, sep = "||")
  )
}

extract_temporal_icamp <- function(basin_name, bin_size) {

  basin_dir <- file.path(output_dir, basin_name)
  input_dir <- file.path(basin_dir, "input")
  final_dir <- file.path(basin_dir, "final")

  metadata <- readRDS(
    file.path(input_dir, paste0("metadata_", basin_name, "_icamp.rds"))
  )

  prefix <- paste0(
    "iCAMP_",
    basin_name,
    "_bMPD_Confidence_rand999_binsize",
    bin_size
  )

  icamp_result <- readRDS(
    file.path(final_dir, paste0(prefix, "_complete_object.rds"))
  )

  process_table <- icamp_result$CbMPDiCBraya

  meta <- data.frame(
    SampleID = rownames(metadata),
    Sampling_Point = as.character(metadata$Sampling_Point),
    Phase = as.character(metadata$Phase),
    stringsAsFactors = FALSE
  )

  meta$TimeIndex <- as.integer(
    sub(".*_T([0-9]+)$", "\\1", meta$SampleID)
  )

  if (anyNA(meta$TimeIndex)) {
    stop("Could not extract temporal index in ", basin_name)
  }

  meta_split <- split(meta, meta$Sampling_Point)

  allowed_pairs <- lapply(
    meta_split,
    function(x) {

      x <- x[order(x$TimeIndex), , drop = FALSE]

      if (nrow(x) < 2) {
        return(NULL)
      }

      data.frame(
        Sampling_Point = x$Sampling_Point[-nrow(x)],
        sample_early = x$SampleID[-nrow(x)],
        sample_late = x$SampleID[-1],
        Time_early = x$TimeIndex[-nrow(x)],
        Time_late = x$TimeIndex[-1],
        Phase_early = x$Phase[-nrow(x)],
        Phase_late = x$Phase[-1],
        stringsAsFactors = FALSE
      )
    }
  ) |>
    dplyr::bind_rows()

  allowed_pairs$PairKey <- make_pair_key(
    allowed_pairs$sample_early,
    allowed_pairs$sample_late
  )

  process_table$PairKey <- make_pair_key(
    process_table$sample1,
    process_table$sample2
  )

  pair_match <- match(
    allowed_pairs$PairKey,
    process_table$PairKey
  )

  if (anyNA(pair_match)) {
    stop("Some temporal pairs were not found in iCAMP for ", basin_name)
  }

  temporal <- cbind(
    data.frame(
      Basin = basin_name,
      allowed_pairs[, c(
        "Sampling_Point",
        "sample_early",
        "sample_late",
        "Time_early",
        "Time_late",
        "Phase_early",
        "Phase_late"
      )],
      stringsAsFactors = FALSE
    ),
    process_table[
      pair_match,
      process_cols,
      drop = FALSE
    ]
  )

  temporal$Transition <- paste(
    temporal$Phase_early,
    temporal$Phase_late,
    sep = " -> "
  )

  dominant_index <- max.col(
    temporal[, process_cols, drop = FALSE],
    ties.method = "first"
  )

  temporal$Dominant_Process <- process_cols[dominant_index]

  write.table(
    temporal,
    file.path(
      final_dir,
      paste0("iCAMP_", basin_name, "_temporal_same_sampling_point.tsv")
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  temporal
}


temporal_icamp <- purrr::map2_dfr(
  final_settings$Basin,
  final_settings$BinSize,
  extract_temporal_icamp
)

write.table(
  temporal_icamp,
  file.path(table_dir, "iCAMP_temporal_same_sampling_point_ALL.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  temporal_icamp,
  file.path(output_dir, "iCAMP_temporal_same_sampling_point_ALL.rds")
)


# ============================================================================ #
# Overall process composition by basin
# ============================================================================ #

process_order_overall <- c(
  "Drift.and.Others",
  "Heterogeneous.Selection",
  "Homogeneous.Selection",
  "Dispersal.Limitation",
  "Homogenizing.Dispersal"
)

point_overall <- aggregate(
  temporal_icamp[, process_order_overall, drop = FALSE],
  by = list(
    Basin = temporal_icamp$Basin,
    Sampling_Point = temporal_icamp$Sampling_Point
  ),
  FUN = mean
)

basin_overall <- aggregate(
  point_overall[, process_order_overall, drop = FALSE],
  by = list(Basin = point_overall$Basin),
  FUN = mean
)

overall_long <- purrr::map_dfr(
  process_order_overall,
  function(p) {
    data.frame(
      Basin = basin_overall$Basin,
      Process = p,
      Importance = basin_overall[[p]] * 100,
      stringsAsFactors = FALSE
    )
  }
) |>
  mutate(
    Basin = factor(Basin, levels = c("HCP", "CHF", "YSB")),
    Process = factor(Process, levels = process_order_overall),
    Label = sprintf("%.1f%%", Importance)
  )

write.table(
  overall_long[, c("Basin", "Process", "Importance")],
  file.path(table_dir, "iCAMP_overall_process_composition.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

p_icamp_overall <- ggplot(
  overall_long,
  aes(
    x = 0.55,
    y = Importance,
    fill = Process
  )
) +
  geom_col(
    width = 1,
    position = position_stack(reverse = TRUE)
  ) +
  geom_vline(
    xintercept = 0.45,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey55"
  ) +
  geom_vline(
    xintercept = 1.12,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey55"
  ) +
  geom_vline(
    xintercept = 1.05,
    linewidth = 0.6,
    color = "black"
  ) +
  geom_text(
    aes(
      x = 1.28,
      label = Label,
      group = Process
    ),
    position = position_stack(
      vjust = 0.5,
      reverse = TRUE
    ),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  coord_polar(theta = "y") +
  facet_wrap(
    ~ Basin,
    nrow = 1,
    labeller = labeller(Basin = basin_labels)
  ) +
  xlim(0, 1.45) +
  scale_fill_manual(
    name = "Assembly process",
    values = process_colors,
    breaks = process_order_overall,
    labels = process_labels[process_order_overall]
  ) +
  theme_void(base_size = 12) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 14,
      margin = margin(t = 6, b = 8)
    ),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.spacing = grid::unit(1.2, "lines")
  )

ggsave(
  file.path(figure_dir, "Figure8_component_overall_process_composition.png"),
  p_icamp_overall,
  width = 10,
  height = 4.2,
  dpi = 300
)

ggsave(
  file.path(figure_dir, "Figure8_component_overall_process_composition.svg"),
  p_icamp_overall,
  width = 10,
  height = 4.2,
  device = grDevices::svg
)


# ============================================================================ #
# Temporal trajectory table
# ============================================================================ #

metadata_all <- purrr::map_dfr(
  final_settings$Basin,
  function(b) {
    md <- readRDS(
      file.path(
        output_dir,
        b,
        "input",
        paste0("metadata_", b, "_icamp.rds")
      )
    )

    data.frame(
      Basin = b,
      SampleID = rownames(md),
      Date = as.Date(as.character(md$Date)),
      Phase = as.character(md$Phase),
      stringsAsFactors = FALSE
    )
  }
)

metadata_key <- paste(
  metadata_all$Basin,
  metadata_all$SampleID,
  sep = "||"
)

early_key <- paste(
  temporal_icamp$Basin,
  temporal_icamp$sample_early,
  sep = "||"
)

late_key <- paste(
  temporal_icamp$Basin,
  temporal_icamp$sample_late,
  sep = "||"
)

idx_early <- match(early_key, metadata_key)
idx_late <- match(late_key, metadata_key)

if (anyNA(idx_early) || anyNA(idx_late)) {
  stop("Some temporal samples could not be matched to metadata dates.")
}

date_early <- metadata_all$Date[idx_early]
date_late <- metadata_all$Date[idx_late]

mid_date_numeric <- (
  as.numeric(date_early) +
    as.numeric(date_late)
) / 2

temporal_icamp$Days_after_rain <- (
  mid_date_numeric -
    as.numeric(rain_date)
)

temporal_icamp$Interval_gap <- (
  temporal_icamp$Time_late -
    temporal_icamp$Time_early
)

trajectory_long <- purrr::map_dfr(
  process_cols,
  function(p) {
    data.frame(
      Basin = temporal_icamp$Basin,
      Sampling_Point = temporal_icamp$Sampling_Point,
      Time_early = temporal_icamp$Time_early,
      Time_late = temporal_icamp$Time_late,
      Days_after_rain = temporal_icamp$Days_after_rain,
      Transition = temporal_icamp$Transition,
      Interval_gap = temporal_icamp$Interval_gap,
      Process = p,
      Importance = temporal_icamp[[p]] * 100,
      stringsAsFactors = FALSE
    )
  }
)

write.table(
  trajectory_long,
  file.path(table_dir, "iCAMP_temporal_trajectory_sampling_points.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================================ #
# Figure 8 temporal component
# ============================================================================ #

trajectory_long <- trajectory_long |>
  mutate(
    Basin = factor(Basin, levels = c("HCP", "CHF", "YSB")),
    Process = factor(Process, levels = process_cols),
    Process_Group = ifelse(
      Process %in% c(
        "Heterogeneous.Selection",
        "Homogeneous.Selection"
      ),
      "Selection",
      "Dispersal / Drift and Others"
    ),
    Process_Group = factor(
      Process_Group,
      levels = c(
        "Selection",
        "Dispersal / Drift and Others"
      )
    )
  )

trajectory_mean <- aggregate(
  Importance ~
    Basin +
    Time_early +
    Time_late +
    Days_after_rain +
    Process +
    Process_Group,
  data = trajectory_long,
  FUN = mean
)

trajectory_mean <- trajectory_mean[
  order(
    trajectory_mean$Basin,
    trajectory_mean$Process_Group,
    trajectory_mean$Process,
    trajectory_mean$Days_after_rain
  ),
  ,
  drop = FALSE
]

phase_boundary <- purrr::map_dfr(
  final_settings$Basin,
  function(b) {

    md <- metadata_all |>
      filter(Basin == b)

    last_wet <- max(
      md$Date[md$Phase == "Wet"]
    )

    first_dry <- min(
      md$Date[md$Phase == "Dry"]
    )

    boundary_date <- mean(
      as.numeric(c(last_wet, first_dry))
    )

    data.frame(
      Basin = b,
      Boundary_Days = boundary_date - as.numeric(rain_date)
    )
  }
) |>
  mutate(
    Basin = factor(Basin, levels = c("HCP", "CHF", "YSB"))
  )

write.table(
  phase_boundary,
  file.path(table_dir, "Figure8_wet_dry_boundaries.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

wet_background <- tidyr::crossing(
  Basin = phase_boundary$Basin,
  Process_Group = levels(trajectory_long$Process_Group)
) |>
  left_join(
    phase_boundary,
    by = "Basin"
  ) |>
  transmute(
    Basin,
    Process_Group = factor(
      Process_Group,
      levels = levels(trajectory_long$Process_Group)
    ),
    xmin = 0,
    xmax = Boundary_Days,
    ymin = -Inf,
    ymax = Inf
  )

max_day <- ceiling(
  max(trajectory_long$Days_after_rain, na.rm = TRUE) / 50
) * 50

p_icamp_temporal <- ggplot() +
  geom_rect(
    data = wet_background,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = "grey85",
    alpha = 0.35
  ) +
  geom_vline(
    data = phase_boundary,
    aes(xintercept = Boundary_Days),
    linetype = "dashed",
    linewidth = 0.55
  ) +
  geom_point(
    data = trajectory_long,
    aes(
      x = Days_after_rain,
      y = Importance,
      color = Process
    ),
    size = 1.6,
    alpha = 0.25
  ) +
  geom_line(
    data = trajectory_mean,
    aes(
      x = Days_after_rain,
      y = Importance,
      color = Process,
      group = Process
    ),
    linewidth = 0.9
  ) +
  geom_point(
    data = trajectory_mean,
    aes(
      x = Days_after_rain,
      y = Importance,
      color = Process
    ),
    size = 2.2
  ) +
  facet_grid(
    rows = vars(Process_Group),
    cols = vars(Basin),
    scales = "free_y",
    labeller = labeller(
      Basin = basin_labels,
      Process_Group = c(
        "Selection" = "Selection",
        "Dispersal / Drift and Others" = "Dispersal / Drift and Others"
      )
    ),
    axes = "all",
    axis.labels = "all"
  ) +
  scale_color_manual(
    name = "Assembly process",
    values = process_colors,
    breaks = names(process_labels),
    labels = process_labels
  ) +
  scale_x_continuous(
    limits = c(0, max_day),
    breaks = seq(0, max_day, by = 100),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = seq(0, 100, by = 20),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = "Days after rainfall",
    y = "Relative importance (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(
      face = "bold",
      size = 14,
      margin = margin(t = 6, b = 8)
    ),
    strip.text.y = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9.5),
    panel.spacing.x = grid::unit(1, "lines"),
    panel.spacing.y = grid::unit(0.8, "lines")
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        alpha = 1,
        size = 2.4
      )
    )
  )

ggsave(
  file.path(figure_dir, "Figure8_component_temporal_trajectory.png"),
  p_icamp_temporal,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(figure_dir, "Figure8_component_temporal_trajectory.svg"),
  p_icamp_temporal,
  width = 12,
  height = 7,
  device = grDevices::svg
)


# ============================================================================ #
# Wet-Dry summaries | sampling point as replicate
# ============================================================================ #

within_phase <- temporal_icamp[
  temporal_icamp$Transition %in% c(
    "Wet -> Wet",
    "Dry -> Dry"
  ),
  ,
  drop = FALSE
]

within_phase$Phase <- ifelse(
  within_phase$Transition == "Wet -> Wet",
  "Wet",
  "Dry"
)

point_phase_means <- aggregate(
  within_phase[, process_cols, drop = FALSE],
  by = list(
    Basin = within_phase$Basin,
    Sampling_Point = within_phase$Sampling_Point,
    Phase = within_phase$Phase
  ),
  FUN = mean
)

write.table(
  point_phase_means,
  file.path(table_dir, "iCAMP_WetDry_means_by_SamplingPoint.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

summary_mean <- aggregate(
  point_phase_means[, process_cols, drop = FALSE],
  by = list(
    Basin = point_phase_means$Basin,
    Phase = point_phase_means$Phase
  ),
  FUN = mean
)

summary_sd <- aggregate(
  point_phase_means[, process_cols, drop = FALSE],
  by = list(
    Basin = point_phase_means$Basin,
    Phase = point_phase_means$Phase
  ),
  FUN = sd
)

colnames(summary_sd)[-(1:2)] <- paste0(
  colnames(summary_sd)[-(1:2)],
  "_SD"
)

wetdry_final_summary <- merge(
  summary_mean,
  summary_sd,
  by = c("Basin", "Phase")
)

write.table(
  wetdry_final_summary,
  file.path(table_dir, "iCAMP_WetDry_final_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================================ #
# Dry minus Wet changes by sampling point
# ============================================================================ #

point_groups <- split(
  point_phase_means,
  interaction(
    point_phase_means$Basin,
    point_phase_means$Sampling_Point,
    drop = TRUE
  )
)

point_deltas <- lapply(
  point_groups,
  function(x) {

    if (!all(c("Wet", "Dry") %in% x$Phase)) {
      return(NULL)
    }

    wet <- x[
      x$Phase == "Wet",
      process_cols,
      drop = FALSE
    ]

    dry <- x[
      x$Phase == "Dry",
      process_cols,
      drop = FALSE
    ]

    delta <- dry - wet

    data.frame(
      Basin = x$Basin[1],
      Sampling_Point = x$Sampling_Point[1],
      delta,
      stringsAsFactors = FALSE
    )
  }
) |>
  dplyr::bind_rows()

write.table(
  point_deltas,
  file.path(table_dir, "iCAMP_Dry_minus_Wet_by_SamplingPoint.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================================ #
# Exact two-sided sign-flip tests
# ============================================================================ #

point_phase_stats <- point_phase_means

point_phase_stats$Deterministic <- (
  point_phase_stats$Heterogeneous.Selection +
    point_phase_stats$Homogeneous.Selection
)

point_phase_stats$Stochastic <- (
  point_phase_stats$Dispersal.Limitation +
    point_phase_stats$Homogenizing.Dispersal +
    point_phase_stats$Drift.and.Others
)

metrics <- c(
  process_cols,
  "Deterministic",
  "Stochastic"
)

point_ids <- unique(
  paste(
    point_phase_stats$Basin,
    point_phase_stats$Sampling_Point,
    sep = "||"
  )
)

wetdry_delta <- lapply(
  point_ids,
  function(id) {

    idx <- paste(
      point_phase_stats$Basin,
      point_phase_stats$Sampling_Point,
      sep = "||"
    ) == id

    x <- point_phase_stats[idx, , drop = FALSE]

    if (!all(c("Wet", "Dry") %in% x$Phase)) {
      return(NULL)
    }

    wet <- x[x$Phase == "Wet", metrics, drop = FALSE]
    dry <- x[x$Phase == "Dry", metrics, drop = FALSE]

    delta <- dry - wet

    data.frame(
      Basin = x$Basin[1],
      Sampling_Point = x$Sampling_Point[1],
      delta,
      stringsAsFactors = FALSE
    )
  }
) |>
  dplyr::bind_rows()

exact_signflip <- function(x) {

  x <- as.numeric(x)
  x <- x[!is.na(x)]

  n <- length(x)
  observed <- mean(x)

  sign_matrix <- expand.grid(
    rep(list(c(-1, 1)), n)
  ) |>
    as.matrix()

  permuted_means <- apply(
    sign_matrix,
    1,
    function(s) mean(x * s)
  )

  p_value <- mean(
    abs(permuted_means) >=
      abs(observed) - 1e-12
  )

  c(
    N = n,
    WetDry_delta_mean = observed,
    WetDry_delta_median = median(x),
    Positive = sum(x > 0),
    Negative = sum(x < 0),
    Zero = sum(x == 0),
    P_perm = p_value
  )
}

wetdry_stats <- purrr::map_dfr(
  metrics,
  function(metric) {
    result <- exact_signflip(
      wetdry_delta[[metric]]
    )

    data.frame(
      Metric = metric,
      t(result),
      stringsAsFactors = FALSE
    )
  }
)

numeric_cols <- setdiff(
  colnames(wetdry_stats),
  "Metric"
)

wetdry_stats[, numeric_cols] <- lapply(
  wetdry_stats[, numeric_cols, drop = FALSE],
  as.numeric
)

wetdry_stats$P_BH <- NA_real_

process_rows <- wetdry_stats$Metric %in% process_cols

wetdry_stats$P_BH[process_rows] <- p.adjust(
  wetdry_stats$P_perm[process_rows],
  method = "BH"
)

aggregate_rows <- wetdry_stats$Metric %in% c(
  "Deterministic",
  "Stochastic"
)

wetdry_stats$P_BH[aggregate_rows] <- p.adjust(
  wetdry_stats$P_perm[aggregate_rows],
  method = "BH"
)

wetdry_stats$Mean_change_pp <- (
  wetdry_stats$WetDry_delta_mean * 100
)

wetdry_stats$Median_change_pp <- (
  wetdry_stats$WetDry_delta_median * 100
)

metric_labels <- c(
  process_labels,
  "Deterministic" = "Deterministic",
  "Stochastic" = "Stochastic"
)

wetdry_stats$Metric_label <- metric_labels[
  wetdry_stats$Metric
]

wetdry_stats <- wetdry_stats[, c(
  "Metric",
  "Metric_label",
  "N",
  "Mean_change_pp",
  "Median_change_pp",
  "Positive",
  "Negative",
  "Zero",
  "P_perm",
  "P_BH"
)]

write.table(
  wetdry_stats,
  file.path(table_dir, "Table_Supplementary9_iCAMP_WetDry_statistics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

openxlsx::write.xlsx(
  list(
    Overall_processes = overall_long,
    WetDry_summary = wetdry_final_summary,
    SamplingPoint_phase_means = point_phase_means,
    Dry_minus_Wet = wetdry_delta,
    Exact_signflip_tests = wetdry_stats
  ),
  file.path(table_dir, "Table_Supplementary9_iCAMP.xlsx"),
  overwrite = TRUE
)


# ============================================================================ #
# Final summary
# ============================================================================ #

cat("\niCAMP analysis completed.\n\n")
print(final_summary)

cat("\nOverall basin process composition (%):\n")
print(overall_long)

cat("\nWet-Dry exact sign-flip tests:\n")
print(wetdry_stats)
