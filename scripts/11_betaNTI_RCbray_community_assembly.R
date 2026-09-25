# ============================================================================ #
# 11_betaNTI_RCbray_community_assembly.R
# ============================================================================ #
#
# Purpose:
# Reproduce the whole-community betaNTI-RCbray assembly framework used in
# Figure 9, Figure Supplementary 8, and associated supplementary tables.
#
# Workflow:
#
#   final non-rarefied phyloseq object
#       -> basin-specific ASV filtering (>=5 total reads per basin)
#       -> short-distance phylogenetic-signal diagnostics
#       -> abundance-weighted betaMNTD
#       -> betaNTI from 999 phylogenetic-label permutations
#       -> temporal comparisons within permanent sampling points
#       -> RCbray for pairs with |betaNTI| <= 2 using 999 null communities
#       -> ecological-process classification
#       -> consecutive temporal transitions
#       -> Figure 9A-C
#       -> Figure Supplementary 8
#
# Final process thresholds:
#
#   betaNTI > +2                  -> Variable selection
#   betaNTI < -2                  -> Homogeneous selection
#   |betaNTI| <= 2 & RCbray > .95 -> Dispersal limitation
#   |betaNTI| <= 2 & RCbray < -.95 -> Homogenizing dispersal
#   |betaNTI| <= 2 & |RCbray| <= .95 -> Undominated
#
# Important:
# - This analysis uses NON-RAREFIED counts.
# - betaNTI and RCbray null models use 999 permutations and seed 123.
# - Null-model functions include checkpointing because they are computationally
#   expensive.
# - Figure 9 analyses use consecutive comparisons within the same permanent
#   sampling point.
# - Wet-Dry boundary-crossing transitions are excluded from phase summaries.
#
# ============================================================================ #


# ============================================================================ #
# Packages
# ============================================================================ #

library(phyloseq)
library(ape)
library(picante)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(cowplot)
library(gghalves)
library(scales)
library(lme4)
library(emmeans)
library(openxlsx)


# ============================================================================ #
# Paths and reproducibility
# ============================================================================ #

phyloseq_path <- file.path(
  "output_phyloseq_tree",
  "ps_mdecon_v2_FTree.rds"
)

output_dir <- "output_betaNTI_RCbray"
object_dir <- file.path(output_dir, "objects")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
diagnostic_dir <- file.path(output_dir, "diagnostics")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(object_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(diagnostic_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(123)


# ============================================================================ #
# Shared settings
# ============================================================================ #

basins <- c(
  "HCP",
  "CHF",
  "YSB"
)

basin_labels <- c(
  "HCP" = "Herradura Clay Pan",
  "CHF" = "Ckoirama Halite Field",
  "YSB" = "Yungay Station Basin"
)

basin_colors <- c(
  "HCP" = "#1b9e77",
  "CHF" = "#d95f02",
  "YSB" = "#7570b3"
)

phase_colors <- c(
  "Wet" = "#4C78A8",
  "Dry" = "#E3B967"
)

assembly_colors <- c(
  "Variable selection" = "#ca0020",
  "Homogeneous selection" = "#f4a582",
  "Undominated or Drift" = "#f7f7f7",
  "Homogenizing dispersal" = "#92c5de",
  "Dispersal limitation" = "#0571b0"
)

process_order <- c(
  "Variable selection",
  "Homogeneous selection",
  "Undominated or Drift",
  "Homogenizing dispersal",
  "Dispersal limitation"
)

env_vars <- c(
  "Wc",
  "Salinity",
  "pH",
  "ORP"
)

N_PERM <- 999
NULL_SEED <- 123
CHECKPOINT_EVERY <- 25


# ============================================================================ #
# Load final phyloseq object
# ============================================================================ #

ps <- readRDS(
  phyloseq_path
)

if (!"Basin" %in% sample_variables(ps)) {
  stop(
    "The phyloseq metadata must contain a column named 'Basin' ",
    "with values HCP, CHF, and YSB."
  )
}

if (!all(basins %in% unique(as.character(sample_data(ps)$Basin)))) {
  stop(
    "Not all expected basin codes HCP, CHF, and YSB were found."
  )
}


# ============================================================================ #
# Basin-specific ASV filtering
# ============================================================================ #

filter_basin <- function(
    ps_object,
    basin_name,
    min_reads = 5
) {

  meta_b <- data.frame(
    sample_data(ps_object)
  )

  keep_samples <- rownames(meta_b)[
    as.character(meta_b$Basin) ==
      basin_name
  ]

  ps_b <- prune_samples(
    keep_samples,
    ps_object
  )

  ps_b <- prune_taxa(
    taxa_sums(ps_b) > 0,
    ps_b
  )

  keep_taxa <- taxa_sums(ps_b) >=
    min_reads

  ps_b <- prune_taxa(
    keep_taxa,
    ps_b
  )

  ps_b <- prune_samples(
    sample_sums(ps_b) > 0,
    ps_b
  )

  ps_b
}

ps_list <- list(
  HCP = filter_basin(
    ps,
    "HCP",
    min_reads = 5
  ),
  CHF = filter_basin(
    ps,
    "CHF",
    min_reads = 5
  ),
  YSB = filter_basin(
    ps,
    "YSB",
    min_reads = 5
  )
)

filter_summary <- purrr::imap_dfr(
  ps_list,
  function(
      ps_b,
      basin_name
  ) {
    tibble(
      Basin = basin_name,
      Samples = nsamples(ps_b),
      ASVs = ntaxa(ps_b),
      Reads = sum(sample_sums(ps_b)),
      Tree_tips = length(
        phy_tree(ps_b)$tip.label
      )
    )
  }
)

openxlsx::write.xlsx(
  filter_summary,
  file.path(
    table_dir,
    "basin_ASV_filter_summary.xlsx"
  ),
  overwrite = TRUE
)

saveRDS(
  ps_list,
  file.path(
    object_dir,
    "basin_phyloseq_filtered_min5reads.rds"
  )
)


# ============================================================================ #
# Phylogenetic distance matrices
# ============================================================================ #

get_phylogenetic_distances <- function(
    ps_object
) {

  tree_b <- phy_tree(
    ps_object
  )

  phylo_dist <- ape::cophenetic.phylo(
    tree_b
  )

  taxa_b <- taxa_names(
    ps_object
  )

  phylo_dist[
    taxa_b,
    taxa_b,
    drop = FALSE
  ]
}

phylo_list <- purrr::map(
  ps_list,
  get_phylogenetic_distances
)

saveRDS(
  phylo_list,
  file.path(
    object_dir,
    "basin_phylogenetic_distance_matrices.rds"
  )
)


# ============================================================================ #
# Short-distance phylogenetic-signal diagnostics
# ============================================================================ #

calculate_niche_optima <- function(
    ps_object,
    environmental_variables
) {

  otu_b <- as(
    otu_table(ps_object),
    "matrix"
  )

  if (taxa_are_rows(ps_object)) {
    otu_b <- t(
      otu_b
    )
  }

  md <- data.frame(
    sample_data(ps_object),
    check.names = FALSE
  )

  md <- md[
    rownames(otu_b),
    ,
    drop = FALSE
  ]

  otu_rel <- sweep(
    otu_b,
    1,
    rowSums(otu_b),
    "/"
  )

  result <- data.frame(
    ASV = colnames(otu_rel),
    stringsAsFactors = FALSE
  )

  for (v in environmental_variables) {

    env_value <- as.numeric(
      md[[v]]
    )

    valid_samples <- !is.na(env_value) &
      is.finite(env_value)

    weights <- otu_rel[
      valid_samples,
      ,
      drop = FALSE
    ]

    env_valid <- env_value[
      valid_samples
    ]

    weighted_values <- sweep(
      weights,
      1,
      env_valid,
      "*"
    )

    denominator <- colSums(
      weights
    )

    optimum <- colSums(
      weighted_values
    ) /
      denominator

    optimum[
      denominator == 0
    ] <- NA_real_

    result[[v]] <- optimum
  }

  result
}


calculate_niche_distances <- function(
    niche_df,
    environmental_variables
) {

  niche_matrix <- niche_df[
    ,
    environmental_variables,
    drop = FALSE
  ]

  rownames(niche_matrix) <- niche_df$ASV

  complete_taxa <- complete.cases(
    niche_matrix
  )

  niche_matrix <- niche_matrix[
    complete_taxa,
    ,
    drop = FALSE
  ]

  niche_scaled <- scale(
    niche_matrix
  )

  combined_dist <- as.matrix(
    dist(
      niche_scaled,
      method = "euclidean"
    )
  )

  individual_dist <- lapply(
    environmental_variables,
    function(v) {

      as.matrix(
        dist(
          niche_scaled[, v],
          method = "euclidean"
        )
      )
    }
  )

  names(
    individual_dist
  ) <- environmental_variables

  list(
    scaled_optima = niche_scaled,
    combined = combined_dist,
    individual = individual_dist
  )
}


run_mantel_correlog <- function(
    niche_matrix,
    phylo_matrix,
    n_classes = 50,
    permutations = 999
) {

  common_taxa <- intersect(
    rownames(niche_matrix),
    rownames(phylo_matrix)
  )

  niche_matrix <- niche_matrix[
    common_taxa,
    common_taxa,
    drop = FALSE
  ]

  phylo_matrix <- phylo_matrix[
    common_taxa,
    common_taxa,
    drop = FALSE
  ]

  vegan::mantel.correlog(
    D.eco = as.dist(niche_matrix),
    D.geo = as.dist(phylo_matrix),
    n.class = n_classes,
    cutoff = FALSE,
    r.type = "pearson",
    nperm = permutations,
    mult = "holm",
    progressive = TRUE
  )
}


niche_list <- purrr::map(
  ps_list,
  calculate_niche_optima,
  environmental_variables = env_vars
)

niche_distance_list <- purrr::map(
  niche_list,
  calculate_niche_distances,
  environmental_variables = env_vars
)

mantel_combined <- list()
mantel_individual <- list()

for (b in basins) {

  set.seed(NULL_SEED)

  mantel_combined[[b]] <- run_mantel_correlog(
    niche_distance_list[[b]]$combined,
    phylo_list[[b]],
    n_classes = 50,
    permutations = N_PERM
  )

  mantel_individual[[b]] <- list()

  for (v in env_vars) {

    set.seed(NULL_SEED)

    mantel_individual[[b]][[v]] <-
      run_mantel_correlog(
        niche_distance_list[[b]]$individual[[v]],
        phylo_list[[b]],
        n_classes = 50,
        permutations = N_PERM
      )
  }
}


extract_mantel_classes <- function(
    mantel_object,
    basin_name,
    variable_name
) {

  x <- as.data.frame(
    mantel_object$mantel.res
  )

  tibble(
    Basin = basin_name,
    Variable = variable_name,
    Class = rownames(x),
    Phylo_distance = x[["class.index"]],
    Mantel_r = x[["Mantel.cor"]],
    P_value = x[["Pr(Mantel)"]],
    P_Holm = x[["Pr(corrected)"]]
  )
}

phylo_signal_table <- bind_rows(
  lapply(
    basins,
    function(b) {

      bind_rows(
        extract_mantel_classes(
          mantel_combined[[b]],
          b,
          "Combined"
        ),
        bind_rows(
          lapply(
            env_vars,
            function(v) {
              extract_mantel_classes(
                mantel_individual[[b]][[v]],
                b,
                v
              )
            }
          )
        )
      )
    }
  )
)

write.table(
  phylo_signal_table,
  file.path(
    table_dir,
    "phylogenetic_signal_mantel_correlograms.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  list(
    combined = mantel_combined,
    individual = mantel_individual
  ),
  file.path(
    object_dir,
    "phylogenetic_signal_mantel_correlograms.rds"
  )
)


# ============================================================================ #
# Prepare abundance-weighted betaMNTD inputs
# ============================================================================ #

prepare_bmntd_input <- function(
    ps_object,
    phylo_matrix
) {

  comm <- as(
    otu_table(ps_object),
    "matrix"
  )

  if (taxa_are_rows(ps_object)) {
    comm <- t(
      comm
    )
  }

  comm <- sweep(
    comm,
    1,
    rowSums(comm),
    "/"
  )

  taxa <- colnames(
    comm
  )

  phylo_matrix <- phylo_matrix[
    taxa,
    taxa,
    drop = FALSE
  ]

  stopifnot(
    identical(
      colnames(comm),
      rownames(phylo_matrix)
    )
  )

  list(
    comm = comm,
    phylo = phylo_matrix
  )
}

bmntd_inputs <- purrr::map2(
  ps_list,
  phylo_list,
  prepare_bmntd_input
)


# ============================================================================ #
# Observed betaMNTD
# ============================================================================ #

calculate_observed_bmntd <- function(
    bmntd_input
) {

  picante::comdistnt(
    comm = bmntd_input$comm,
    dis = bmntd_input$phylo,
    abundance.weighted = TRUE,
    exclude.conspecifics = FALSE
  )
}

betaMNTD_observed <- purrr::map(
  bmntd_inputs,
  calculate_observed_bmntd
)

saveRDS(
  betaMNTD_observed,
  file.path(
    object_dir,
    "betaMNTD_observed_all_basins.rds"
  )
)


# ============================================================================ #
# betaNTI null model with checkpointing
# ============================================================================ #

run_bnti_null_checkpoint <- function(
    bmntd_input,
    observed_bmntd,
    basin_name,
    nperm = 999,
    checkpoint_every = 25,
    seed = 123
) {

  comm <- bmntd_input$comm
  phylo <- bmntd_input$phylo

  taxa <- colnames(
    comm
  )

  obs_vec <- as.numeric(
    observed_bmntd
  )

  n_pairs <- length(
    obs_vec
  )

  checkpoint_file <- file.path(
    object_dir,
    paste0(
      "betaNTI_",
      basin_name,
      "_",
      nperm,
      "perm_CHECKPOINT.rds"
    )
  )

  final_file <- file.path(
    object_dir,
    paste0(
      "betaNTI_",
      basin_name,
      "_",
      nperm,
      "perm_FINAL.rds"
    )
  )

  if (file.exists(final_file)) {
    message(
      basin_name,
      ": loading existing final betaNTI object."
    )
    return(
      readRDS(final_file)
    )
  }

  if (file.exists(checkpoint_file)) {

    checkpoint <- readRDS(
      checkpoint_file
    )

    null_matrix <- checkpoint$null_matrix
    completed <- checkpoint$completed

    stopifnot(
      nrow(null_matrix) == n_pairs,
      ncol(null_matrix) == nperm
    )

    assign(
      ".Random.seed",
      checkpoint$rng_state,
      envir = .GlobalEnv
    )

  } else {

    set.seed(seed)

    null_matrix <- matrix(
      NA_real_,
      nrow = n_pairs,
      ncol = nperm
    )

    completed <- 0L
  }

  if (completed < nperm) {

    for (
      i in seq.int(
        completed + 1,
        nperm
      )
    ) {

      perm_index <- sample.int(
        length(taxa)
      )

      phylo_null <- phylo[
        perm_index,
        perm_index,
        drop = FALSE
      ]

      rownames(phylo_null) <- taxa
      colnames(phylo_null) <- taxa

      null_dist <- picante::comdistnt(
        comm = comm,
        dis = phylo_null,
        abundance.weighted = TRUE,
        exclude.conspecifics = FALSE
      )

      null_matrix[, i] <- as.numeric(
        null_dist
      )

      if (
        i %% checkpoint_every == 0 ||
          i == nperm
      ) {

        saveRDS(
          list(
            basin = basin_name,
            nperm = nperm,
            completed = i,
            null_matrix = null_matrix,
            rng_state = get(
              ".Random.seed",
              envir = .GlobalEnv
            )
          ),
          checkpoint_file
        )

        message(
          basin_name,
          ": betaNTI checkpoint ",
          i,
          "/",
          nperm
        )
      }
    }
  }

  null_mean <- rowMeans(
    null_matrix
  )

  null_sd <- apply(
    null_matrix,
    1,
    sd
  )

  betaNTI <- (
    obs_vec -
      null_mean
  ) /
    null_sd

  obs_matrix <- as.matrix(
    observed_bmntd
  )

  pair_index <- which(
    lower.tri(obs_matrix),
    arr.ind = TRUE
  )

  stopifnot(
    isTRUE(
      all.equal(
        obs_vec,
        as.numeric(
          obs_matrix[
            pair_index
          ]
        )
      )
    )
  )

  pairwise <- tibble(
    Sample_1 = rownames(obs_matrix)[
      pair_index[, 1]
    ],
    Sample_2 = colnames(obs_matrix)[
      pair_index[, 2]
    ],
    betaMNTD_observed = obs_vec,
    betaMNTD_null_mean = null_mean,
    betaMNTD_null_sd = null_sd,
    betaNTI = betaNTI
  )

  result <- list(
    Basin = basin_name,
    permutations = nperm,
    observed_betaMNTD = observed_bmntd,
    null_betaMNTD = null_matrix,
    pairwise = pairwise
  )

  saveRDS(
    result,
    final_file
  )

  write.table(
    pairwise,
    file.path(
      table_dir,
      paste0(
        "betaNTI_pairwise_",
        basin_name,
        "_",
        nperm,
        "perm.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  result
}


# ============================================================================ #
# Run betaNTI for all basins
# ============================================================================ #

betaNTI_results <- list()

for (b in basins) {

  betaNTI_results[[b]] <-
    run_bnti_null_checkpoint(
      bmntd_input = bmntd_inputs[[b]],
      observed_bmntd = betaMNTD_observed[[b]],
      basin_name = b,
      nperm = N_PERM,
      checkpoint_every = CHECKPOINT_EVERY,
      seed = NULL_SEED
    )
}


# ============================================================================ #
# Add metadata to all betaNTI pairs
# ============================================================================ #

build_bnti_pair_metadata <- function(
    bnti_result,
    ps_object,
    basin_name
) {

  meta <- data.frame(
    sample_data(ps_object),
    check.names = FALSE
  )

  meta$Sample <- rownames(
    meta
  )

  meta <- meta |>
    select(
      Sample,
      Sampling_Point,
      Time,
      Time_Point,
      Time_Interval,
      Date,
      Phase,
      Wc,
      Salinity,
      pH,
      ORP
    )

  meta$Date <- as.Date(
    as.character(
      meta$Date
    )
  )

  meta_1 <- meta |>
    rename(
      Sample_1 = Sample
    ) |>
    rename_with(
      ~ paste0(
        .x,
        "_1"
      ),
      -Sample_1
    )

  meta_2 <- meta |>
    rename(
      Sample_2 = Sample
    ) |>
    rename_with(
      ~ paste0(
        .x,
        "_2"
      ),
      -Sample_2
    )

  bnti_result$pairwise |>
    left_join(
      meta_1,
      by = "Sample_1"
    ) |>
    left_join(
      meta_2,
      by = "Sample_2"
    ) |>
    mutate(
      Basin = basin_name,
      Same_sampling_point =
        Sampling_Point_1 ==
        Sampling_Point_2,
      Time_gap = abs(
        Time_Point_2 -
          Time_Point_1
      )
    )
}

bnti_meta <- purrr::map2(
  betaNTI_results,
  ps_list,
  function(
      bnti_result,
      ps_object
  ) {
    basin_name <- bnti_result$Basin

    build_bnti_pair_metadata(
      bnti_result,
      ps_object,
      basin_name
    )
  }
)

bnti_all_pairs <- bind_rows(
  bnti_meta
)

saveRDS(
  bnti_all_pairs,
  file.path(
    object_dir,
    "betaNTI_all_pairwise_comparisons.rds"
  )
)

write.table(
  bnti_all_pairs,
  file.path(
    table_dir,
    "betaNTI_all_pairwise_comparisons.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================================ #
# Within-sampling-point temporal betaNTI pairs
# ============================================================================ #

prepare_temporal_pairs <- function(
    x
) {

  x |>
    filter(
      Same_sampling_point
    ) |>
    mutate(
      Early_Sample = if_else(
        Time_Point_1 <
          Time_Point_2,
        Sample_1,
        Sample_2
      ),
      Late_Sample = if_else(
        Time_Point_1 <
          Time_Point_2,
        Sample_2,
        Sample_1
      ),
      Early_Time = if_else(
        Time_Point_1 <
          Time_Point_2,
        Time_1,
        Time_2
      ),
      Late_Time = if_else(
        Time_Point_1 <
          Time_Point_2,
        Time_2,
        Time_1
      ),
      Early_Time_Point = pmin(
        Time_Point_1,
        Time_Point_2
      ),
      Late_Time_Point = pmax(
        Time_Point_1,
        Time_Point_2
      ),
      Early_Time_Interval = if_else(
        Time_Point_1 <
          Time_Point_2,
        Time_Interval_1,
        Time_Interval_2
      ),
      Late_Time_Interval = if_else(
        Time_Point_1 <
          Time_Point_2,
        Time_Interval_2,
        Time_Interval_1
      ),
      Early_Date = pmin(
        Date_1,
        Date_2
      ),
      Late_Date = pmax(
        Date_1,
        Date_2
      ),
      Delta_days = as.numeric(
        Late_Date -
          Early_Date
      ),
      Consecutive =
        (
          Late_Time_Point -
            Early_Time_Point
        ) ==
        1
    )
}

bnti_temporal <- purrr::map(
  bnti_meta,
  prepare_temporal_pairs
)

bnti_temporal_all <- bind_rows(
  bnti_temporal
)

saveRDS(
  bnti_temporal_all,
  file.path(
    object_dir,
    "betaNTI_temporal_pairs_ALL.rds"
  )
)


# ============================================================================ #
# Prepare RCbray input
# ============================================================================ #

prepare_rcbray_input <- function(
    ps_object,
    bnti_temporal_data,
    basin_name
) {

  comm <- as(
    otu_table(ps_object),
    "matrix"
  )

  if (taxa_are_rows(ps_object)) {
    comm <- t(
      comm
    )
  }

  target_pairs <- bnti_temporal_data |>
    filter(
      abs(betaNTI) <=
        2
    )

  list(
    Basin = basin_name,
    comm = comm,
    target_pairs = target_pairs,
    richness = rowSums(
      comm > 0
    ),
    total_abundance = rowSums(
      comm
    ),
    occurrence = colSums(
      comm > 0
    ),
    regional_abundance = colSums(
      comm
    )
  )
}

rc_inputs <- list()

for (b in basins) {
  rc_inputs[[b]] <-
    prepare_rcbray_input(
      ps_list[[b]],
      bnti_temporal[[b]],
      b
    )
}


# ============================================================================ #
# RCbray null model with checkpointing
# ============================================================================ #

run_rcbray_checkpoint <- function(
    rc_input,
    basin_name,
    nperm = 999,
    checkpoint_every = 25,
    seed = 123
) {

  comm <- rc_input$comm
  target <- rc_input$target_pairs

  richness <- rc_input$richness
  total_abundance <- rc_input$total_abundance
  occurrence <- rc_input$occurrence
  regional_abundance <- rc_input$regional_abundance

  nsamp <- nrow(
    comm
  )

  ntaxa <- ncol(
    comm
  )

  npairs <- nrow(
    target
  )

  sample_names <- rownames(
    comm
  )

  taxa_names_local <- colnames(
    comm
  )

  idx1 <- match(
    target$Sample_1,
    sample_names
  )

  idx2 <- match(
    target$Sample_2,
    sample_names
  )

  stopifnot(
    !anyNA(idx1),
    !anyNA(idx2)
  )

  observed_matrix <- as.matrix(
    vegan::vegdist(
      comm,
      method = "bray"
    )
  )

  observed_bray <- observed_matrix[
    cbind(
      idx1,
      idx2
    )
  ]

  checkpoint_file <- file.path(
    object_dir,
    paste0(
      "RCbray_",
      basin_name,
      "_",
      nperm,
      "perm_CHECKPOINT.rds"
    )
  )

  final_file <- file.path(
    object_dir,
    paste0(
      "RCbray_",
      basin_name,
      "_",
      nperm,
      "perm_FINAL.rds"
    )
  )

  if (file.exists(final_file)) {

    message(
      basin_name,
      ": loading existing final RCbray object."
    )

    return(
      readRDS(final_file)
    )
  }

  if (file.exists(checkpoint_file)) {

    checkpoint <- readRDS(
      checkpoint_file
    )

    null_bray <- checkpoint$null_bray
    completed <- checkpoint$completed

    assign(
      ".Random.seed",
      checkpoint$rng_state,
      envir = .GlobalEnv
    )

  } else {

    set.seed(
      seed
    )

    null_bray <- matrix(
      NA_real_,
      nrow = npairs,
      ncol = nperm
    )

    completed <- 0L
  }

  if (completed < nperm) {

    for (
      iteration in seq.int(
        completed + 1,
        nperm
      )
    ) {

      null_comm <- matrix(
        0,
        nrow = nsamp,
        ncol = ntaxa,
        dimnames = list(
          sample_names,
          taxa_names_local
        )
      )

      for (s in seq_len(nsamp)) {

        S <- richness[s]
        N <- total_abundance[s]

        selected <- sample.int(
          ntaxa,
          size = S,
          replace = FALSE,
          prob = occurrence
        )

        null_comm[
          s,
          selected
        ] <- 1

        remaining <- N -
          S

        if (remaining > 0) {

          extra_counts <- as.vector(
            rmultinom(
              n = 1,
              size = remaining,
              prob = regional_abundance[
                selected
              ]
            )
          )

          null_comm[
            s,
            selected
          ] <- null_comm[
            s,
            selected
          ] +
            extra_counts
        }
      }

      null_matrix <- as.matrix(
        vegan::vegdist(
          null_comm,
          method = "bray"
        )
      )

      null_bray[
        ,
        iteration
      ] <- null_matrix[
        cbind(
          idx1,
          idx2
        )
      ]

      if (
        iteration %% checkpoint_every == 0 ||
          iteration == nperm
      ) {

        saveRDS(
          list(
            Basin = basin_name,
            nperm = nperm,
            completed = iteration,
            null_bray = null_bray,
            rng_state = get(
              ".Random.seed",
              envir = .GlobalEnv
            )
          ),
          checkpoint_file
        )

        message(
          basin_name,
          ": RCbray checkpoint ",
          iteration,
          "/",
          nperm
        )
      }
    }
  }

  num_less <- rowSums(
    null_bray <
      observed_bray
  )

  num_equal <- rowSums(
    null_bray ==
      observed_bray
  )

  RCbray <- (
    (
      num_less +
        num_equal /
        2
    ) /
      nperm -
      0.5
  ) *
    2

  pairwise <- target |>
    mutate(
      Bray_observed = observed_bray,
      RCbray = RCbray
    )

  result <- list(
    Basin = basin_name,
    permutations = nperm,
    observed_Bray = observed_bray,
    null_Bray = null_bray,
    pairwise = pairwise
  )

  saveRDS(
    result,
    final_file
  )

  write.table(
    pairwise,
    file.path(
      table_dir,
      paste0(
        "RCbray_temporal_",
        basin_name,
        "_",
        nperm,
        "perm.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  result
}


# ============================================================================ #
# Run RCbray
# ============================================================================ #

RCbray_results <- list()

for (b in basins) {

  RCbray_results[[b]] <-
    run_rcbray_checkpoint(
      rc_input = rc_inputs[[b]],
      basin_name = b,
      nperm = N_PERM,
      checkpoint_every = CHECKPOINT_EVERY,
      seed = NULL_SEED
    )
}


# ============================================================================ #
# Combine betaNTI and RCbray
# ============================================================================ #

rcbray_all <- bind_rows(
  lapply(
    basins,
    function(b) {

      RCbray_results[[b]]$pairwise |>
        select(
          Basin,
          Sample_1,
          Sample_2,
          RCbray
        )
    }
  )
)

assembly_process_all <- bnti_temporal_all |>
  left_join(
    rcbray_all,
    by = c(
      "Basin",
      "Sample_1",
      "Sample_2"
    )
  ) |>
  mutate(
    Assembly_Process = case_when(
      betaNTI > 2 ~
        "Variable selection",
      betaNTI < -2 ~
        "Homogeneous selection",
      abs(betaNTI) <= 2 &
        RCbray > 0.95 ~
        "Dispersal limitation",
      abs(betaNTI) <= 2 &
        RCbray < -0.95 ~
        "Homogenizing dispersal",
      abs(betaNTI) <= 2 &
        abs(RCbray) <= 0.95 ~
        "Undominated",
      TRUE ~ NA_character_
    )
  )

assembly_process_all$Assembly_Process <- factor(
  assembly_process_all$Assembly_Process,
  levels = c(
    "Variable selection",
    "Homogeneous selection",
    "Dispersal limitation",
    "Homogenizing dispersal",
    "Undominated"
  )
)

saveRDS(
  assembly_process_all,
  file.path(
    object_dir,
    "assembly_process_temporal_ALL.rds"
  )
)


# ============================================================================ #
# Consecutive temporal transitions
# ============================================================================ #

assembly_consecutive <- assembly_process_all |>
  filter(
    Consecutive
  ) |>
  mutate(
    Sampling_Point = Sampling_Point_1,
    Early_Phase = if_else(
      Time_Point_1 <
        Time_Point_2,
      as.character(Phase_1),
      as.character(Phase_2)
    ),
    Late_Phase = if_else(
      Time_Point_1 <
        Time_Point_2,
      as.character(Phase_2),
      as.character(Phase_1)
    ),
    Transition = paste0(
      Early_Time,
      " → ",
      Late_Time
    ),
    Phase_Transition = paste0(
      Early_Phase,
      " → ",
      Late_Phase
    ),
    Phase_State = case_when(
      Early_Phase == "Wet" &
        Late_Phase == "Wet" ~
        "Wet",
      Early_Phase == "Dry" &
        Late_Phase == "Dry" ~
        "Dry",
      TRUE ~
        "Wet-Dry transition"
    ),
    Assembly_Category = case_when(
      Assembly_Process %in%
        c(
          "Variable selection",
          "Homogeneous selection"
        ) ~
        "Deterministic",
      Assembly_Process %in%
        c(
          "Dispersal limitation",
          "Homogenizing dispersal",
          "Undominated"
        ) ~
        "Stochastic-associated",
      TRUE ~
        NA_character_
    ),
    Deterministic = if_else(
      Assembly_Category ==
        "Deterministic",
      1L,
      0L
    )
  )

saveRDS(
  assembly_consecutive,
  file.path(
    object_dir,
    "assembly_process_consecutive_transitions.rds"
  )
)

write.table(
  assembly_consecutive,
  file.path(
    table_dir,
    "assembly_process_consecutive_transitions.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================================ #
# Figure 9A | Temporal betaNTI trajectories
# ============================================================================ #

bnti_time_summary <- assembly_consecutive |>
  mutate(
    Early_Date = as.Date(
      Early_Date
    ),
    Late_Date = as.Date(
      Late_Date
    ),
    Transition_Date = Early_Date +
      round(
        as.numeric(
          Late_Date -
            Early_Date
        ) /
          2
      )
  ) |>
  group_by(
    Basin,
    Early_Time_Point,
    Late_Time_Point,
    Transition,
    Transition_Date
  ) |>
  summarise(
    N = n(),
    Mean_betaNTI = mean(
      betaNTI,
      na.rm = TRUE
    ),
    Median_betaNTI = median(
      betaNTI,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

bnti_plot_data <- bnti_time_summary |>
  mutate(
    Basin = factor(
      Basin,
      levels = c(
        "YSB",
        "HCP",
        "CHF"
      )
    )
  )

figure9a <- ggplot(
  bnti_plot_data,
  aes(
    x = Transition_Date,
    y = Mean_betaNTI
  )
) +
  geom_smooth(
    aes(
      colour = Basin,
      fill = Basin,
      group = Basin
    ),
    method = "loess",
    se = TRUE,
    span = 0.8,
    alpha = 0.15,
    linewidth = 1.5
  ) +
  geom_hline(
    yintercept = c(
      -2,
      2
    ),
    linetype = "dashed",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = 0,
    colour = "black",
    linewidth = 0.35
  ) +
  geom_point(
    aes(
      fill = Basin
    ),
    shape = 21,
    size = 2.5,
    colour = "black",
    stroke = 0.8
  ) +
  scale_colour_manual(
    values = basin_colors,
    breaks = basins,
    labels = basin_labels,
    name = "Study basin"
  ) +
  scale_fill_manual(
    values = basin_colors,
    breaks = basins,
    labels = basin_labels,
    name = "Study basin"
  ) +
  scale_x_date(
    breaks = as.Date(
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
    ),
    labels = scales::date_format(
      "%d-%b"
    ),
    expand = expansion(
      mult = 0.01
    )
  ) +
  scale_y_continuous(
    breaks = c(
      -5,
      -2,
      0,
      2,
      5
    ),
    expand = expansion(
      mult = c(
        0.02,
        0.02
      )
    )
  ) +
  labs(
    x = NULL,
    y = expression(
      beta * "NTI"
    )
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
    panel.grid.major.x = element_line(
      colour = "grey70",
      linetype = "dotted",
      linewidth = 0.4
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  coord_cartesian(
    ylim = c(
      -5,
      5.7
    ),
    clip = "on"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure9A_betaNTI_temporal_all_basins.tiff"
  ),
  figure9a,
  width = 12,
  height = 5.5,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure9A_betaNTI_temporal_all_basins.svg"
  ),
  figure9a,
  width = 12,
  height = 5.5,
  units = "in"
)


# ============================================================================ #
# Figure 9B | Assembly processes by basin and phase
# ============================================================================ #

assembly_plot_df <- assembly_consecutive |>
  filter(
    Phase_State %in%
      c(
        "Wet",
        "Dry"
      )
  ) |>
  mutate(
    Process_plot = recode(
      as.character(
        Assembly_Process
      ),
      "Undominated" =
        "Undominated or Drift"
    )
  ) |>
  count(
    Basin,
    Phase_State,
    Process_plot,
    name = "N"
  ) |>
  complete(
    Basin,
    Phase_State,
    Process_plot = process_order,
    fill = list(
      N = 0
    )
  ) |>
  group_by(
    Basin,
    Phase_State
  ) |>
  mutate(
    Proportion = N /
      sum(N)
  ) |>
  ungroup() |>
  mutate(
    Basin = factor(
      Basin,
      levels = basins
    ),
    Phase_State = factor(
      Phase_State,
      levels = c(
        "Wet",
        "Dry"
      )
    ),
    Process_plot = factor(
      Process_plot,
      levels = process_order
    )
  )

wet_x <- 1
dry_x <- 1.65
bar_halfwidth <- 0.20

assembly_connected_stack <- assembly_plot_df |>
  mutate(
    Stack_order = factor(
      as.character(
        Process_plot
      ),
      levels = rev(
        process_order
      )
    )
  ) |>
  arrange(
    Basin,
    Phase_State,
    Stack_order
  ) |>
  group_by(
    Basin,
    Phase_State
  ) |>
  mutate(
    ymax = cumsum(
      Proportion
    ),
    ymin = lag(
      ymax,
      default = 0
    )
  ) |>
  ungroup()

assembly_connected_wide <- assembly_connected_stack |>
  select(
    Basin,
    Phase_State,
    Process_plot,
    ymin,
    ymax
  ) |>
  pivot_wider(
    names_from = Phase_State,
    values_from = c(
      ymin,
      ymax
    )
  )

assembly_ribbons <- assembly_connected_wide |>
  group_by(
    Basin,
    Process_plot
  ) |>
  group_modify(
    ~ tibble(
      x = c(
        wet_x +
          bar_halfwidth,
        dry_x -
          bar_halfwidth,
        dry_x -
          bar_halfwidth,
        wet_x +
          bar_halfwidth
      ),
      y = c(
        .x$ymin_Wet,
        .x$ymin_Dry,
        .x$ymax_Dry,
        .x$ymax_Wet
      )
    )
  ) |>
  ungroup()

assembly_connected_bars <- assembly_connected_stack |>
  mutate(
    Phase_x = case_when(
      Phase_State ==
        "Wet" ~
        wet_x,
      Phase_State ==
        "Dry" ~
        dry_x
    )
  )

figure9b <- ggplot() +
  geom_polygon(
    data = assembly_ribbons,
    aes(
      x = x,
      y = y,
      group = interaction(
        Basin,
        Process_plot
      ),
      fill = Process_plot
    ),
    alpha = 0.50,
    colour = NA
  ) +
  geom_rect(
    data = assembly_connected_bars,
    aes(
      xmin = Phase_x -
        bar_halfwidth,
      xmax = Phase_x +
        bar_halfwidth,
      ymin = ymin,
      ymax = ymax,
      fill = Process_plot
    ),
    colour = "black",
    linewidth = 0.35
  ) +
  facet_wrap(
    ~ Basin,
    nrow = 1,
    labeller = as_labeller(
      basin_labels
    )
  ) +
  scale_fill_manual(
    values = assembly_colors,
    breaks = process_order,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = c(
      wet_x,
      dry_x
    ),
    labels = c(
      "Wet",
      "Dry"
    ),
    limits = c(
      wet_x - 0.35,
      dry_x + 0.35
    ),
    expand = c(
      0,
      0
    )
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    ),
    labels = percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  labs(
    x = NULL,
    y = "Relative contribution (%)",
    fill = "Assembly process"
  ) +
  theme_cowplot() +
  theme(
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure9B_assembly_processes_wet_dry.tiff"
  ),
  figure9b,
  width = 8.5,
  height = 4.5,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure9B_assembly_processes_wet_dry.svg"
  ),
  figure9b,
  width = 8.5,
  height = 4.5,
  units = "in"
)


# ============================================================================ #
# Deterministic-assembly mixed model
# ============================================================================ #

assembly_model_df <- assembly_consecutive |>
  filter(
    Phase_State %in%
      c(
        "Wet",
        "Dry"
      )
  ) |>
  mutate(
    Phase_State = factor(
      Phase_State,
      levels = c(
        "Wet",
        "Dry"
      )
    ),
    Basin = factor(
      Basin,
      levels = basins
    )
  )

model_det_phase <- lme4::glmer(
  Deterministic ~
    Basin *
    Phase_State +
    (
      1 |
        Sampling_Point
    ),
  data = assembly_model_df,
  family = binomial(
    link = "logit"
  )
)

phase_emmeans <- emmeans::emmeans(
  model_det_phase,
  ~ Phase_State |
    Basin,
  type = "response"
)

phase_contrasts <- emmeans::contrast(
  phase_emmeans,
  method = "revpairwise",
  adjust = "BH"
)

saveRDS(
  model_det_phase,
  file.path(
    object_dir,
    "model_deterministic_phase_glmer.rds"
  )
)


# ============================================================================ #
# Figure 9C | RCbray by basin and phase
# ============================================================================ #

rcbray_plot_df <- assembly_consecutive |>
  filter(
    Phase_State %in%
      c(
        "Wet",
        "Dry"
      ),
    !is.na(
      RCbray
    )
  ) |>
  mutate(
    Basin = factor(
      Basin,
      levels = basins
    ),
    Phase_State = factor(
      Phase_State,
      levels = c(
        "Wet",
        "Dry"
      )
    )
  )

rcbray_phase_stats <- lapply(
  basins,
  function(
      basin_i
  ) {

    dat_i <- rcbray_plot_df |>
      filter(
        Basin ==
          basin_i
      )

    wet_values <- dat_i |>
      filter(
        Phase_State ==
          "Wet"
      ) |>
      pull(
        RCbray
      )

    dry_values <- dat_i |>
      filter(
        Phase_State ==
          "Dry"
      ) |>
      pull(
        RCbray
      )

    test_i <- wilcox.test(
      wet_values,
      dry_values,
      paired = FALSE,
      exact = FALSE,
      conf.int = FALSE
    )

    pairwise_diff <- outer(
      wet_values,
      dry_values,
      FUN = "-"
    )

    cliffs_delta <- (
      sum(
        pairwise_diff > 0,
        na.rm = TRUE
      ) -
        sum(
          pairwise_diff < 0,
          na.rm = TRUE
        )
    ) /
      (
        length(wet_values) *
          length(dry_values)
      )

    tibble(
      Basin = basin_i,
      N_Wet = length(
        wet_values
      ),
      N_Dry = length(
        dry_values
      ),
      Median_Wet = median(
        wet_values,
        na.rm = TRUE
      ),
      IQR_Wet = IQR(
        wet_values,
        na.rm = TRUE
      ),
      Median_Dry = median(
        dry_values,
        na.rm = TRUE
      ),
      IQR_Dry = IQR(
        dry_values,
        na.rm = TRUE
      ),
      Wilcoxon_W = unname(
        test_i$statistic
      ),
      P_value = test_i$p.value,
      Cliffs_delta_Wet_vs_Dry =
        cliffs_delta
    )
  }
) |>
  bind_rows() |>
  mutate(
    P_adjust_BH = p.adjust(
      P_value,
      method = "BH"
    )
  )

figure9c <- ggplot(
  rcbray_plot_df,
  aes(
    x = Basin,
    y = RCbray
  )
) +
  geom_hline(
    yintercept = c(
      -0.95,
      0.95
    ),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35
  ) +
  gghalves::geom_half_violin(
    data = filter(
      rcbray_plot_df,
      Phase_State ==
        "Wet"
    ),
    aes(
      fill = Phase_State
    ),
    side = "l",
    trim = TRUE,
    alpha = 0.60,
    colour = "black",
    linewidth = 0.7
  ) +
  gghalves::geom_half_violin(
    data = filter(
      rcbray_plot_df,
      Phase_State ==
        "Dry"
    ),
    aes(
      fill = Phase_State
    ),
    side = "r",
    trim = TRUE,
    alpha = 0.60,
    colour = "black",
    linewidth = 0.7
  ) +
  geom_point(
    aes(
      colour = Phase_State,
      group = Phase_State
    ),
    position = position_jitterdodge(
      jitter.width = 0.06,
      jitter.height = 0,
      dodge.width = 0.22
    ),
    size = 1.5,
    alpha = 0.65
  ) +
  geom_boxplot(
    aes(
      fill = Phase_State
    ),
    width = 0.14,
    outlier.shape = NA,
    linewidth = 0.7,
    position = position_dodge(
      width = 0.18
    )
  ) +
  scale_fill_manual(
    values = phase_colors,
    name = "Phase"
  ) +
  scale_colour_manual(
    values = phase_colors,
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = c(
      -1,
      -0.95,
      -0.5,
      0,
      0.5,
      0.95,
      1
    ),
    expand = expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +
  labs(
    x = NULL,
    y = expression(
      RC[Bray]
    )
  ) +
  theme_cowplot() +
  theme(
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.5
    ),
    legend.position = "right"
  ) +
  coord_cartesian(
    ylim = c(
      -1.08,
      1.13
    ),
    clip = "off"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure9C_RCbray_wet_dry.tiff"
  ),
  figure9c,
  width = 7.5,
  height = 5.2,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure9C_RCbray_wet_dry.svg"
  ),
  figure9c,
  width = 7.5,
  height = 5.2,
  units = "in"
)


# ============================================================================ #
# Figure Supplementary 8 | Time-by-time mean betaNTI
# ============================================================================ #

timepoint_to_numeric <- function(
    x
) {

  as.integer(
    gsub(
      "^T",
      "",
      as.character(
        x
      )
    )
  )
}

time_interval_lookup <- bind_rows(
  lapply(
    basins,
    function(b) {

      md <- data.frame(
        sample_data(
          ps_list[[b]]
        ),
        check.names = FALSE
      )

      md |>
        transmute(
          Basin = b,
          Time_Point = as.integer(
            Time_Point
          ),
          Time_Interval = as.numeric(
            Time_Interval
          )
        ) |>
        distinct()
    }
  )
)

bnti_heatmap_df <- assembly_process_all |>
  mutate(
    Early_Time_Point_num =
      timepoint_to_numeric(
        Early_Time
      ),
    Late_Time_Point_num =
      timepoint_to_numeric(
        Late_Time
      )
  ) |>
  group_by(
    Basin,
    Early_Time_Point_num,
    Late_Time_Point_num
  ) |>
  summarise(
    N = n(),
    Mean_betaNTI = mean(
      betaNTI,
      na.rm = TRUE
    ),
    Median_betaNTI = median(
      betaNTI,
      na.rm = TRUE
    ),
    Min_betaNTI = min(
      betaNTI,
      na.rm = TRUE
    ),
    Max_betaNTI = max(
      betaNTI,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  left_join(
    time_interval_lookup |>
      rename(
        Early_Time_Point_num =
          Time_Point,
        Early_Time_Interval =
          Time_Interval
      ),
    by = c(
      "Basin",
      "Early_Time_Point_num"
    )
  ) |>
  left_join(
    time_interval_lookup |>
      rename(
        Late_Time_Point_num =
          Time_Point,
        Late_Time_Interval =
          Time_Interval
      ),
    by = c(
      "Basin",
      "Late_Time_Point_num"
    )
  ) |>
  mutate(
    Basin_Name = factor(
      Basin,
      levels = basins,
      labels = unname(
        basin_labels[
          basins
        ]
      )
    )
  )

time_interval_levels <- sort(
  unique(
    c(
      bnti_heatmap_df$Early_Time_Interval,
      bnti_heatmap_df$Late_Time_Interval
    )
  )
)

bnti_heatmap_df <- bnti_heatmap_df |>
  mutate(
    Early_Time_Interval_factor =
      factor(
        Early_Time_Interval,
        levels =
          time_interval_levels
      ),
    Late_Time_Interval_factor =
      factor(
        Late_Time_Interval,
        levels =
          time_interval_levels
      )
  )

figure_s8 <- ggplot(
  bnti_heatmap_df,
  aes(
    x = Early_Time_Interval_factor,
    y = Late_Time_Interval_factor,
    fill = Mean_betaNTI
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.12
  ) +
  facet_wrap(
    ~ Basin_Name,
    nrow = 1,
    scales = "free",
    space = "free_x"
  ) +
  scale_fill_gradient2(
    low = "#543005",
    mid = "#FFFFFF",
    high = "#015147",
    midpoint = 0,
    name = expression(
      mean ~ beta * "NTI"
    )
  ) +
  labs(
    x = "Earlier sampling time (days)",
    y = "Later sampling time (days)"
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 8
    ),
    axis.text.y = element_text(
      size = 8
    )
  )

ggsave(
  file.path(
    figure_dir,
    "Figure_Supplementary8_TimeTime_mean_betaNTI.tiff"
  ),
  figure_s8,
  width = 13,
  height = 5.2,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "Figure_Supplementary8_TimeTime_mean_betaNTI.svg"
  ),
  figure_s8,
  width = 13,
  height = 5.2,
  units = "in"
)


# ============================================================================ #
# Supplementary tables and statistical exports
# ============================================================================ #

assembly_phase_table <- assembly_plot_df |>
  transmute(
    Basin = as.character(
      Basin
    ),
    Basin_name = basin_labels[
      as.character(
        Basin
      )
    ],
    Phase = as.character(
      Phase_State
    ),
    Assembly_Process = as.character(
      Process_plot
    ),
    N = N,
    Proportion = Proportion,
    Percentage = Proportion *
      100
  )

openxlsx::write.xlsx(
  list(
    Deterministic_probability =
      as.data.frame(
        phase_emmeans
      ),
    Deterministic_phase_contrasts =
      as.data.frame(
        phase_contrasts
      ),
    RCbray_WetDry =
      rcbray_phase_stats,
    Assembly_process_proportions =
      assembly_phase_table
  ),
  file.path(
    table_dir,
    "Table_Supplementary9_betaNTI_RCbray_statistics.xlsx"
  ),
  overwrite = TRUE
)

openxlsx::write.xlsx(
  bnti_heatmap_df |>
    select(
      Basin,
      Basin_Name,
      Early_Time_Point_num,
      Late_Time_Point_num,
      Early_Time_Interval,
      Late_Time_Interval,
      N,
      Mean_betaNTI,
      Median_betaNTI,
      Min_betaNTI,
      Max_betaNTI
    ),
  file.path(
    table_dir,
    "Table_Supplementary10_pairwise_betaNTI.xlsx"
  ),
  overwrite = TRUE
)

write.table(
  bnti_time_summary,
  file.path(
    table_dir,
    "Figure9A_temporal_betaNTI_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  assembly_phase_table,
  file.path(
    table_dir,
    "Figure9B_assembly_processes_Wet_Dry_by_Basin.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  rcbray_phase_stats,
  file.path(
    table_dir,
    "Figure9C_RCbray_Wet_vs_Dry_statistics.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================================ #
# Final validation summary
# ============================================================================ #

validation_summary <- bind_rows(
  lapply(
    basins,
    function(b) {

      bnti_x <- betaNTI_results[[b]]$pairwise$betaNTI
      rc_x <- RCbray_results[[b]]$pairwise$RCbray

      tibble(
        Basin = b,
        betaNTI_pairs = length(
          bnti_x
        ),
        betaNTI_variable_selection =
          sum(
            bnti_x > 2,
            na.rm = TRUE
          ),
        betaNTI_homogeneous_selection =
          sum(
            bnti_x < -2,
            na.rm = TRUE
          ),
        betaNTI_null_range =
          sum(
            abs(bnti_x) <= 2,
            na.rm = TRUE
          ),
        RCbray_pairs = length(
          rc_x
        ),
        RCbray_dispersal_limitation =
          sum(
            rc_x > 0.95,
            na.rm = TRUE
          ),
        RCbray_homogenizing_dispersal =
          sum(
            rc_x < -0.95,
            na.rm = TRUE
          ),
        RCbray_undominated =
          sum(
            abs(rc_x) <= 0.95,
            na.rm = TRUE
          )
      )
    }
  )
)

openxlsx::write.xlsx(
  validation_summary,
  file.path(
    table_dir,
    "betaNTI_RCbray_validation_summary.xlsx"
  ),
  overwrite = TRUE
)

cat(
  "\nbetaNTI-RCbray community-assembly analysis completed.\n"
)

print(
  validation_summary
)

cat(
  "\nFigure 9A-C and Figure Supplementary 8 were exported.\n"
)
