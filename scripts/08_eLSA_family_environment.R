# ============================================================================ #
# 08_eLSA_family_environment.R
# ============================================================================ #
#
# Purpose:
# Reproduce the family-level Extended Local Similarity Analysis (eLSA) workflow
# used for Figure 6D and Table Supplementary 6.
#
# Workflow:
#
#   final non-rarefied phyloseq object
#       -> Family-level agglomeration
#       -> basin-specific prevalence/abundance filtering
#       -> eLSA input files for HCP, CHF, and YSB
#       -> external eLSA analysis
#       -> import eLSA result files
#       -> statistical-support summaries
#       -> Figure 6D basin-specific family-environment networks
#
# IMPORTANT:
# - The eLSA computation itself was run outside R. The original master R script
#   prepares the input files and imports the resulting eLSA tables, but does
#   not contain the external command used to execute eLSA.
# - Therefore this curated script does NOT invent an external command.
# - The eLSA settings documented for the final analysis were:
#       normalization: robustZ
#       interpolation: none
#       maximum delay: 1 sampling step
#       minimum occurrence: 10
#       bootstrap iterations: 1,000
#       permutations: 10,000
# - Figure 6D displays associations with nominal permutation P < 0.05.
#   Associations with native eLSA Q < 0.05 are shown at full opacity.
#
# Input:
# - output_phyloseq_tree/ps_mdecon_v2_FTree.rds
#
# External eLSA result files expected after the external step:
# - eLSA_HCP_Family_vs_Environment_results_10k.txt
# - eLSA_CHF_Family_vs_Environment_results_10k.txt
# - eLSA_YSB_Family_vs_Environment_results_10k.txt
#
# Main outputs:
# - eLSA Family and Environment input tables for each basin
# - Figure6D_eLSA_family_environment.*
# - Figure6D_eLSA_edges_ALL.tsv
# - Figure6D_eLSA_nodes_ALL.tsv
# - Figure6D_eLSA_network_summary.tsv
# - Table_Supplementary6_eLSA_all_results.tsv
#
# ============================================================================ #


# ============================================================================ #
# Packages
# ============================================================================ #

library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(readr)
library(stringr)
library(igraph)
library(tidygraph)
library(ggraph)
library(graphlayouts)
library(ggplot2)
library(patchwork)


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

output_dir <- "output_eLSA"
input_dir <- file.path(
  output_dir,
  "inputs"
)
results_dir <- file.path(
  output_dir,
  "external_results"
)
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
  input_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  results_dir,
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
# Analysis settings
# ============================================================================ #

prev_threshold <- 10
ra_threshold <- 0.01

environment_vars <- c(
  "Wc",
  "Salinity",
  "pH",
  "ORP"
)

# Final Figure 6D criterion.
EDGE_CRITERION <- "P"
EDGE_THRESHOLD <- 0.05
ABS_LS_MIN <- 0

# Native eLSA FDR threshold used to distinguish stronger statistical support.
Q_THRESHOLD <- 0.05

# External eLSA analysis settings documented for the final analysis.
ELSA_NORMALIZATION <- "robustZ"
ELSA_MAX_DELAY <- 1
ELSA_MIN_OCCURRENCE <- 10
ELSA_BOOTSTRAPS <- 1000
ELSA_PERMUTATIONS <- 10000
ELSA_INTERPOLATION <- "none"

# Network layout.
LAYOUT_METHOD <- "stress"
NETWORK_SEED <- 2026

# Visual settings.
EDGE_COLORS <- c(
  "Negative" = "#CB181D",
  "Positive" = "#238B45"
)

EDGE_WIDTH_RANGE <- c(
  0.5,
  2.5
)

NODE_SIZE_RANGE <- c(
  1.5,
  10
)


# ============================================================================ #
# Load final non-rarefied phyloseq object
# ============================================================================ #

ps <- readRDS(
  phyloseq_path
)

metadata_df <- data.frame(
  sample_data(
    ps
  ),
  check.names = FALSE
)

if (
  !"SampleID" %in%
    colnames(
      metadata_df
    )
) {
  metadata_df$SampleID <- rownames(
    metadata_df
  )
}

required_metadata <- c(
  "SampleID",
  "Basin",
  "Sampling_Point",
  "Time",
  "Time_Point",
  "Time_Interval",
  environment_vars
)

missing_metadata_columns <- setdiff(
  required_metadata,
  colnames(
    metadata_df
  )
)

if (
  length(
    missing_metadata_columns
  ) > 0
) {
  stop(
    "Missing metadata columns required for eLSA: ",
    paste(
      missing_metadata_columns,
      collapse = ", "
    )
  )
}


# ============================================================================ #
# Aggregate ASVs at Family level
# ============================================================================ #

family_values <- as.vector(
  tax_table(
    ps
  )[
    ,
    "Family"
  ]
)

keep_family <- !is.na(
  family_values
) &
  family_values != ""

ps_family <- prune_taxa(
  keep_family,
  ps
)

ps_family <- tax_glom(
  ps_family,
  taxrank = "Family",
  NArm = TRUE
)


# ============================================================================ #
# Create unique Family labels when names are ambiguous
# ============================================================================ #

tax_family <- as.data.frame(
  tax_table(
    ps_family
  )
) |>
  tibble::rownames_to_column(
    "TaxonID"
  )

duplicated_families <- tax_family |>
  dplyr::filter(
    !is.na(
      Family
    ),
    Family != ""
  ) |>
  dplyr::count(
    Family
  ) |>
  dplyr::filter(
    n > 1
  ) |>
  dplyr::pull(
    Family
  )

tax_family <- tax_family |>
  dplyr::mutate(
    Parent_taxon = dplyr::coalesce(
      as.character(
        Order
      ),
      as.character(
        Class
      ),
      as.character(
        Phylum
      ),
      as.character(
        Kingdom
      )
    ),
    Family = dplyr::if_else(
      Family %in%
        duplicated_families,
      paste0(
        Family,
        "_",
        Parent_taxon
      ),
      as.character(
        Family
      )
    ),
    Family = make.unique(
      Family,
      sep = "_"
    )
  ) |>
  dplyr::select(
    -Parent_taxon
  )

tax_table(
  ps_family
) <- tax_table(
  as.matrix(
    tax_family |>
      tibble::column_to_rownames(
        "TaxonID"
      )
  )
)


# ============================================================================ #
# Basin-specific Family filtering
# ============================================================================ #

otu_family <- as(
  otu_table(
    ps_family
  ),
  "matrix"
)

if (
  !taxa_are_rows(
    ps_family
  )
) {
  otu_family <- t(
    otu_family
  )
}

meta_family <- data.frame(
  sample_data(
    ps_family
  ),
  check.names = FALSE
)

basin_filter_stats <- purrr::map_dfr(
  unique(
    as.character(
      meta_family$Basin
    )
  ),
  function(
      basin_i
  ) {

    samples_i <- rownames(
      meta_family
    )[
      as.character(
        meta_family$Basin
      ) ==
        basin_i
    ]

    mat_i <- otu_family[
      ,
      samples_i,
      drop = FALSE
    ]

    tibble::tibble(
      TaxonID = rownames(
        mat_i
      ),
      Basin = basin_i,
      Prevalence_pct =
        100 *
        rowSums(
          mat_i > 0
        ) /
        ncol(
          mat_i
        ),
      Relative_abundance_pct =
        100 *
        rowSums(
          mat_i
        ) /
        sum(
          mat_i
        )
    )
  }
) |>
  dplyr::left_join(
    as.data.frame(
      tax_table(
        ps_family
      )
    ) |>
      tibble::rownames_to_column(
        "TaxonID"
      ) |>
      dplyr::select(
        TaxonID,
        Kingdom,
        Phylum,
        Class,
        Order,
        Family
      ),
    by = "TaxonID"
  )

families_keep <- basin_filter_stats |>
  dplyr::filter(
    Prevalence_pct >=
      prev_threshold,
    Relative_abundance_pct >=
      ra_threshold
  ) |>
  dplyr::distinct(
    TaxonID,
    Family
  )

ps_family_filt <- prune_taxa(
  taxa_names(
    ps_family
  ) %in%
    families_keep$TaxonID,
  ps_family
)


# ============================================================================ #
# Relative abundance using original library sizes
# ============================================================================ #

otu_family_filt <- as(
  otu_table(
    ps_family_filt
  ),
  "matrix"
)

if (
  !taxa_are_rows(
    ps_family_filt
  )
) {
  otu_family_filt <- t(
    otu_family_filt
  )
}

library_sizes <- sample_sums(
  ps
)

library_sizes <- library_sizes[
  colnames(
    otu_family_filt
  )
]

family_rel <- sweep(
  otu_family_filt,
  2,
  library_sizes,
  FUN = "/"
) *
  100


# ============================================================================ #
# Long-format Family dataset
# ============================================================================ #

family_taxonomy <- as.data.frame(
  tax_table(
    ps_family_filt
  )
) |>
  tibble::rownames_to_column(
    "TaxonID"
  ) |>
  dplyr::select(
    TaxonID,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family
  )

family_long <- as.data.frame(
  family_rel
) |>
  tibble::rownames_to_column(
    "TaxonID"
  ) |>
  tidyr::pivot_longer(
    cols = -TaxonID,
    names_to = "SampleID",
    values_to = "Relative_abundance"
  ) |>
  dplyr::left_join(
    family_taxonomy,
    by = "TaxonID"
  ) |>
  dplyr::left_join(
    metadata_df |>
      dplyr::select(
        SampleID,
        Basin,
        Sampling_Point,
        Time,
        Time_Point,
        Time_Interval,
        Wc,
        Salinity,
        pH,
        ORP
      ),
    by = "SampleID"
  )

family_filter_summary <- basin_filter_stats |>
  dplyr::filter(
    Prevalence_pct >=
      prev_threshold,
    Relative_abundance_pct >=
      ra_threshold
  ) |>
  dplyr::count(
    Basin,
    name = "Families_retained"
  )

readr::write_tsv(
  basin_filter_stats,
  file.path(
    table_dir,
    "eLSA_family_filter_statistics.tsv"
  )
)

readr::write_tsv(
  family_filter_summary,
  file.path(
    table_dir,
    "eLSA_family_filter_summary.tsv"
  )
)


# ============================================================================ #
# Helper | Create eLSA input matrices for one basin
# ============================================================================ #

create_elsa_input <- function(
    basin_code,
    point_ids,
    time_offset_mode = c(
      "minus_one",
      "from_minimum"
    ),
    last_elsa_time
) {

  time_offset_mode <- match.arg(
    time_offset_mode
  )

  basin_families <- basin_filter_stats |>
    dplyr::filter(
      Basin == basin_code,
      Prevalence_pct >=
        prev_threshold,
      Relative_abundance_pct >=
        ra_threshold
    ) |>
    dplyr::pull(
      TaxonID
    )

  taxa_basin <- family_long |>
    dplyr::filter(
      Basin == basin_code,
      TaxonID %in%
        basin_families
    )

  if (
    time_offset_mode ==
      "minus_one"
  ) {

    taxa_basin <- taxa_basin |>
      dplyr::mutate(
        Time_eLSA =
          Time_Point -
          1
      )

  } else {

    min_time <- min(
      taxa_basin$Time_Point,
      na.rm = TRUE
    )

    taxa_basin <- taxa_basin |>
      dplyr::mutate(
        Time_eLSA =
          Time_Point -
          min_time
      )
  }

  taxa_basin <- taxa_basin |>
    dplyr::mutate(
      Replicate = dplyr::case_when(
        Sampling_Point ==
          point_ids[1] ~ 0,
        Sampling_Point ==
          point_ids[2] ~ 1,
        Sampling_Point ==
          point_ids[3] ~ 2,
        TRUE ~ NA_real_
      ),
      TR = paste0(
        "T",
        Time_eLSA,
        "R",
        Replicate
      )
    )

  elsa_columns <- unlist(
    lapply(
      0:last_elsa_time,
      function(
          t
      ) {
        paste0(
          "T",
          t,
          "R",
          0:2
        )
      }
    )
  )

  family_matrix <- taxa_basin |>
    dplyr::select(
      Family,
      TR,
      Relative_abundance
    ) |>
    tidyr::pivot_wider(
      names_from = TR,
      values_from = Relative_abundance
    )

  for (
    col_i in setdiff(
      elsa_columns,
      names(
        family_matrix
      )
    )
  ) {
    family_matrix[
      [
        col_i
      ]
    ] <- NA_real_
  }

  family_matrix <- family_matrix |>
    dplyr::select(
      Family,
      dplyr::all_of(
        elsa_columns
      )
    )

  names(
    family_matrix
  )[1] <- "#F"

  environment_matrix <- metadata_df |>
    dplyr::filter(
      Basin == basin_code
    )

  if (
    time_offset_mode ==
      "minus_one"
  ) {

    environment_matrix <- environment_matrix |>
      dplyr::mutate(
        Time_eLSA =
          Time_Point -
          1
      )

  } else {

    min_time <- min(
      environment_matrix$Time_Point,
      na.rm = TRUE
    )

    environment_matrix <- environment_matrix |>
      dplyr::mutate(
        Time_eLSA =
          Time_Point -
          min_time
      )
  }

  environment_matrix <- environment_matrix |>
    dplyr::mutate(
      Replicate = dplyr::case_when(
        Sampling_Point ==
          point_ids[1] ~ 0,
        Sampling_Point ==
          point_ids[2] ~ 1,
        Sampling_Point ==
          point_ids[3] ~ 2,
        TRUE ~ NA_real_
      ),
      TR = paste0(
        "T",
        Time_eLSA,
        "R",
        Replicate
      )
    ) |>
    dplyr::select(
      TR,
      dplyr::all_of(
        environment_vars
      )
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        environment_vars
      ),
      names_to = "#F",
      values_to = "Value"
    ) |>
    tidyr::pivot_wider(
      names_from = TR,
      values_from = Value
    )

  for (
    col_i in setdiff(
      elsa_columns,
      names(
        environment_matrix
      )
    )
  ) {
    environment_matrix[
      [
        col_i
      ]
    ] <- NA_real_
  }

  environment_matrix <- environment_matrix |>
    dplyr::select(
      `#F`,
      dplyr::all_of(
        elsa_columns
      )
    )

  list(
    family = family_matrix,
    environment = environment_matrix
  )
}


# ============================================================================ #
# Create eLSA input files
# ============================================================================ #

elsa_input_HCP <- create_elsa_input(
  basin_code = "HCP",
  point_ids = c(
    "HCP1",
    "HCP2",
    "HCP3"
  ),
  time_offset_mode = "minus_one",
  last_elsa_time = 19
)

elsa_input_CHF <- create_elsa_input(
  basin_code = "CHF",
  point_ids = c(
    "CHF1",
    "CHF2",
    "CHF3"
  ),
  time_offset_mode = "minus_one",
  last_elsa_time = 19
)

elsa_input_YSB <- create_elsa_input(
  basin_code = "YSB",
  point_ids = c(
    "YSB1",
    "YSB2",
    "YSB3"
  ),
  time_offset_mode = "from_minimum",
  last_elsa_time = 18
)

write.table(
  elsa_input_HCP$family,
  file.path(
    input_dir,
    "eLSA_HCP_Family.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "na"
)

write.table(
  elsa_input_HCP$environment,
  file.path(
    input_dir,
    "eLSA_HCP_Environment.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "na"
)

write.table(
  elsa_input_CHF$family,
  file.path(
    input_dir,
    "eLSA_CHF_Family.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "na"
)

write.table(
  elsa_input_CHF$environment,
  file.path(
    input_dir,
    "eLSA_CHF_Environment.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "na"
)

write.table(
  elsa_input_YSB$family,
  file.path(
    input_dir,
    "eLSA_YSB_Family.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "na"
)

write.table(
  elsa_input_YSB$environment,
  file.path(
    input_dir,
    "eLSA_YSB_Environment.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "na"
)


# ============================================================================ #
# External eLSA step
# ============================================================================ #

# The original R master script does not contain the command used to run eLSA.
# Therefore the execution command is intentionally not reconstructed here.
#
# Run eLSA externally using the Family and Environment input files created
# above and the settings documented for the final analysis:
#
#   normalization        = robustZ
#   interpolation        = none
#   maximum delay        = 1
#   minimum occurrence   = 10
#   bootstrap iterations = 1,000
#   permutations         = 10,000
#
# Place the resulting files in:
#
#   output_eLSA/external_results/
#
# with the exact names:
#
#   eLSA_HCP_Family_vs_Environment_results_10k.txt
#   eLSA_CHF_Family_vs_Environment_results_10k.txt
#   eLSA_YSB_Family_vs_Environment_results_10k.txt
#
# Missing spatial replicates remain encoded as "na"; no interpolation is
# introduced by this R script.


# ============================================================================ #
# Expected external result files
# ============================================================================ #

elsa_result_files <- c(
  HCP = file.path(
    results_dir,
    "eLSA_HCP_Family_vs_Environment_results_10k.txt"
  ),
  CHF = file.path(
    results_dir,
    "eLSA_CHF_Family_vs_Environment_results_10k.txt"
  ),
  YSB = file.path(
    results_dir,
    "eLSA_YSB_Family_vs_Environment_results_10k.txt"
  )
)

missing_result_files <- elsa_result_files[
  !file.exists(
    elsa_result_files
  )
]

if (
  length(
    missing_result_files
  ) > 0
) {

  message(
    "\neLSA input files were created successfully.\n",
    "The external eLSA result files are not yet present, so the script ",
    "will stop before Figure 6D.\n\n",
    "Expected files:\n",
    paste(
      missing_result_files,
      collapse = "\n"
    ),
    "\n"
  )

} else {


  # ========================================================================== #
  # Import eLSA results
  # ========================================================================== #

  read_elsa_results <- function(
      file,
      basin,
      permutations = 10000
  ) {

    read.delim(
      file,
      header = TRUE,
      sep = "\t",
      check.names = FALSE
    ) |>
      tibble::as_tibble() |>
      dplyr::mutate(
        Basin = basin,
        Permutations_extreme = round(
          P *
            permutations
        ),
        P_MC = (
          Permutations_extreme +
            1
        ) /
          (
            permutations +
              1
          ),
        BH_FDR = p.adjust(
          P_MC,
          method = "BH"
        ),
        abs_LS = abs(
          LS
        ),
        Sign = dplyr::case_when(
          LS > 0 ~ "Positive",
          LS < 0 ~ "Negative",
          TRUE ~ "Zero"
        )
      )
  }

  elsa_hcp_results <- read_elsa_results(
    elsa_result_files[
      "HCP"
    ],
    basin = "HCP",
    permutations = ELSA_PERMUTATIONS
  )

  elsa_chf_results <- read_elsa_results(
    elsa_result_files[
      "CHF"
    ],
    basin = "CHF",
    permutations = ELSA_PERMUTATIONS
  )

  elsa_ysb_results <- read_elsa_results(
    elsa_result_files[
      "YSB"
    ],
    basin = "YSB",
    permutations = ELSA_PERMUTATIONS
  )

  elsa_all_results <- dplyr::bind_rows(
    elsa_hcp_results,
    elsa_chf_results,
    elsa_ysb_results
  ) |>
    dplyr::mutate(
      BH_FDR_global = p.adjust(
        P_MC,
        method = "BH"
      )
    )


  # ========================================================================== #
  # Statistical diagnostics
  # ========================================================================== #

  elsa_test_counts <- elsa_all_results |>
    dplyr::count(
      Basin,
      name = "Total_tests"
    )

  elsa_significance_summary <- elsa_all_results |>
    dplyr::group_by(
      Basin
    ) |>
    dplyr::summarise(
      Total = dplyr::n(),
      P_005 = sum(
        P < 0.05,
        na.rm = TRUE
      ),
      P_MC_005 = sum(
        P_MC < 0.05,
        na.rm = TRUE
      ),
      eLSA_Q_005 = sum(
        Q < 0.05,
        na.rm = TRUE
      ),
      BH_FDR_005 = sum(
        BH_FDR < 0.05,
        na.rm = TRUE
      ),
      Global_BH_005 = sum(
        BH_FDR_global < 0.05,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  elsa_min_stats <- elsa_all_results |>
    dplyr::group_by(
      Basin
    ) |>
    dplyr::summarise(
      min_P = min(
        P,
        na.rm = TRUE
      ),
      min_P_MC = min(
        P_MC,
        na.rm = TRUE
      ),
      min_Q = min(
        Q,
        na.rm = TRUE
      ),
      min_BH = min(
        BH_FDR,
        na.rm = TRUE
      ),
      min_BH_global = min(
        BH_FDR_global,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  readr::write_tsv(
    elsa_all_results,
    file.path(
      table_dir,
      "Table_Supplementary6_eLSA_all_results.tsv"
    )
  )

  readr::write_tsv(
    elsa_test_counts,
    file.path(
      table_dir,
      "eLSA_test_counts.tsv"
    )
  )

  readr::write_tsv(
    elsa_significance_summary,
    file.path(
      table_dir,
      "eLSA_significance_summary.tsv"
    )
  )

  readr::write_tsv(
    elsa_min_stats,
    file.path(
      table_dir,
      "eLSA_min_statistics.tsv"
    )
  )


  # ========================================================================== #
  # Helper | Select Figure 6D edges
  # ========================================================================== #

  select_network_edges <- function(
      results,
      basin
  ) {

    results |>
      dplyr::filter(
        Basin == basin,
        .data[
          [
            EDGE_CRITERION
          ]
        ] <
          EDGE_THRESHOLD,
        abs_LS >=
          ABS_LS_MIN
      ) |>
      dplyr::transmute(
        Basin = Basin,
        from = X,
        to = Y,
        Family = X,
        Environment = Y,
        LS = LS,
        abs_LS = abs(
          LS
        ),
        Sign = dplyr::if_else(
          LS > 0,
          "Positive",
          "Negative"
        ),
        Delay = Delay,
        Temporal_type = dplyr::if_else(
          Delay == 0,
          "Synchronous",
          "Lagged"
        ),
        Xs = Xs,
        Ys = Ys,
        Len = Len,
        lowCI = lowCI,
        upCI = upCI,
        CI_support = dplyr::case_when(
          lowCI > 0 &
            upCI > 0 ~ "Supported",
          lowCI < 0 &
            upCI < 0 ~ "Supported",
          TRUE ~ "Crosses zero"
        ),
        P = P,
        P_MC = P_MC,
        Q = Q,
        BH_FDR = BH_FDR,
        BH_FDR_global = BH_FDR_global,
        Edge_criterion = EDGE_CRITERION,
        Edge_threshold = EDGE_THRESHOLD,
        abs_LS_threshold = ABS_LS_MIN
      ) |>
      dplyr::mutate(
        Sign = factor(
          Sign,
          levels = c(
            "Negative",
            "Positive"
          )
        ),
        Temporal_type = factor(
          Temporal_type,
          levels = c(
            "Lagged",
            "Synchronous"
          )
        )
      ) |>
      dplyr::arrange(
        Environment,
        dplyr::desc(
          abs_LS
        )
      )
  }


  # ========================================================================== #
  # Family abundance and domain metadata
  # ========================================================================== #

  family_node_metadata <- family_long |>
    dplyr::group_by(
      Basin,
      Family
    ) |>
    dplyr::summarise(
      Mean_RA_pct = mean(
        Relative_abundance,
        na.rm = TRUE
      ),
      Median_RA_pct = median(
        Relative_abundance,
        na.rm = TRUE
      ),
      Domain = {
        vals <- unique(
          stats::na.omit(
            as.character(
              Kingdom
            )
          )
        )

        if (
          length(
            vals
          ) == 0
        ) {
          NA_character_
        } else {
          vals[1]
        }
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Domain_group = dplyr::case_when(
        stringr::str_detect(
          stringr::str_to_lower(
            dplyr::coalesce(
              Domain,
              ""
            )
          ),
          "archaea"
        ) ~ "Archaea",
        stringr::str_detect(
          stringr::str_to_lower(
            dplyr::coalesce(
              Domain,
              ""
            )
          ),
          "bacteria"
        ) ~ "Bacteria",
        TRUE ~ "Other"
      )
    )


  # ========================================================================== #
  # Helper | Build basin network
  # ========================================================================== #

  build_network_data <- function(
      basin
  ) {

    edges <- select_network_edges(
      elsa_all_results,
      basin
    )

    if (
      nrow(
        edges
      ) == 0
    ) {
      return(
        list(
          basin = basin,
          edges = edges,
          nodes = tibble::tibble(),
          graph = NULL
        )
      )
    }

    nodes <- tibble::tibble(
      name = unique(
        c(
          edges$from,
          edges$to
        )
      )
    ) |>
      dplyr::mutate(
        Basin = basin,
        Node_type = dplyr::if_else(
          name %in%
            environment_vars,
          "Environment",
          "Microbial"
        )
      ) |>
      dplyr::left_join(
        family_node_metadata |>
          dplyr::filter(
            Basin ==
              basin
          ) |>
          dplyr::select(
            Family,
            Mean_RA_pct,
            Median_RA_pct,
            Domain,
            Domain_group
          ),
        by = c(
          "name" = "Family"
        )
      ) |>
      dplyr::mutate(
        Domain_group = dplyr::case_when(
          Node_type ==
            "Environment" ~
            "Environment",
          TRUE ~
            Domain_group
        ),
        Domain_group = factor(
          Domain_group,
          levels = c(
            "Bacteria",
            "Archaea",
            "Other",
            "Environment"
          )
        )
      )

    graph <- tidygraph::tbl_graph(
      nodes = nodes,
      edges = edges,
      directed = FALSE
    )

    list(
      basin = basin,
      edges = edges,
      nodes = nodes,
      graph = graph
    )
  }


  # ========================================================================== #
  # Build networks
  # ========================================================================== #

  network_data <- list(
    HCP = build_network_data(
      "HCP"
    ),
    CHF = build_network_data(
      "CHF"
    ),
    YSB = build_network_data(
      "YSB"
    )
  )

  all_network_edges <- dplyr::bind_rows(
    network_data$HCP$edges,
    network_data$CHF$edges,
    network_data$YSB$edges
  )

  all_network_nodes <- dplyr::bind_rows(
    network_data$HCP$nodes,
    network_data$CHF$nodes,
    network_data$YSB$nodes
  )

  if (
    nrow(
      all_network_edges
    ) == 0
  ) {
    stop(
      "No associations pass nominal P < 0.05."
    )
  }


  # ========================================================================== #
  # Shared visual scales
  # ========================================================================== #

  EDGE_MAX_GLOBAL <- ceiling(
    max(
      all_network_edges$abs_LS,
      na.rm = TRUE
    ) *
      10
  ) /
    10

  EDGE_LIMITS_GLOBAL <- c(
    0,
    EDGE_MAX_GLOBAL
  )

  EDGE_BREAKS_GLOBAL <- pretty(
    c(
      0,
      EDGE_MAX_GLOBAL
    ),
    n = 5
  )

  EDGE_BREAKS_GLOBAL <- EDGE_BREAKS_GLOBAL[
    EDGE_BREAKS_GLOBAL > 0 &
      EDGE_BREAKS_GLOBAL <=
      EDGE_MAX_GLOBAL
  ]

  if (
    length(
      EDGE_BREAKS_GLOBAL
    ) > 5
  ) {
    keep_positions <- unique(
      round(
        seq(
          1,
          length(
            EDGE_BREAKS_GLOBAL
          ),
          length.out = 5
        )
      )
    )

    EDGE_BREAKS_GLOBAL <- EDGE_BREAKS_GLOBAL[
      keep_positions
    ]
  }

  microbial_abundances <- all_network_nodes |>
    dplyr::filter(
      Node_type ==
        "Microbial"
    ) |>
    dplyr::pull(
      Mean_RA_pct
    )

  NODE_MAX_GLOBAL <- ceiling(
    max(
      microbial_abundances,
      na.rm = TRUE
    )
  )

  NODE_LIMITS_GLOBAL <- c(
    0,
    NODE_MAX_GLOBAL
  )

  if (
    NODE_MAX_GLOBAL <= 1
  ) {
    NODE_BREAKS_GLOBAL <- c(
      0.1,
      0.25,
      0.5,
      1
    )
  } else if (
    NODE_MAX_GLOBAL <= 5
  ) {
    NODE_BREAKS_GLOBAL <- c(
      0.1,
      0.5,
      1,
      2.5,
      5
    )
  } else if (
    NODE_MAX_GLOBAL <= 10
  ) {
    NODE_BREAKS_GLOBAL <- c(
      0.1,
      1,
      2.5,
      5,
      10
    )
  } else if (
    NODE_MAX_GLOBAL <= 20
  ) {
    NODE_BREAKS_GLOBAL <- c(
      0.1,
      1,
      5,
      10,
      20
    )
  } else if (
    NODE_MAX_GLOBAL <= 40
  ) {
    NODE_BREAKS_GLOBAL <- c(
      0.1,
      1,
      5,
      20,
      30
    )
  } else if (
    NODE_MAX_GLOBAL <= 50
  ) {
    NODE_BREAKS_GLOBAL <- c(
      0.1,
      1,
      5,
      20,
      50
    )
  } else {
    NODE_BREAKS_GLOBAL <- c(
      0.1,
      1,
      10,
      50,
      NODE_MAX_GLOBAL
    )
  }

  NODE_BREAKS_GLOBAL <- unique(
    NODE_BREAKS_GLOBAL[
      NODE_BREAKS_GLOBAL > 0 &
        NODE_BREAKS_GLOBAL <=
        NODE_MAX_GLOBAL
    ]
  )


  # ========================================================================== #
  # Network summary
  # ========================================================================== #

  network_summary_edges <- all_network_edges |>
    dplyr::group_by(
      Basin
    ) |>
    dplyr::summarise(
      Associations = dplyr::n(),
      Q_supported = sum(
        Q <
          Q_THRESHOLD,
        na.rm = TRUE
      ),
      Positive = sum(
        Sign ==
          "Positive",
        na.rm = TRUE
      ),
      Negative = sum(
        Sign ==
          "Negative",
        na.rm = TRUE
      ),
      Synchronous = sum(
        Temporal_type ==
          "Synchronous",
        na.rm = TRUE
      ),
      Lagged = sum(
        Temporal_type ==
          "Lagged",
        na.rm = TRUE
      ),
      CI_supported = sum(
        CI_support ==
          "Supported",
        na.rm = TRUE
      ),
      CI_crosses_zero = sum(
        CI_support ==
          "Crosses zero",
        na.rm = TRUE
      ),
      Min_abs_LS = min(
        abs_LS,
        na.rm = TRUE
      ),
      Median_abs_LS = median(
        abs_LS,
        na.rm = TRUE
      ),
      Max_abs_LS = max(
        abs_LS,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  network_summary_nodes <- all_network_nodes |>
    dplyr::group_by(
      Basin
    ) |>
    dplyr::summarise(
      Microbial_nodes = sum(
        Node_type ==
          "Microbial"
      ),
      Environmental_nodes = sum(
        Node_type ==
          "Environment"
      ),
      Bacteria_nodes = sum(
        Domain_group ==
          "Bacteria",
        na.rm = TRUE
      ),
      Archaea_nodes = sum(
        Domain_group ==
          "Archaea",
        na.rm = TRUE
      ),
      Other_microbe_nodes = sum(
        Domain_group ==
          "Other",
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  network_summary <- dplyr::full_join(
    network_summary_edges,
    network_summary_nodes,
    by = "Basin"
  ) |>
    dplyr::mutate(
      Edge_criterion =
        EDGE_CRITERION,
      Edge_threshold =
        EDGE_THRESHOLD,
      Q_support_threshold =
        Q_THRESHOLD
    )


  # ========================================================================== #
  # Helper | Basin title
  # ========================================================================== #

  basin_title <- function(
      basin
  ) {

    dplyr::case_when(
      basin ==
        "HCP" ~
        "Herradura Clay Pan",
      basin ==
        "CHF" ~
        "Ckoirama Halite Field",
      basin ==
        "YSB" ~
        "Yungay Station Basin",
      TRUE ~
        basin
    )
  }


  # ========================================================================== #
  # Helper | Plot standardized network
  # ========================================================================== #

  plot_network <- function(
      network_object,
      basin,
      show_legend = TRUE
  ) {

    if (
      is.null(
        network_object$graph
      ) ||
        nrow(
          network_object$edges
        ) ==
          0
    ) {
      return(
        ggplot() +
          annotate(
            "text",
            x = 0,
            y = 0,
            label = "No nominal P < 0.05 associations",
            size = 5
          ) +
          xlim(
            -1,
            1
          ) +
          ylim(
            -1,
            1
          ) +
          labs(
            title = basin_title(
              basin
            )
          ) +
          theme_void()
      )
    }

    set.seed(
      NETWORK_SEED
    )

    network_layout <- ggraph::create_layout(
      network_object$graph,
      layout = LAYOUT_METHOD
    )

    ggraph::ggraph(
      network_layout
    ) +
      ggraph::geom_edge_link(
        aes(
          edge_colour = Sign,
          edge_width = abs_LS,
          edge_linetype =
            Temporal_type,
          edge_alpha = if_else(
            Q <
              Q_THRESHOLD,
            "eLSA Q < 0.05",
            "Nominal P < 0.05"
          )
        ),
        lineend = "round"
      ) +
      ggraph::scale_edge_colour_manual(
        values = EDGE_COLORS,
        name = "Association",
        drop = FALSE
      ) +
      ggraph::scale_edge_width(
        range = EDGE_WIDTH_RANGE,
        limits = EDGE_LIMITS_GLOBAL,
        breaks = EDGE_BREAKS_GLOBAL,
        name = "|LS|"
      ) +
      ggraph::scale_edge_linetype_manual(
        values = c(
          "Lagged" =
            "dashed",
          "Synchronous" =
            "solid"
        ),
        name = "Temporal relationship",
        drop = FALSE
      ) +
      ggraph::scale_edge_alpha_manual(
        values = c(
          "Nominal P < 0.05" =
            0.55,
          "eLSA Q < 0.05" =
            1
        ),
        breaks = c(
          "eLSA Q < 0.05",
          "Nominal P < 0.05"
        ),
        name = "Statistical support"
      ) +
      ggraph::geom_node_point(
        data = function(
            x
        ) {
          x |>
            dplyr::filter(
              Node_type ==
                "Microbial"
            )
        },
        aes(
          size = Mean_RA_pct,
          shape = Domain_group
        ),
        stroke = 1,
        fill = "white",
        colour = "black"
      ) +
      ggraph::geom_node_point(
        data = function(
            x
        ) {
          x |>
            dplyr::filter(
              Node_type ==
                "Environment"
            )
        },
        shape = 22,
        size = 8,
        stroke = 1,
        fill = "grey30",
        colour = "black",
        show.legend = FALSE
      ) +
      scale_shape_manual(
        values = c(
          "Bacteria" = 21,
          "Archaea" = 24,
          "Other" = 21,
          "Environment" = 22
        ),
        name = "Node type",
        drop = FALSE
      ) +
      scale_size_continuous(
        name = "Mean relative\nabundance (%)",
        range = NODE_SIZE_RANGE,
        trans = "sqrt",
        limits = NODE_LIMITS_GLOBAL,
        breaks = NODE_BREAKS_GLOBAL
      ) +
      ggraph::geom_node_text(
        aes(
          label = name,
          fontface = ifelse(
            Node_type ==
              "Environment",
            "bold",
            "plain"
          )
        ),
        repel = TRUE,
        size = 4,
        max.overlaps = Inf
      ) +
      labs(
        title = basin_title(
          basin
        )
      ) +
      theme_void() +
      theme(
        plot.title = element_text(
          size = 14,
          face = "bold",
          hjust = 0.5
        ),
        legend.position = if (
          show_legend
        ) {
          "right"
        } else {
          "none"
        },
        legend.title = element_text(
          face = "bold"
        ),
        plot.margin = margin(
          10,
          10,
          10,
          10
        )
      )
  }


  # ========================================================================== #
  # Figure 6D
  # ========================================================================== #

  p_hcp <- plot_network(
    network_data$HCP,
    basin = "HCP",
    show_legend = FALSE
  )

  p_chf <- plot_network(
    network_data$CHF,
    basin = "CHF",
    show_legend = TRUE
  )

  p_ysb <- plot_network(
    network_data$YSB,
    basin = "YSB",
    show_legend = FALSE
  )

  figure6d <- (
    p_hcp +
      p_chf +
      p_ysb
  ) +
    patchwork::plot_layout(
      ncol = 3,
      guides = "collect"
    ) &
    theme(
      legend.position = "right"
    )

  ggsave(
    file.path(
      figure_dir,
      "Figure6D_eLSA_family_environment.tiff"
    ),
    plot = figure6d,
    width = 20,
    height = 8,
    units = "in",
    dpi = 300,
    device = "tiff",
    compression = "lzw"
  )

  ggsave(
    file.path(
      figure_dir,
      "Figure6D_eLSA_family_environment.svg"
    ),
    plot = figure6d,
    width = 20,
    height = 8,
    units = "in",
    device = "svg"
  )


  # ========================================================================== #
  # Export Figure 6D tables
  # ========================================================================== #

  readr::write_tsv(
    all_network_edges,
    file.path(
      table_dir,
      "Figure6D_eLSA_edges_ALL.tsv"
    )
  )

  readr::write_tsv(
    all_network_nodes,
    file.path(
      table_dir,
      "Figure6D_eLSA_nodes_ALL.tsv"
    )
  )

  readr::write_tsv(
    network_summary,
    file.path(
      table_dir,
      "Figure6D_eLSA_network_summary.tsv"
    )
  )


  # ========================================================================== #
  # Final summary
  # ========================================================================== #

  cat(
    "\neLSA analysis import and Figure 6D completed.\n"
  )

  cat(
    "Total family-environment tests: ",
    nrow(
      elsa_all_results
    ),
    "\n",
    sep = ""
  )

  cat(
    "Nominal P < 0.05 associations: ",
    sum(
      elsa_all_results$P <
        0.05,
      na.rm = TRUE
    ),
    "\n",
    sep = ""
  )

  cat(
    "Native eLSA Q < 0.05 associations: ",
    sum(
      elsa_all_results$Q <
        0.05,
      na.rm = TRUE
    ),
    "\n",
    sep = ""
  )
}
