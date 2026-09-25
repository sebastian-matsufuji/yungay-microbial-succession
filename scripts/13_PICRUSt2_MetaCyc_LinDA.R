# ============================================================================ #
# 13_PICRUSt2_MetaCyc_LinDA.R
# ============================================================================ #
#
# Purpose:
# Reproduce the PICRUSt2 + MetaCyc + LinDA workflow supporting
# Figure Supplementary 9.
#
# Workflow:
#
#   final non-rarefied phyloseq object
#       -> export ASV counts + representative sequences
#       -> external PICRUSt2 v2.6.3
#       -> unstratified MetaCyc pathway table
#       -> pooled LinDA model
#       -> basin-specific LinDA models
#       -> manual selection of 28 representative significant pathways
#       -> Figure Supplementary 9
#
# Important:
# - PICRUSt2 is run on NON-RAREFIED ASV counts and representative sequences.
# - The external pipeline was run with PICRUSt2 v2.6.3.
# - LinDA uses count-type pathway input, adaptive compositional-bias
#   correction, a 0.5 pseudocount, no winsorization, and Sampling_Point as
#   a random intercept.
# - Pooled and basin-specific BH correction families are kept separate.
# - For the pooled model, the historical code retains pathways present in
#   >=3 samples across the complete dataset.
# - For basin-specific models, the historical code applies the >=3-sample
#   prevalence filter independently within each basin.
# - The 28 pathways displayed in Figure Supplementary 9 were selected after
#   biological review from BH-significant basin-specific candidates. The
#   selection affects visualization only; all LinDA results are exported.
#
# Input:
# - output_phyloseq_tree/ps_mdecon_v2_FTree.rds
#
# External PICRUSt2 output expected:
# - output_PICRUSt2/external_results/picrust2_yungay_out/
#     pathways_out/path_abun_unstrat_descrip.tsv
#
# Main outputs:
# - raw_counts.tsv
# - seqs.fna
# - sample_metadata.tsv
# - PICRUSt2_MetaCyc_pooled_LinDA_results.csv
# - PICRUSt2_MetaCyc_basin_LinDA_all_results.csv
# - Figure_Supplementary9_MetaCyc_Wet_vs_Dry_3Basins.*
# - Figure_Supplementary9_selected_28_pathways.csv
#
# ============================================================================ #


# ============================================================================ #
# Packages
# ============================================================================ #

library(phyloseq)
library(Biostrings)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(MicrobiomeStat)
library(patchwork)
library(openxlsx)


# ============================================================================ #
# Paths
# ============================================================================ #

phyloseq_path <- file.path(
  "output_phyloseq_tree",
  "ps_mdecon_v2_FTree.rds"
)

output_dir <- "output_PICRUSt2"

input_dir <- file.path(
  output_dir,
  "inputs"
)

external_dir <- file.path(
  output_dir,
  "external_results",
  "picrust2_yungay_out"
)

figure_dir <- file.path(
  output_dir,
  "figures"
)

table_dir <- file.path(
  output_dir,
  "tables"
)

model_dir <- file.path(
  output_dir,
  "models"
)

dir.create(
  input_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  external_dir,
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
  model_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

counts_file <- file.path(
  input_dir,
  "raw_counts.tsv"
)

seqs_file <- file.path(
  input_dir,
  "seqs.fna"
)

metadata_file <- file.path(
  input_dir,
  "sample_metadata.tsv"
)

metacyc_file <- file.path(
  external_dir,
  "pathways_out",
  "path_abun_unstrat_descrip.tsv"
)


# ============================================================================ #
# Shared settings
# ============================================================================ #

set.seed(123)

basin_levels <- c(
  "Herradura Clay Pan",
  "Ckoirama Halite Field",
  "Yungay Station Basin"
)

basin_short_map <- c(
  "Herradura Clay Pan" = "HCP",
  "Ckoirama Halite Field" = "CHF",
  "Yungay Station Basin" = "YSB"
)

phase_colors <- c(
  "Wet" = "steelblue4",
  "Dry" = "#d8b365"
)


# ============================================================================ #
# Load final phyloseq object
# ============================================================================ #

ps <- readRDS(
  phyloseq_path
)

sample_data(
  ps
)$Ephemeral_Basin <- factor(
  as.character(
    sample_data(
      ps
    )$Ephemeral_Basin
  ),
  levels = basin_levels
)

sample_data(
  ps
)$Sampling_Point <- factor(
  as.character(
    sample_data(
      ps
    )$Sampling_Point
  ),
  levels = c(
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
)

sample_data(
  ps
)$Time <- factor(
  as.character(
    sample_data(
      ps
    )$Time
  ),
  levels = paste0(
    "T",
    1:20
  )
)


# ============================================================================ #
# Export PICRUSt2 count table
# ============================================================================ #

otu_mat <- as(
  otu_table(
    ps
  ),
  "matrix"
)

if (!taxa_are_rows(ps)) {
  otu_mat <- t(
    otu_mat
  )
}

picrust_counts <- data.frame(
  "#OTU ID" = rownames(
    otu_mat
  ),
  otu_mat,
  check.names = FALSE
)

write.table(
  picrust_counts,
  file = counts_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  col.names = TRUE
)


# ============================================================================ #
# Export representative ASV sequences
# ============================================================================ #

seqs <- refseq(
  ps
)

seqs <- seqs[
  taxa_names(
    ps
  )
]

names(
  seqs
) <- taxa_names(
  ps
)

Biostrings::writeXStringSet(
  x = seqs,
  filepath = seqs_file,
  format = "fasta"
)


# ============================================================================ #
# Export metadata for downstream reference
# ============================================================================ #

metadata_export <- data.frame(
  sample_data(
    ps
  ),
  check.names = FALSE
)

metadata_export <- cbind(
  SampleID_picrust2 = rownames(
    metadata_export
  ),
  metadata_export,
  check.names = FALSE
)

write.table(
  metadata_export,
  file = metadata_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# ============================================================================ #
# External PICRUSt2 step
# ============================================================================ #

# The following shell commands document the actual external workflow.
# Run them in Ubuntu/WSL, not inside R.
#
# Installation used:
#
# mamba create -n picrust2 -c conda-forge -c bioconda \
#   --override-channels picrust2=2.6.3
#
# conda activate picrust2
#
# Run PICRUSt2:
#
# picrust2_pipeline.py \
#   -s output_PICRUSt2/inputs/seqs.fna \
#   -i output_PICRUSt2/inputs/raw_counts.tsv \
#   -o output_PICRUSt2/external_results/picrust2_yungay_out \
#   -p 8 \
#   --stratified \
#   --verbose
#
# Add MetaCyc descriptions:
#
# cd output_PICRUSt2/external_results/picrust2_yungay_out
#
# add_descriptions.py \
#   -i pathways_out/path_abun_unstrat.tsv.gz \
#   -m METACYC \
#   -o pathways_out/path_abun_unstrat_descrip.tsv.gz
#
# Decompress for R:
#
# gunzip -c pathways_out/path_abun_unstrat_descrip.tsv.gz \
#   > pathways_out/path_abun_unstrat_descrip.tsv


# ============================================================================ #
# Stop cleanly if external PICRUSt2 output is absent
# ============================================================================ #

if (!file.exists(
  metacyc_file
)) {

  message(
    "\nPICRUSt2 inputs were exported successfully.\n\n",
    "Run the external PICRUSt2 workflow and place the resulting MetaCyc ",
    "table at:\n",
    metacyc_file,
    "\n"
  )

} else {


  # ========================================================================== #
  # Load metadata for LinDA
  # ========================================================================== #

  metadata <- data.frame(
    sample_data(
      ps
    ),
    check.names = FALSE
  )

  metadata$sample_name <- rownames(
    metadata
  )

  metadata$Phase <- factor(
    as.character(
      metadata$Phase
    ),
    levels = c(
      "Dry",
      "Wet"
    )
  )

  metadata$Ephemeral_Basin <- factor(
    as.character(
      metadata$Ephemeral_Basin
    ),
    levels = basin_levels
  )


  # ========================================================================== #
  # Load MetaCyc pathway predictions
  # ========================================================================== #

  p2path <- as.data.frame(
    data.table::fread(
      metacyc_file
    )
  )

  required_metacyc_columns <- c(
    "pathway",
    "description"
  )

  missing_metacyc_columns <- setdiff(
    required_metacyc_columns,
    colnames(
      p2path
    )
  )

  if (
    length(
      missing_metacyc_columns
    ) >
      0
  ) {
    stop(
      "Missing expected MetaCyc columns: ",
      paste(
        missing_metacyc_columns,
        collapse = ", "
      )
    )
  }

  path_annotation <- p2path[
    ,
    c(
      "pathway",
      "description"
    ),
    drop = FALSE
  ]

  colnames(
    path_annotation
  )[1] <- "feature"

  path_abundance <- p2path[
    ,
    !colnames(
      p2path
    ) %in%
      c(
        "pathway",
        "description"
      ),
    drop = FALSE
  ]

  rownames(
    path_abundance
  ) <- p2path[
    [
      "pathway"
    ]
  ]


  # ========================================================================== #
  # Match pathway table with phyloseq metadata
  # ========================================================================== #

  common_samples <- intersect(
    metadata$sample_name,
    colnames(
      path_abundance
    )
  )

  path_abundance <- path_abundance[
    ,
    common_samples,
    drop = FALSE
  ]

  metadata_path <- metadata[
    match(
      common_samples,
      metadata$sample_name
    ),
    ,
    drop = FALSE
  ]

  rownames(
    metadata_path
  ) <- metadata_path$sample_name

  if (
    !identical(
      colnames(
        path_abundance
      ),
      metadata_path$sample_name
    )
  ) {
    stop(
      "MetaCyc samples and metadata are not in identical order."
    )
  }

  if (
    length(
      common_samples
    ) !=
      nsamples(
        ps
      )
  ) {
    warning(
      "The MetaCyc table does not contain all phyloseq samples."
    )
  }


  # ========================================================================== #
  # Relative abundance for descriptive figure bars
  # ========================================================================== #

  path_rel_abundance <- sweep(
    as.matrix(
      path_abundance
    ),
    2,
    colSums(
      path_abundance,
      na.rm = TRUE
    ),
    "/"
  )

  path_rel_abundance[
    !is.finite(
      path_rel_abundance
    )
  ] <- 0

  path_rel_long <- path_rel_abundance |>
    as.data.frame() |>
    rownames_to_column(
      "feature"
    ) |>
    pivot_longer(
      cols = -feature,
      names_to = "sample_name",
      values_to = "Relative_abundance"
    )


  # ========================================================================== #
  # Pooled MetaCyc model
  # ========================================================================== #

  # Historical pooled analysis:
  # retain pathways detected in >=3 samples across the complete dataset.

  pooled_prevalence <- rowSums(
    path_abundance >
      0
  )

  path_pooled_filt <- path_abundance[
    pooled_prevalence >=
      3,
    ,
    drop = FALSE
  ]

  finite_pooled <- apply(
    path_pooled_filt,
    1,
    function(x) {
      all(
        is.finite(
          x
        )
      )
    }
  )

  path_pooled_filt <- path_pooled_filt[
    finite_pooled,
    ,
    drop = FALSE
  ]

  metadata_linda_pooled <- data.frame(
    Phase = factor(
      as.character(
        metadata_path$Phase
      ),
      levels = c(
        "Dry",
        "Wet"
      )
    ),
    Ephemeral_Basin = factor(
      as.character(
        metadata_path$Ephemeral_Basin
      ),
      levels = basin_levels
    ),
    Sampling_Point = factor(
      as.character(
        metadata_path$Sampling_Point
      )
    ),
    row.names = colnames(
      path_pooled_filt
    )
  )

  set.seed(
    123
  )

  linda_pooled <- MicrobiomeStat::linda(
    feature.dat = as.matrix(
      path_pooled_filt
    ),
    meta.dat = metadata_linda_pooled,
    formula =
      "~ Phase + Ephemeral_Basin + (1 | Sampling_Point)",
    feature.dat.type = "count",
    prev.filter = 0,
    mean.abund.filter = 0,
    adaptive = TRUE,
    zero.handling = "pseudo-count",
    pseudo.cnt = 0.5,
    is.winsor = FALSE,
    p.adj.method = "BH",
    alpha = 0.05,
    n.cores = 1,
    verbose = TRUE
  )

  if (
    !"PhaseWet" %in%
      names(
        linda_pooled$output
      )
  ) {
    stop(
      "PhaseWet was not found in the pooled LinDA output."
    )
  }

  pooled_results <- linda_pooled$output[
    [
      "PhaseWet"
    ]
  ] |>
    as.data.frame() |>
    rownames_to_column(
      "feature"
    ) |>
    mutate(
      CI_low =
        log2FoldChange -
        qt(
          0.975,
          df
        ) *
        lfcSE,
      CI_high =
        log2FoldChange +
        qt(
          0.975,
          df
        ) *
        lfcSE,
      Direction =
        case_when(
          is.na(
            log2FoldChange
          ) ~
            NA_character_,
          log2FoldChange >
            0 ~
            "Wet",
          TRUE ~
            "Dry"
        )
    ) |>
    left_join(
      path_annotation,
      by = "feature"
    ) |>
    mutate(
      Display_name =
        if_else(
          is.na(
            description
          ) |
            description ==
            "",
          feature,
          description
        )
    )

  write.csv(
    pooled_results,
    file.path(
      table_dir,
      "PICRUSt2_MetaCyc_pooled_LinDA_results.csv"
    ),
    row.names = FALSE
  )

  saveRDS(
    linda_pooled,
    file.path(
      model_dir,
      "PICRUSt2_MetaCyc_pooled_LinDA_model.rds"
    )
  )


  # ========================================================================== #
  # Basin-specific LinDA helper
  # ========================================================================== #

  run_basin_metacyc <- function(
      basin_name,
      basin_short
  ) {

    metadata_basin <- metadata_path |>
      filter(
        Ephemeral_Basin ==
          basin_name
      ) |>
      droplevels()

    basin_samples <- metadata_basin$sample_name

    path_basin <- path_abundance[
      ,
      basin_samples,
      drop = FALSE
    ]

    # Historical basin-specific filter:
    # retain pathways present in >=3 samples in the focal basin.
    basin_prevalence <- rowSums(
      path_basin >
        0
    )

    path_basin_filt <- path_basin[
      basin_prevalence >=
        3,
      ,
      drop = FALSE
    ]

    finite_features <- apply(
      path_basin_filt,
      1,
      function(x) {
        all(
          is.finite(
            x
          )
        )
      }
    )

    path_basin_filt <- path_basin_filt[
      finite_features,
      ,
      drop = FALSE
    ]

    metadata_linda_basin <- data.frame(
      Phase = factor(
        as.character(
          metadata_basin$Phase
        ),
        levels = c(
          "Dry",
          "Wet"
        )
      ),
      Sampling_Point = factor(
        as.character(
          metadata_basin$Sampling_Point
        )
      ),
      row.names = colnames(
        path_basin_filt
      )
    )

    set.seed(
      123
    )

    linda_basin <- MicrobiomeStat::linda(
      feature.dat = as.matrix(
        path_basin_filt
      ),
      meta.dat = metadata_linda_basin,
      formula =
        "~ Phase + (1 | Sampling_Point)",
      feature.dat.type = "count",
      prev.filter = 0,
      mean.abund.filter = 0,
      adaptive = TRUE,
      zero.handling = "pseudo-count",
      pseudo.cnt = 0.5,
      is.winsor = FALSE,
      p.adj.method = "BH",
      alpha = 0.05,
      n.cores = 1,
      verbose = TRUE
    )

    if (
      !"PhaseWet" %in%
        names(
          linda_basin$output
        )
    ) {
      stop(
        "PhaseWet was not found in LinDA output for ",
        basin_name
      )
    }

    linda_phase <- linda_basin$output[
      [
        "PhaseWet"
      ]
    ]

    basin_results <- linda_phase |>
      as.data.frame() |>
      rownames_to_column(
        "feature"
      ) |>
      mutate(
        CI_low =
          log2FoldChange -
          qt(
            0.975,
            df
          ) *
          lfcSE,
        CI_high =
          log2FoldChange +
          qt(
            0.975,
            df
          ) *
          lfcSE,
        Direction =
          case_when(
            is.na(
              log2FoldChange
            ) ~
              NA_character_,
            log2FoldChange >
              0 ~
              "Wet",
            TRUE ~
              "Dry"
          ),
        Significance =
          case_when(
            is.na(
              padj
            ) ~
              "",
            padj <
              0.0001 ~
              "****",
            padj <
              0.001 ~
              "***",
            padj <
              0.01 ~
              "**",
            padj <
              0.05 ~
              "*",
            TRUE ~
              ""
          )
      ) |>
      left_join(
        path_annotation,
        by = "feature"
      ) |>
      mutate(
        Display_name =
          if_else(
            is.na(
              description
            ) |
              description ==
              "",
            feature,
            description
          )
      )

    metadata_phase_basin <- data.frame(
      sample_name =
        as.character(
          metadata_basin$sample_name
        ),
      Phase =
        factor(
          as.character(
            metadata_basin$Phase
          ),
          levels = c(
            "Dry",
            "Wet"
          )
        ),
      stringsAsFactors = FALSE
    )

    abundance_all <- path_rel_long |>
      filter(
        sample_name %in%
          basin_samples
      ) |>
      left_join(
        metadata_phase_basin,
        by = "sample_name"
      ) |>
      group_by(
        feature,
        Phase
      ) |>
      summarise(
        Mean_abundance =
          mean(
            Relative_abundance,
            na.rm = TRUE
          ) *
          100,
        SD_abundance =
          sd(
            Relative_abundance,
            na.rm = TRUE
          ) *
          100,
        .groups = "drop"
      )

    abundance_summary <- abundance_all |>
      pivot_wider(
        names_from = Phase,
        values_from = c(
          Mean_abundance,
          SD_abundance
        ),
        names_glue =
          "{Phase}_{.value}"
      )

    results_final <- basin_results |>
      left_join(
        abundance_summary,
        by = "feature"
      ) |>
      mutate(
        Basin =
          basin_name,
        Basin_short =
          basin_short,
        Pathways_before_filter =
          nrow(
            path_basin
          ),
        Pathways_retained =
          nrow(
            path_basin_filt
          ),
        .before = 1
      )

    saveRDS(
      linda_basin,
      file.path(
        model_dir,
        paste0(
          "PICRUSt2_MetaCyc_",
          basin_short,
          "_LinDA_model.rds"
        )
      )
    )

    list(
      results =
        results_final,
      abundance =
        abundance_all,
      linda =
        linda_basin,
      n_before =
        nrow(
          path_basin
        ),
      n_retained =
        nrow(
          path_basin_filt
        )
    )
  }


  # ========================================================================== #
  # Run basin-specific models
  # ========================================================================== #

  res_HCP <- run_basin_metacyc(
    "Herradura Clay Pan",
    "HCP"
  )

  res_CHF <- run_basin_metacyc(
    "Ckoirama Halite Field",
    "CHF"
  )

  res_YSB <- run_basin_metacyc(
    "Yungay Station Basin",
    "YSB"
  )

  metacyc_basin_results <- bind_rows(
    res_HCP$results,
    res_CHF$results,
    res_YSB$results
  )

  write.csv(
    metacyc_basin_results,
    file.path(
      table_dir,
      "PICRUSt2_MetaCyc_basin_LinDA_all_results.csv"
    ),
    row.names = FALSE
  )


  # ========================================================================== #
  # Significant basin-specific pathways
  # ========================================================================== #

  metacyc_sig <- metacyc_basin_results |>
    filter(
      !is.na(
        padj
      ),
      padj <
        0.05
    ) |>
    mutate(
      abs_log2FC =
        abs(
          log2FoldChange
        ),
      Max_mean_abundance =
        pmax(
          Wet_Mean_abundance,
          Dry_Mean_abundance,
          na.rm = TRUE
        ),
      Archaeal_related_name =
        grepl(
          paste(
            c(
              "archaea",
              "archaebacteria",
              "haloarchaea",
              "archaeol",
              "archaeatidyl"
            ),
            collapse = "|"
          ),
          Display_name,
          ignore.case = TRUE
        )
    )

  write.csv(
    metacyc_sig,
    file.path(
      table_dir,
      "PICRUSt2_MetaCyc_basin_all_significant_pathways.csv"
    ),
    row.names = FALSE
  )


  # ========================================================================== #
  # Candidate screening metadata
  # ========================================================================== #

  pathway_cross_basin <- metacyc_sig |>
    group_by(
      feature
    ) |>
    summarise(
      N_basins_significant =
        n_distinct(
          Basin_short
        ),
      Significant_basins =
        paste(
          sort(
            unique(
              Basin_short
            )
          ),
          collapse = ";"
        ),
      N_significant_directions =
        n_distinct(
          Direction
        ),
      Direction_pattern =
        paste(
          paste0(
            Basin_short,
            ":",
            Direction
          ),
          collapse = ";"
        ),
      .groups = "drop"
    ) |>
    mutate(
      Specificity =
        case_when(
          N_basins_significant ==
            1 ~
            "Basin-specific",
          N_basins_significant ==
            2 ~
            "Shared in 2 basins",
          N_basins_significant ==
            3 ~
            "Shared in 3 basins",
          TRUE ~
            NA_character_
        ),
      Direction_relationship =
        case_when(
          N_basins_significant ==
            1 ~
            "Basin-specific",
          N_significant_directions >
            1 ~
            "Opposite directions among basins",
          TRUE ~
            "Same direction among significant basins"
        )
    )

  metacyc_candidates <- metacyc_sig |>
    left_join(
      pathway_cross_basin,
      by = "feature"
    ) |>
    group_by(
      feature
    ) |>
    mutate(
      Effect_rank_across_sig_basins =
        min_rank(
          desc(
            abs_log2FC
          )
        ),
      Strongest_effect_in_this_basin =
        Effect_rank_across_sig_basins ==
        1
    ) |>
    ungroup() |>
    mutate(
      Priority_group =
        case_when(
          Specificity ==
            "Basin-specific" ~
            "1 | Basin-specific",
          Direction_relationship ==
            "Opposite directions among basins" ~
            "2 | Opposite basin response",
          Strongest_effect_in_this_basin ~
            "3 | Shared, strongest effect here",
          TRUE ~
            "4 | Other shared significant"
        ),
      Priority_order =
        case_when(
          Priority_group ==
            "1 | Basin-specific" ~
            1L,
          Priority_group ==
            "2 | Opposite basin response" ~
            2L,
          Priority_group ==
            "3 | Shared, strongest effect here" ~
            3L,
          TRUE ~
            4L
        )
    ) |>
    arrange(
      Basin_short,
      Priority_order,
      padj,
      desc(
        abs_log2FC
      )
    )

  write.csv(
    metacyc_candidates,
    file.path(
      table_dir,
      "PICRUSt2_MetaCyc_basin_significant_candidates.csv"
    ),
    row.names = FALSE
  )


  # ========================================================================== #
  # Final manual pathway selection for Supplementary Figure 9
  # ========================================================================== #

  # This is the historical final curated display:
  #
  # HCP = 8 pathways  | 4 Wet + 4 Dry
  # CHF = 10 pathways | 9 Wet + 1 Dry
  # YSB = 10 pathways | 4 Wet + 6 Dry
  #
  # Selection affects visualization only.

  pathway_selection <- tribble(
    ~Basin_short, ~feature, ~Selection_order,

    # HCP | Wet
    "HCP", "PWY-7294", 1,
    "HCP", "PWY-5130", 2,
    "HCP", "PWY-8011", 3,
    "HCP", "PWY-1622", 4,

    # HCP | Dry
    "HCP", "PWY0-42", 5,
    "HCP", "POLYAMINSYN3-PWY", 6,
    "HCP", "PWY-7783", 7,
    "HCP", "PWY-6760", 8,

    # CHF | Wet
    "CHF", "GALLATE-DEGRADATION-II-PWY", 1,
    "CHF", "PWY3O-4107", 2,
    "CHF", "PWY66-399", 3,
    "CHF", "PWY0-1337", 4,
    "CHF", "NONOXIPENT-PWY", 5,
    "CHF", "PWY-801", 6,
    "CHF", "HEME-BIOSYNTHESIS-II", 7,
    "CHF", "PPGPPMET-PWY", 8,
    "CHF", "GLYCOCAT-PWY", 9,

    # CHF | Dry
    "CHF", "P101-PWY", 10,

    # YSB | Wet
    "YSB", "PWY0-321", 1,
    "YSB", "P185-PWY", 2,
    "YSB", "PWY-7052", 3,
    "YSB", "PWY0-1586", 4,

    # YSB | Dry
    "YSB", "PWY-6174", 5,
    "YSB", "PWY-6165", 6,
    "YSB", "PWY490-3", 7,
    "YSB", "P23-PWY", 8,
    "YSB", "PWY-7094", 9,
    "YSB", "PWY-822", 10
  )


  # ========================================================================== #
  # Figure helper
  # ========================================================================== #

  build_basin_metacyc_plot <- function(
      basin_name,
      basin_short,
      basin_results,
      abundance_data
  ) {

    selected_basin <- pathway_selection |>
      filter(
        Basin_short ==
          basin_short
      )

    selected_results <- basin_results |>
      inner_join(
        selected_basin,
        by = c(
          "Basin_short",
          "feature"
        )
      )

    missing_selected <- setdiff(
      selected_basin$feature,
      selected_results$feature
    )

    if (
      length(
        missing_selected
      ) >
        0
    ) {
      stop(
        "Selected pathways missing for ",
        basin_short,
        ": ",
        paste(
          missing_selected,
          collapse = ", "
        )
      )
    }

    if (
      any(
        is.na(
          selected_results$padj
        ) |
          selected_results$padj >=
          0.05
      )
    ) {
      warning(
        "At least one selected pathway is not BH-significant in ",
        basin_short,
        "."
      )
    }

    metabolism_plot <- selected_results |>
      mutate(
        Direction = factor(
          Direction,
          levels = c(
            "Wet",
            "Dry"
          )
        )
      ) |>
      arrange(
        Direction,
        desc(
          abs(
            log2FoldChange
          )
        )
      ) |>
      mutate(
        plot_log2FC =
          -log2FoldChange,
        plot_CI_low =
          -CI_high,
        plot_CI_high =
          -CI_low,
        row_id =
          rev(
            seq_len(
              n()
            )
          )
      )

    abundance_plot_data <- abundance_data |>
      filter(
        feature %in%
          metabolism_plot$feature
      ) |>
      inner_join(
        metabolism_plot |>
          select(
            feature,
            row_id
          ),
        by = "feature"
      ) |>
      mutate(
        Phase = factor(
          Phase,
          levels = c(
            "Wet",
            "Dry"
          )
        ),
        y_plot =
          if_else(
            Phase ==
              "Wet",
            row_id +
              0.16,
            row_id -
              0.16
          )
      )

    p_abundance <- ggplot() +
      geom_rect(
        data = metabolism_plot |>
          filter(
            row_id %%
              2 ==
              0
          ),
        aes(
          xmin = -Inf,
          xmax = Inf,
          ymin = row_id -
            0.5,
          ymax = row_id +
            0.5
        ),
        fill = "grey94"
      ) +
      geom_col(
        data = abundance_plot_data,
        aes(
          x =
            Mean_abundance,
          y =
            y_plot,
          fill =
            Phase
        ),
        orientation = "y",
        width = 0.32,
        color = "black",
        linewidth = 0.35
      ) +
      scale_fill_manual(
        values =
          phase_colors,
        breaks = c(
          "Wet",
          "Dry"
        )
      ) +
      scale_y_continuous(
        breaks =
          metabolism_plot$row_id,
        labels =
          metabolism_plot$Display_name,
        limits = c(
          0.5,
          nrow(
            metabolism_plot
          ) +
            0.5
        ),
        expand = c(
          0,
          0
        )
      ) +
      scale_x_continuous(
        expand = expansion(
          mult = c(
            0,
            0.04
          )
        )
      ) +
      labs(
        x =
          "Mean relative abundance (%)",
        y =
          NULL,
        fill =
          basin_name
      ) +
      guides(
        fill = guide_legend(
          direction =
            "horizontal",
          title.position =
            "right",
          title.hjust =
            0
        )
      ) +
      theme_classic() +
      theme(
        axis.text.y =
          element_text(
            size = 12,
            color = "black"
          ),
        axis.text.x =
          element_text(
            size = 10,
            color = "black"
          ),
        axis.title.x =
          element_text(
            size = 11,
            color = "black"
          ),
        axis.ticks.y =
          element_blank(),
        legend.position =
          "top",
        legend.justification =
          "left",
        legend.direction =
          "horizontal",
        legend.title =
          element_text(
            face = "bold"
          )
      )

    p_effect <- ggplot() +
      geom_rect(
        data = metabolism_plot |>
          filter(
            row_id %%
              2 ==
              0
          ),
        aes(
          xmin = -Inf,
          xmax = Inf,
          ymin = row_id -
            0.5,
          ymax = row_id +
            0.5
        ),
        fill = "grey94"
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.5
      ) +
      geom_errorbarh(
        data = metabolism_plot,
        aes(
          xmin =
            plot_CI_low,
          xmax =
            plot_CI_high,
          y =
            row_id
        ),
        height = 0.20,
        linewidth = 0.60,
        color = "black"
      ) +
      geom_point(
        data = metabolism_plot,
        aes(
          x =
            plot_log2FC,
          y =
            row_id,
          fill =
            Direction
        ),
        shape = 21,
        size = 3.5,
        stroke = 0.60,
        color = "black"
      ) +
      scale_fill_manual(
        values =
          phase_colors
      ) +
      scale_y_continuous(
        limits = c(
          0.5,
          nrow(
            metabolism_plot
          ) +
            0.5
        ),
        expand = c(
          0,
          0
        )
      ) +
      scale_x_continuous(
        expand = expansion(
          mult = c(
            0.05,
            0.07
          )
        )
      ) +
      labs(
        title =
          "95% confidence intervals",
        x =
          "LinDA log2 fold change",
        y =
          NULL
      ) +
      guides(
        fill = "none"
      ) +
      theme_classic() +
      theme(
        plot.title =
          element_text(
            hjust = 0.5,
            size = 11
          ),
        axis.text.y =
          element_blank(),
        axis.ticks.y =
          element_blank(),
        axis.line.y =
          element_blank()
      )

    p_fdr <- ggplot(
      metabolism_plot,
      aes(
        x = 1,
        y = row_id,
        label = Significance
      )
    ) +
      geom_rect(
        data = metabolism_plot |>
          filter(
            row_id %%
              2 ==
              0
          ),
        aes(
          xmin = 0,
          xmax = 2,
          ymin = row_id -
            0.5,
          ymax = row_id +
            0.5
        ),
        inherit.aes = FALSE,
        fill = "grey94"
      ) +
      geom_text(
        fontface = "bold",
        size = 4.2
      ) +
      scale_x_continuous(
        limits = c(
          0,
          2
        ),
        expand = c(
          0,
          0
        )
      ) +
      scale_y_continuous(
        limits = c(
          0.5,
          nrow(
            metabolism_plot
          ) +
            0.5
        ),
        expand = c(
          0,
          0
        )
      ) +
      labs(
        title = "FDR"
      ) +
      theme_void() +
      theme(
        plot.title =
          element_text(
            hjust = 0.5,
            face = "bold",
            size = 11
          )
      )

    p_basin <-
      p_abundance +
      p_effect +
      p_fdr +
      patchwork::plot_layout(
        widths = c(
          5.0,
          2.3,
          0.50
        ),
        guides = "keep"
      )

    list(
      plot =
        p_basin,
      selected =
        metabolism_plot
    )
  }


  # ========================================================================== #
  # Build final three-basin figure
  # ========================================================================== #

  plot_HCP <- build_basin_metacyc_plot(
    basin_name =
      "Herradura Clay Pan",
    basin_short =
      "HCP",
    basin_results =
      res_HCP$results,
    abundance_data =
      res_HCP$abundance
  )

  plot_CHF <- build_basin_metacyc_plot(
    basin_name =
      "Ckoirama Halite Field",
    basin_short =
      "CHF",
    basin_results =
      res_CHF$results,
    abundance_data =
      res_CHF$abundance
  )

  plot_YSB <- build_basin_metacyc_plot(
    basin_name =
      "Yungay Station Basin",
    basin_short =
      "YSB",
    basin_results =
      res_YSB$results,
    abundance_data =
      res_YSB$abundance
  )

  figure_s9 <-
    patchwork::wrap_plots(
      list(
        HCP =
          plot_HCP$plot,
        CHF =
          plot_CHF$plot,
        YSB =
          plot_YSB$plot
      ),
      ncol = 1,
      guides = "keep"
    )

  selected_28 <- bind_rows(
    plot_HCP$selected,
    plot_CHF$selected,
    plot_YSB$selected
  )

  write.csv(
    selected_28,
    file.path(
      table_dir,
      "Figure_Supplementary9_selected_28_pathways.csv"
    ),
    row.names = FALSE
  )

  ggsave(
    plot = figure_s9,
    filename = file.path(
      figure_dir,
      "Figure_Supplementary9_MetaCyc_Wet_vs_Dry_3Basins.tiff"
    ),
    width = 16,
    height = 19.5,
    units = "in",
    dpi = 300,
    device = "tiff",
    compression = "lzw"
  )

  ggsave(
    plot = figure_s9,
    filename = file.path(
      figure_dir,
      "Figure_Supplementary9_MetaCyc_Wet_vs_Dry_3Basins.svg"
    ),
    width = 16,
    height = 19.5,
    units = "in"
  )


  # ========================================================================== #
  # Validation summary
  # ========================================================================== #

  validation_summary <- tibble(
    Total_MetaCyc_pathways =
      nrow(
        path_abundance
      ),
    Pooled_pathways_retained =
      nrow(
        path_pooled_filt
      ),
    Pooled_BH_significant =
      sum(
        pooled_results$padj <
          0.05,
        na.rm = TRUE
      ),
    HCP_pathways_retained =
      res_HCP$n_retained,
    CHF_pathways_retained =
      res_CHF$n_retained,
    YSB_pathways_retained =
      res_YSB$n_retained,
    HCP_BH_significant =
      sum(
        res_HCP$results$padj <
          0.05,
        na.rm = TRUE
      ),
    CHF_BH_significant =
      sum(
        res_CHF$results$padj <
          0.05,
        na.rm = TRUE
      ),
    YSB_BH_significant =
      sum(
        res_YSB$results$padj <
          0.05,
        na.rm = TRUE
      ),
    Selected_for_figure =
      nrow(
        selected_28
      )
  )

  openxlsx::write.xlsx(
    validation_summary,
    file.path(
      table_dir,
      "PICRUSt2_MetaCyc_validation_summary.xlsx"
    ),
    overwrite = TRUE
  )

  cat(
    "\nPICRUSt2 MetaCyc + LinDA analysis completed.\n"
  )

  print(
    validation_summary
  )

  cat(
    "\nFigure Supplementary 9 was exported.\n"
  )
}
