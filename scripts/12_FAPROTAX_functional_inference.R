# ============================================================================ #
# 12_FAPROTAX_functional_inference.R
# ============================================================================ #
#
# Purpose:
# Reproduce the FAPROTAX workflow supporting Figure 10 and Supplementary
# Table 11.
#
# Workflow:
#
#   final non-rarefied phyloseq object
#       -> export ASV count table + taxonomy for FAPROTAX
#       -> external FAPROTAX step
#       -> import raw integer function-count table
#       -> Figure 10: seven selected functions through time
#       -> beta-binomial mixed models for eligible functions
#       -> within-basin Wet-Dry contrasts + BH correction
#
# Important:
# - FAPROTAX groups are potentially overlapping annotations. They are NOT
#   mutually exclusive components of the community.
# - Figure 10 renormalizes ONLY the seven displayed functions so that their
#   displayed contribution sums to 100% within each sample. This is a
#   visualization choice and is not used for statistical inference.
# - Statistical models use the original non-normalized FAPROTAX integer count
#   for each function and the original 16S library size for that sample.
# - Functions are tested when present in >=10% of samples and non-constant.
# - One beta-binomial mixed model is fitted per function:
#
#     cbind(Function_count, Library_size - Function_count) ~
#       Phase * Ephemeral_Basin + (1 | Sampling_Point)
#
# - Wet-Dry contrasts are estimated within each basin and BH-adjusted jointly
#   across all eligible function-by-basin contrasts.
#
# External-step note:
# - The master R script refers to an Ubuntu instruction file named
#   "ubuntu_FAPROTAX_1.2.11_Yungay2026.txt".
# - The current manuscript states FAPROTAX v1.2.12.
# - The exact external shell command is not present in the master script, so
#   this curated script does not invent one. Reconcile the version before the
#   public repository is finalized.
# - The raw output expected below was generated without internal normalization
#   (collapse_table.py -n none).
#
# Input:
# - output_phyloseq_tree/ps_mdecon_v2_FTree.rds
#
# Expected external FAPROTAX output:
# - output_FAPROTAX/external_results/faprotax_yungay_counts_raw.tsv
#
# Main outputs:
# - otu_table_faprotax.tsv
# - Figure10_FAPROTAX_selected_functions_zscore.*
# - Table_Supplementary11_FAPROTAX_beta_binomial.xlsx
# - model diagnostics and fitted-model RDS
#
# ============================================================================ #


# ============================================================================ #
# Packages
# ============================================================================ #

library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(ggplot2)
library(cowplot)
library(scales)
library(glmmTMB)
library(emmeans)
library(purrr)
library(openxlsx)


# ============================================================================ #
# Paths
# ============================================================================ #

phyloseq_path <- file.path(
  "output_phyloseq_tree",
  "ps_mdecon_v2_FTree.rds"
)

output_dir <- "output_FAPROTAX"

input_dir <- file.path(
  output_dir,
  "inputs"
)

external_results_dir <- file.path(
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
  external_results_dir,
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

faprotax_input_file <- file.path(
  input_dir,
  "otu_table_faprotax.tsv"
)

faprotax_raw_counts_file <- file.path(
  external_results_dir,
  "faprotax_yungay_counts_raw.tsv"
)


# ============================================================================ #
# Shared settings
# ============================================================================ #

Sys.setlocale(
  "LC_TIME",
  "C"
)

basin_levels <- c(
  "Herradura Clay Pan",
  "Ckoirama Halite Field",
  "Yungay Station Basin"
)

selected_functions <- c(
  "phototrophy",
  "aerobic_chemoheterotrophy",
  "fermentation",
  "nitrate_reduction",
  "respiration_of_sulfur_compounds",
  "methylotrophy",
  "manganese_oxidation"
)

function_labels <- c(
  "phototrophy" =
    "Phototrophy",
  "aerobic_chemoheterotrophy" =
    "Aerobic chemoheterotrophy",
  "fermentation" =
    "Fermentation",
  "nitrate_reduction" =
    "Nitrate reduction",
  "respiration_of_sulfur_compounds" =
    "Respiration of sulfur compounds",
  "methylotrophy" =
    "Methylotrophy",
  "manganese_oxidation" =
    "Manganese oxidation"
)

min_prevalence_fraction <- 0.10


# ============================================================================ #
# Load final non-rarefied phyloseq object
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

metadata$Sample_Time <- rownames(
  metadata
)

metadata$Ephemeral_Basin <- factor(
  as.character(
    metadata$Ephemeral_Basin
  ),
  levels = basin_levels
)

metadata$Sampling_Point <- factor(
  as.character(
    metadata$Sampling_Point
  )
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


# ============================================================================ #
# Date parser
# ============================================================================ #

parse_date_column <- function(
    x
) {

  if (inherits(x, "Date")) {
    return(
      as.Date(x)
    )
  }

  if (inherits(x, "POSIXt")) {
    return(
      as.Date(x)
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

  n_valid <- sum(
    !is.na(x_chr) &
      x_chr != ""
  )

  n_num <- sum(
    !is.na(x_num)
  )

  if (
    n_valid > 0 &&
      n_num /
        n_valid >=
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
# Export phyloseq table for external FAPROTAX
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

tax_df <- as.data.frame(
  tax_table(
    ps
  ),
  stringsAsFactors = FALSE
)

tax_df[] <- lapply(
  tax_df,
  as.character
)

taxonomy_string <- apply(
  tax_df,
  1,
  function(x) {

    x <- x[
      !is.na(x) &
        x != ""
    ]

    paste(
      x,
      collapse = "; "
    )
  }
)

faprotax_export <- data.frame(
  ASV_ID = rownames(
    otu_mat
  ),
  otu_mat,
  taxonomy = taxonomy_string[
    rownames(
      otu_mat
    )
  ],
  check.names = FALSE
)

write.table(
  faprotax_export,
  file = faprotax_input_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# ============================================================================ #
# External FAPROTAX step
# ============================================================================ #

# Run FAPROTAX externally using:
#
#   output_FAPROTAX/inputs/otu_table_faprotax.tsv
#
# The final R analysis requires the raw function-level integer count table
# generated WITHOUT FAPROTAX internal normalization:
#
#   collapse_table.py -n none
#
# The exact historical shell command is not preserved in the master R script.
# Do not substitute a guessed command in the public repository.
#
# Place the resulting table at:
#
#   output_FAPROTAX/external_results/faprotax_yungay_counts_raw.tsv


# ============================================================================ #
# Stop cleanly if external output is not yet available
# ============================================================================ #

if (!file.exists(
  faprotax_raw_counts_file
)) {

  message(
    "\nFAPROTAX input was exported successfully:\n",
    faprotax_input_file,
    "\n\nThe external raw FAPROTAX count table is not present yet.\n",
    "Expected file:\n",
    faprotax_raw_counts_file,
    "\n\nRun the external FAPROTAX step before continuing.\n"
  )

} else {


  # ========================================================================== #
  # Import raw FAPROTAX function counts
  # ========================================================================== #

  faprotax_primary <- read.delim(
    faprotax_raw_counts_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    comment.char = "",
    quote = "",
    stringsAsFactors = FALSE
  )

  names(
    faprotax_primary
  )[1] <- "Functions"

  faprotax_primary$Functions <- as.character(
    faprotax_primary$Functions
  )

  faprotax_primary[
    ,
    -1,
    drop = FALSE
  ] <- lapply(
    faprotax_primary[
      ,
      -1,
      drop = FALSE
    ],
    function(x) {
      suppressWarnings(
        as.numeric(
          x
        )
      )
    }
  )

  faprotax_count_matrix <- as.matrix(
    faprotax_primary[
      ,
      -1,
      drop = FALSE
    ]
  )

  storage.mode(
    faprotax_count_matrix
  ) <- "numeric"

  faprotax_count_check <- tibble(
    Functions =
      nrow(
        faprotax_count_matrix
      ),
    Samples =
      ncol(
        faprotax_count_matrix
      ),
    Minimum =
      min(
        faprotax_count_matrix,
        na.rm = TRUE
      ),
    Maximum =
      max(
        faprotax_count_matrix,
        na.rm = TRUE
      ),
    Non_integer_values =
      sum(
        abs(
          faprotax_count_matrix -
            round(
              faprotax_count_matrix
            )
        ) >
          1e-8,
        na.rm = TRUE
      ),
    Missing_values =
      sum(
        is.na(
          faprotax_count_matrix
        )
      ),
    Negative_values =
      sum(
        faprotax_count_matrix <
          0,
        na.rm = TRUE
      )
  )

  print(
    faprotax_count_check
  )

  if (
    faprotax_count_check$Samples !=
      nsamples(
        ps
      )
  ) {
    stop(
      "The FAPROTAX sample count does not match the phyloseq object."
    )
  }

  if (
    faprotax_count_check$Non_integer_values >
      0 ||
      faprotax_count_check$Missing_values >
      0 ||
      faprotax_count_check$Negative_values >
      0
  ) {
    stop(
      "The FAPROTAX table is not a valid non-negative integer count table."
    )
  }

  if (
    anyDuplicated(
      faprotax_primary$Functions
    ) >
      0
  ) {
    stop(
      "Duplicated function names were found in the FAPROTAX table."
    )
  }


  # ========================================================================== #
  # FAPROTAX function occurrence summary
  # ========================================================================== #

  function_occurrence_summary <- tibble(
    Total_rows_in_FAPROTAX_table =
      nrow(
        faprotax_count_matrix
      ),
    Functions_present_in_dataset =
      sum(
        rowSums(
          faprotax_count_matrix
        ) >
          0
      )
  )

  openxlsx::write.xlsx(
    function_occurrence_summary,
    file.path(
      table_dir,
      "FAPROTAX_function_occurrence_summary.xlsx"
    ),
    overwrite = TRUE
  )


  # ========================================================================== #
  # Figure 10 | Prepare seven selected functions
  # ========================================================================== #

  missing_selected_functions <- setdiff(
    selected_functions,
    faprotax_primary$Functions
  )

  if (
    length(
      missing_selected_functions
    ) >
      0
  ) {
    stop(
      "Selected Figure 10 functions missing from FAPROTAX output: ",
      paste(
        missing_selected_functions,
        collapse = ", "
      )
    )
  }

  faprotax_primary_long <- faprotax_primary |>
    pivot_longer(
      cols = -Functions,
      names_to = "Sample_Time",
      values_to = "Function_count"
    ) |>
    filter(
      Functions %in%
        selected_functions
    ) |>
    left_join(
      metadata,
      by = "Sample_Time"
    ) |>
    group_by(
      Sample_Time
    ) |>
    mutate(
      Sum_selected = sum(
        Function_count,
        na.rm = TRUE
      ),
      Relative_contribution_selected =
        ifelse(
          Sum_selected >
            0,
          100 *
            Function_count /
            Sum_selected,
          0
        )
    ) |>
    ungroup()

  if (
    anyNA(
      faprotax_primary_long$Ephemeral_Basin
    )
  ) {
    stop(
      "FAPROTAX sample names did not match the phyloseq metadata."
    )
  }

  faprotax_primary_long$Functions <- factor(
    faprotax_primary_long$Functions,
    levels = selected_functions
  )

  faprotax_primary_mean <- faprotax_primary_long |>
    group_by(
      Functions,
      Ephemeral_Basin,
      Date
    ) |>
    summarise(
      Mean_relative_contribution =
        mean(
          Relative_contribution_selected,
          na.rm = TRUE
        ),
      .groups = "drop"
    )


  # ========================================================================== #
  # Figure 10 | Basin-specific temporal z-scores
  # ========================================================================== #

  faprotax_primary_z <- faprotax_primary_mean |>
    group_by(
      Ephemeral_Basin,
      Functions
    ) |>
    mutate(
      Basin_function_mean =
        mean(
          Mean_relative_contribution,
          na.rm = TRUE
        ),
      Basin_function_sd =
        sd(
          Mean_relative_contribution,
          na.rm = TRUE
        ),
      Z_score =
        ifelse(
          !is.na(
            Basin_function_sd
          ) &&
            Basin_function_sd >
            0,
          (
            Mean_relative_contribution -
              Basin_function_mean
          ) /
            Basin_function_sd,
          0
        )
    ) |>
    ungroup()

  figure10_data <- faprotax_primary_z |>
    filter(
      Mean_relative_contribution >
        0
    )

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

  plot_date_breaks <- as.Date(
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

  figure10 <- ggplot(
    figure10_data,
    aes(
      x = Date,
      y = Functions
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
      fill = "gray",
      alpha = 0.25
    ) +
    geom_vline(
      data = surface_water_periods,
      aes(
        xintercept = xmin
      ),
      inherit.aes = FALSE,
      linetype = "dotted",
      color = "gray40",
      linewidth = 0.4
    ) +
    geom_vline(
      data = surface_water_periods,
      aes(
        xintercept = xmax
      ),
      inherit.aes = FALSE,
      linetype = "dotted",
      color = "gray40",
      linewidth = 0.4
    ) +
    geom_point(
      aes(
        size =
          Mean_relative_contribution,
        fill =
          Z_score
      ),
      shape = 21,
      color = "black",
      stroke = 0.5,
      alpha = 0.95
    ) +
    facet_wrap(
      ~ Ephemeral_Basin,
      dir = "v",
      ncol = 1
    ) +
    scale_fill_gradient2(
      low = "#8c510a",
      mid = "white",
      high = "#01665e",
      midpoint = 0,
      limits = c(
        -2.5,
        2.5
      ),
      breaks = c(
        -2,
        -1,
        0,
        1,
        2
      ),
      oob = scales::squish,
      name = "Temporal\nz-score"
    ) +
    scale_size_continuous(
      name =
        "Relative contribution\nwithin selected\nfunctions (%)",
      range = c(
        0.5,
        12
      ),
      breaks = c(
        1,
        20,
        40,
        60,
        80
      )
    ) +
    scale_y_discrete(
      labels = function_labels
    ) +
    scale_x_date(
      breaks = plot_date_breaks,
      labels = scales::label_date(
        "%d-%b"
      ),
      expand = expansion(
        mult = c(
          0.05,
          0.05
        )
      )
    ) +
    labs(
      x = "Date",
      y = NULL
    ) +
    theme_cowplot() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 8.5
      ),
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.5
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 12,
        face = "bold"
      ),
      panel.grid.major.x = element_line(
        color = "grey80",
        linetype = "dotted"
      )
    )

  ggsave(
    plot = figure10,
    filename = file.path(
      figure_dir,
      "Figure10_FAPROTAX_selected_functions_zscore.tiff"
    ),
    width = 11,
    height = 8,
    device = "tiff",
    dpi = 300,
    compression = "lzw"
  )

  ggsave(
    plot = figure10,
    filename = file.path(
      figure_dir,
      "Figure10_FAPROTAX_selected_functions_zscore.svg"
    ),
    width = 11,
    height = 8,
    device = "svg"
  )

  readr::write_csv(
    faprotax_primary_z,
    file.path(
      table_dir,
      "Figure10_FAPROTAX_zscore_source_data.csv"
    )
  )


  # ========================================================================== #
  # Raw-count table for statistical models
  # ========================================================================== #

  library_size_table <- tibble(
    Sample_Time = names(
      sample_sums(
        ps
      )
    ),
    Library_size = as.numeric(
      sample_sums(
        ps
      )
    )
  ) |>
    mutate(
      Library_size = as.integer(
        round(
          Library_size
        )
      )
    )

  faprotax_stats_long <- faprotax_primary |>
    pivot_longer(
      cols = -Functions,
      names_to = "Sample_Time",
      values_to = "Function_count"
    ) |>
    mutate(
      Function_count = as.integer(
        round(
          Function_count
        )
      )
    ) |>
    left_join(
      metadata |>
        select(
          Sample_Time,
          Phase,
          Ephemeral_Basin,
          Sampling_Point,
          Date
        ),
      by = "Sample_Time"
    ) |>
    left_join(
      library_size_table,
      by = "Sample_Time"
    ) |>
    mutate(
      Predicted_proportion =
        Function_count /
        Library_size,
      Predicted_relative_abundance_percent =
        100 *
        Predicted_proportion
    )

  if (
    anyNA(
      faprotax_stats_long$Phase
    ) ||
      anyNA(
        faprotax_stats_long$Ephemeral_Basin
      ) ||
      anyNA(
        faprotax_stats_long$Sampling_Point
      ) ||
      anyNA(
        faprotax_stats_long$Library_size
      )
  ) {
    stop(
      "Metadata/library-size join produced missing values."
    )
  }

  invalid_function_counts <- faprotax_stats_long |>
    filter(
      Function_count >
        Library_size
    )

  if (
    nrow(
      invalid_function_counts
    ) >
      0
  ) {

    readr::write_tsv(
      invalid_function_counts,
      file.path(
        table_dir,
        "ERROR_FAPROTAX_counts_exceed_library_size.tsv"
      )
    )

    stop(
      "At least one FAPROTAX function count exceeds its original library size."
    )
  }


  # ========================================================================== #
  # Objective prevalence and variability filter
  # ========================================================================== #

  n_analyzed_samples <- n_distinct(
    faprotax_stats_long$Sample_Time
  )

  min_positive_samples <- ceiling(
    min_prevalence_fraction *
      n_analyzed_samples
  )

  function_filter_table <- faprotax_stats_long |>
    group_by(
      Functions
    ) |>
    summarise(
      Positive_samples =
        sum(
          Function_count >
            0
        ),
      Prevalence =
        Positive_samples /
        n_analyzed_samples,
      Mean_predicted_relative_abundance_percent =
        mean(
          Predicted_relative_abundance_percent
        ),
      Maximum_predicted_relative_abundance_percent =
        max(
          Predicted_relative_abundance_percent
        ),
      Variable_abundance =
        n_distinct(
          Predicted_proportion
        ) >
        1,
      .groups = "drop"
    ) |>
    mutate(
      Minimum_positive_samples =
        min_positive_samples,
      Eligible_for_model =
        Positive_samples >=
        min_positive_samples &
        Variable_abundance,
      Selected_for_figure =
        Functions %in%
        selected_functions
    ) |>
    arrange(
      desc(
        Eligible_for_model
      ),
      desc(
        Prevalence
      ),
      Functions
    )

  functions_to_test <- function_filter_table |>
    filter(
      Eligible_for_model
    ) |>
    pull(
      Functions
    )

  readr::write_tsv(
    function_filter_table,
    file.path(
      table_dir,
      "FAPROTAX_function_prevalence_filter.tsv"
    )
  )

  if (
    length(
      functions_to_test
    ) ==
      0
  ) {
    stop(
      "No FAPROTAX function passed the prevalence/variability filter."
    )
  }


  # ========================================================================== #
  # Helper | Find emmeans output columns
  # ========================================================================== #

  find_emmeans_columns <- function(
      result_table
  ) {

    lower_name <- intersect(
      c(
        "lower.CL",
        "asymp.LCL"
      ),
      colnames(
        result_table
      )
    )

    upper_name <- intersect(
      c(
        "upper.CL",
        "asymp.UCL"
      ),
      colnames(
        result_table
      )
    )

    statistic_name <- intersect(
      c(
        "z.ratio",
        "t.ratio"
      ),
      colnames(
        result_table
      )
    )

    if (
      length(
        lower_name
      ) ==
        0 ||
        length(
          upper_name
        ) ==
        0
    ) {
      stop(
        "Confidence-limit columns were not found in emmeans output."
      )
    }

    list(
      lower =
        lower_name[[1]],
      upper =
        upper_name[[1]],
      statistic =
        if (
          length(
            statistic_name
          ) ==
            0
        ) {
          NA_character_
        } else {
          statistic_name[[1]]
        }
    )
  }


  # ========================================================================== #
  # Helper | Fit one beta-binomial model
  # ========================================================================== #

  fit_faprotax_function <- function(
      function_name
  ) {

    function_data <- faprotax_stats_long |>
      filter(
        Functions ==
          function_name
      ) |>
      droplevels()

    fitted_model <- tryCatch(
      glmmTMB::glmmTMB(
        cbind(
          Function_count,
          Library_size -
            Function_count
        ) ~
          Phase *
          Ephemeral_Basin +
          (
            1 |
              Sampling_Point
          ),
        family =
          glmmTMB::betabinomial(
            link = "logit"
          ),
        data =
          function_data
      ),
      error = function(e) {
        e
      }
    )

    if (
      inherits(
        fitted_model,
        "error"
      )
    ) {

      return(
        list(
          model = NULL,
          contrasts = NULL,
          phase_estimates = NULL,
          diagnostics = tibble(
            Functions =
              function_name,
            Model_status =
              "failed",
            Convergence_code =
              NA_integer_,
            Positive_definite_hessian =
              NA,
            N_observations =
              nrow(
                function_data
              ),
            AIC =
              NA_real_,
            Message =
              conditionMessage(
                fitted_model
              )
          )
        )
      )
    }

    extracted <- tryCatch(
      {

        phase_emmeans <- emmeans::emmeans(
          fitted_model,
          specs =
            ~ Phase |
              Ephemeral_Basin
        )

        # Phase order is Dry, Wet.
        # c(-1, 1) therefore represents Wet minus Dry.
        wet_dry <- emmeans::contrast(
          phase_emmeans,
          method = list(
            Wet_vs_Dry =
              c(
                -1,
                1
              )
          ),
          adjust = "none"
        )

        contrast_table <- as.data.frame(
          summary(
            wet_dry,
            infer = c(
              TRUE,
              TRUE
            ),
            level = 0.95,
            adjust = "none",
            type = "link"
          )
        )

        contrast_names <- find_emmeans_columns(
          contrast_table
        )

        contrast_table$lower_log_odds <-
          contrast_table[
            [
              contrast_names$lower
            ]
          ]

        contrast_table$upper_log_odds <-
          contrast_table[
            [
              contrast_names$upper
            ]
          ]

        if (
          is.na(
            contrast_names$statistic
          )
        ) {
          contrast_table$test_statistic <-
            NA_real_
        } else {
          contrast_table$test_statistic <-
            contrast_table[
              [
                contrast_names$statistic
              ]
            ]
        }

        contrast_results <- contrast_table |>
          transmute(
            Functions =
              function_name,
            Ephemeral_Basin,
            Contrast =
              "Wet_vs_Dry",
            Log_odds_ratio_Wet_vs_Dry =
              estimate,
            SE,
            df,
            Lower_log_odds_95CI =
              lower_log_odds,
            Upper_log_odds_95CI =
              upper_log_odds,
            Odds_ratio_Wet_vs_Dry =
              exp(
                estimate
              ),
            Lower_odds_ratio_95CI =
              exp(
                lower_log_odds
              ),
            Upper_odds_ratio_95CI =
              exp(
                upper_log_odds
              ),
            Test_statistic =
              test_statistic,
            P_value =
              p.value
          )

        phase_table <- as.data.frame(
          summary(
            phase_emmeans,
            infer = c(
              TRUE,
              FALSE
            ),
            level = 0.95,
            type = "link"
          )
        )

        phase_names <- find_emmeans_columns(
          phase_table
        )

        phase_table$lower_link <-
          phase_table[
            [
              phase_names$lower
            ]
          ]

        phase_table$upper_link <-
          phase_table[
            [
              phase_names$upper
            ]
          ]

        phase_estimates <- phase_table |>
          transmute(
            Functions =
              function_name,
            Ephemeral_Basin,
            Phase,
            Estimated_logit =
              emmean,
            SE,
            df,
            Estimated_predicted_proportion =
              plogis(
                emmean
              ),
            Lower_predicted_proportion_95CI =
              plogis(
                lower_link
              ),
            Upper_predicted_proportion_95CI =
              plogis(
                upper_link
              ),
            Estimated_predicted_relative_abundance_percent =
              100 *
              plogis(
                emmean
              ),
            Lower_predicted_relative_abundance_percent_95CI =
              100 *
              plogis(
                lower_link
              ),
            Upper_predicted_relative_abundance_percent_95CI =
              100 *
              plogis(
                upper_link
              )
          )

        diagnostics <- tibble(
          Functions =
            function_name,
          Model_status =
            if (
              fitted_model$fit$convergence ==
                0 &&
                isTRUE(
                  fitted_model$sdr$pdHess
                )
            ) {
              "converged"
            } else {
              "check_model"
            },
          Convergence_code =
            fitted_model$fit$convergence,
          Positive_definite_hessian =
            isTRUE(
              fitted_model$sdr$pdHess
            ),
          N_observations =
            stats::nobs(
              fitted_model
            ),
          AIC =
            stats::AIC(
              fitted_model
            ),
          Message =
            if (
              fitted_model$fit$convergence ==
                0 &&
                isTRUE(
                  fitted_model$sdr$pdHess
                )
            ) {
              ""
            } else {
              "Inspect convergence and Hessian before interpretation."
            }
        )

        list(
          contrasts =
            contrast_results,
          phase_estimates =
            phase_estimates,
          diagnostics =
            diagnostics
        )
      },
      error = function(e) {
        e
      }
    )

    if (
      inherits(
        extracted,
        "error"
      )
    ) {

      return(
        list(
          model =
            fitted_model,
          contrasts =
            NULL,
          phase_estimates =
            NULL,
          diagnostics =
            tibble(
              Functions =
                function_name,
              Model_status =
                "extraction_failed",
              Convergence_code =
                fitted_model$fit$convergence,
              Positive_definite_hessian =
                isTRUE(
                  fitted_model$sdr$pdHess
                ),
              N_observations =
                stats::nobs(
                  fitted_model
                ),
              AIC =
                stats::AIC(
                  fitted_model
                ),
              Message =
                conditionMessage(
                  extracted
                )
            )
        )
      )
    }

    list(
      model =
        fitted_model,
      contrasts =
        extracted$contrasts,
      phase_estimates =
        extracted$phase_estimates,
      diagnostics =
        extracted$diagnostics
    )
  }


  # ========================================================================== #
  # Fit all eligible functions
  # ========================================================================== #

  faprotax_model_results <- purrr::map(
    functions_to_test,
    fit_faprotax_function
  )

  names(
    faprotax_model_results
  ) <- functions_to_test

  faprotax_models <- purrr::map(
    faprotax_model_results,
    "model"
  )

  faprotax_model_diagnostics <- bind_rows(
    purrr::map(
      faprotax_model_results,
      "diagnostics"
    )
  )

  faprotax_phase_estimates <- bind_rows(
    purrr::map(
      faprotax_model_results,
      "phase_estimates"
    )
  )

  faprotax_wet_dry_results <- bind_rows(
    purrr::map(
      faprotax_model_results,
      "contrasts"
    )
  )

  if (
    nrow(
      faprotax_wet_dry_results
    ) ==
      0
  ) {
    stop(
      "No Wet-Dry contrasts were estimated. Inspect model diagnostics."
    )
  }

  faprotax_wet_dry_results <-
    faprotax_wet_dry_results |>
    left_join(
      faprotax_model_diagnostics |>
        select(
          Functions,
          Model_status
        ),
      by = "Functions"
    ) |>
    mutate(
      P_value_for_BH =
        if_else(
          Model_status ==
            "converged",
          P_value,
          NA_real_
        ),
      P_adjust_BH =
        p.adjust(
          P_value_for_BH,
          method = "BH"
        ),
      Direction =
        case_when(
          Model_status !=
            "converged" ~
            "Model check required",
          Log_odds_ratio_Wet_vs_Dry >
            0 ~
            "Wet",
          Log_odds_ratio_Wet_vs_Dry <
            0 ~
            "Dry",
          TRUE ~
            "No difference"
        ),
      Significant_FDR_0.05 =
        !is.na(
          P_adjust_BH
        ) &
        P_adjust_BH <
        0.05,
      Selected_for_figure =
        Functions %in%
        selected_functions
    ) |>
    arrange(
      P_adjust_BH,
      Functions,
      Ephemeral_Basin
    )

  faprotax_phase_estimates <-
    faprotax_phase_estimates |>
    left_join(
      faprotax_model_diagnostics |>
        select(
          Functions,
          Model_status
        ),
      by = "Functions"
    )


  # ========================================================================== #
  # Detailed significant contrasts
  # ========================================================================== #

  phase_estimates_wide <-
    faprotax_phase_estimates |>
    filter(
      Model_status ==
        "converged"
    ) |>
    select(
      Functions,
      Ephemeral_Basin,
      Phase,
      Estimated_predicted_relative_abundance_percent
    ) |>
    distinct() |>
    pivot_wider(
      names_from = Phase,
      values_from =
        Estimated_predicted_relative_abundance_percent,
      names_prefix =
        "Estimated_"
    )

  faprotax_significant_details <-
    faprotax_wet_dry_results |>
    filter(
      Model_status ==
        "converged",
      Significant_FDR_0.05
    ) |>
    left_join(
      phase_estimates_wide,
      by = c(
        "Functions",
        "Ephemeral_Basin"
      )
    ) |>
    mutate(
      Wet_minus_Dry_percentage_points =
        Estimated_Wet -
        Estimated_Dry
    ) |>
    select(
      Ephemeral_Basin,
      Functions,
      Estimated_Dry,
      Estimated_Wet,
      Wet_minus_Dry_percentage_points,
      Log_odds_ratio_Wet_vs_Dry,
      Odds_ratio_Wet_vs_Dry,
      Lower_odds_ratio_95CI,
      Upper_odds_ratio_95CI,
      P_value,
      P_adjust_BH,
      Direction
    ) |>
    arrange(
      Ephemeral_Basin,
      P_adjust_BH
    )


  # ========================================================================== #
  # Basin-level significance summary
  # ========================================================================== #

  significance_summary <- faprotax_wet_dry_results |>
    filter(
      Model_status ==
        "converged"
    ) |>
    group_by(
      Ephemeral_Basin
    ) |>
    summarise(
      Tested_contrasts =
        n(),
      Significant_FDR_0.05 =
        sum(
          Significant_FDR_0.05,
          na.rm = TRUE
        ),
      Higher_in_Wet =
        sum(
          Significant_FDR_0.05 &
            Direction ==
            "Wet",
          na.rm = TRUE
        ),
      Higher_in_Dry =
        sum(
          Significant_FDR_0.05 &
            Direction ==
            "Dry",
          na.rm = TRUE
        ),
      .groups = "drop"
    )


  # ========================================================================== #
  # Supplementary Table 11 and model exports
  # ========================================================================== #

  openxlsx::write.xlsx(
    list(
      Function_filter =
        function_filter_table,
      Wet_Dry_all_tested_functions =
        faprotax_wet_dry_results,
      Phase_estimated_proportions =
        faprotax_phase_estimates,
      Significant_details =
        faprotax_significant_details,
      Significance_summary =
        significance_summary,
      Model_diagnostics =
        faprotax_model_diagnostics
    ),
    file.path(
      table_dir,
      "Table_Supplementary11_FAPROTAX_beta_binomial.xlsx"
    ),
    overwrite = TRUE
  )

  saveRDS(
    faprotax_models,
    file.path(
      model_dir,
      "FAPROTAX_beta_binomial_fitted_models.rds"
    )
  )

  writeLines(
    capture.output(
      sessionInfo()
    ),
    file.path(
      model_dir,
      "FAPROTAX_sessionInfo.txt"
    )
  )


  # ========================================================================== #
  # Validation summary
  # ========================================================================== #

  validation_summary <- tibble(
    Samples_analyzed =
      n_analyzed_samples,
    Minimum_positive_samples =
      min_positive_samples,
    Functions_in_external_table =
      nrow(
        faprotax_primary
      ),
    Functions_present =
      sum(
        rowSums(
          faprotax_count_matrix
        ) >
          0
      ),
    Functions_tested =
      length(
        functions_to_test
      ),
    Expected_contrasts =
      3 *
      length(
        functions_to_test
      ),
    Observed_contrasts =
      nrow(
        faprotax_wet_dry_results
      ),
    Converged_models =
      sum(
        faprotax_model_diagnostics$Model_status ==
          "converged"
      ),
    Significant_BH_contrasts =
      sum(
        faprotax_wet_dry_results$Significant_FDR_0.05,
        na.rm = TRUE
      )
  )

  openxlsx::write.xlsx(
    validation_summary,
    file.path(
      table_dir,
      "FAPROTAX_validation_summary.xlsx"
    ),
    overwrite = TRUE
  )

  cat(
    "\nFAPROTAX analysis completed.\n"
  )

  print(
    validation_summary
  )

  cat(
    "\nFigure 10 and Supplementary Table 11 were exported.\n"
  )
}
