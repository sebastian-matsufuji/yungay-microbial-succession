# ============================================================================ #
# 01_DADA2_taxonomy.R
# ============================================================================ #
#
# Purpose:
# Process the four 16S rRNA sequencing batches independently with DADA2,
# merge the resulting sequence tables, remove chimeras, and assign taxonomy.
#
# This script is a curated version of the sequence-processing code used for
# the Yungay study. FASTQ files are assumed to use the final SampleID names,
# e.g. CHF1_T1_1.fastq.gz / CHF1_T1_2.fastq.gz and YB1_T0_1.fastq.gz /
# YB1_T0_2.fastq.gz.
#
# IMPORTANT:
# - Metadata are NOT renamed or modified in this script.
# - The four sequencing batches are processed independently, as in the
#   original analysis.
# - DADA2 filtering/inference parameters are retained from the original
#   analysis.
# - Update only the paths in the "User paths" section before running.
#
# Main outputs:
# - seqtab_batch1.rds ... seqtab_batch4.rds
# - seqtab_all.rds
# - seqtab_nochim.rds
# - taxonomy.rds
# - batch-specific read-tracking tables and QC plots
#
# ============================================================================ #


# ============================================================================ #
# Packages
# ============================================================================ #

library(dada2)
library(ggplot2)


# ============================================================================ #
# Reproducibility
# ============================================================================ #

set.seed(2026)


# ============================================================================ #
# User paths
# ============================================================================ #

# EDIT THESE PATHS LOCALLY.
# Each directory should contain the paired FASTQ files belonging to one
# sequencing batch. The raw FASTQ files are not intended to be stored in GitHub.

batch_dirs <- c(
  batch1 = "raw1_fastq.gz",
  batch2 = "raw2_fastq.gz",
  batch3 = "raw3_fastq.gz",
  batch4 = "raw4_fastq.gz"
)

# SILVA reference files used in the original analysis.
silva_trainset <- "reference_files/silva_nr99_v138.2_toGenus_trainset.fa.gz"
silva_species  <- "reference_files/silva_v138.2_assignSpecies.fa.gz"

# Output directory.
output_dir <- "output_sequence_processing"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "filtered"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "qc"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "batch_objects"), showWarnings = FALSE, recursive = TRUE)


# ============================================================================ #
# DADA2 parameters
# ============================================================================ #

# Parameters retained from the original analysis.
trunc_len <- c(249, 230)
max_ee <- c(1, 2)
trim_left <- c(19, 21)

# Paired-read merging parameters.
max_mismatch <- 0
min_overlap <- 12


# ============================================================================ #
# Helper functions
# ============================================================================ #

get_n <- function(x) {
  sum(getUniques(x))
}


get_paired_fastq <- function(path) {

  seq_Fs <- sort(list.files(
    path,
    pattern = "_1\\.fastq\\.gz$",
    full.names = TRUE
  ))

  seq_Rs <- sort(list.files(
    path,
    pattern = "_2\\.fastq\\.gz$",
    full.names = TRUE
  ))

  if (length(seq_Fs) == 0 || length(seq_Rs) == 0) {
    stop("No paired FASTQ files were found in: ", path)
  }

  sample_F <- sub("_1\\.fastq\\.gz$", "", basename(seq_Fs))
  sample_R <- sub("_2\\.fastq\\.gz$", "", basename(seq_Rs))

  if (!setequal(sample_F, sample_R)) {
    stop(
      "Forward and reverse FASTQ SampleIDs do not match in: ",
      path,
      "\nForward-only IDs: ",
      paste(setdiff(sample_F, sample_R), collapse = ", "),
      "\nReverse-only IDs: ",
      paste(setdiff(sample_R, sample_F), collapse = ", ")
    )
  }

  sample_names <- sort(sample_F)

  seq_Fs <- seq_Fs[match(sample_names, sample_F)]
  seq_Rs <- seq_Rs[match(sample_names, sample_R)]

  stopifnot(
    identical(
      sub("_1\\.fastq\\.gz$", "", basename(seq_Fs)),
      sub("_2\\.fastq\\.gz$", "", basename(seq_Rs))
    )
  )

  list(
    forward = seq_Fs,
    reverse = seq_Rs,
    sample_names = sample_names
  )
}


process_dada_batch <- function(path, batch_number) {

  message("Processing sequencing batch ", batch_number, "...")

  fastq <- get_paired_fastq(path)

  seq_Fs <- fastq$forward
  seq_Rs <- fastq$reverse
  sample_names <- fastq$sample_names

  # -------------------------------------------------------------------------- #
  # Read-quality profiles
  # -------------------------------------------------------------------------- #

  fwd_qual <- plotQualityProfile(seq_Fs)
  rev_qual <- plotQualityProfile(seq_Rs)

  saveRDS(
    fwd_qual,
    file.path(output_dir, "qc", paste0("batch", batch_number, "_forward_quality.rds"))
  )

  saveRDS(
    rev_qual,
    file.path(output_dir, "qc", paste0("batch", batch_number, "_reverse_quality.rds"))
  )

  ggsave(
    filename = file.path(
      output_dir, "qc",
      paste0("batch", batch_number, "_forward_quality.png")
    ),
    plot = fwd_qual,
    width = 10,
    height = 10,
    dpi = 300
  )

  ggsave(
    filename = file.path(
      output_dir, "qc",
      paste0("batch", batch_number, "_reverse_quality.png")
    ),
    plot = rev_qual,
    width = 10,
    height = 10,
    dpi = 300
  )

  # -------------------------------------------------------------------------- #
  # Filtering and trimming
  # -------------------------------------------------------------------------- #

  batch_filtered_dir <- file.path(
    output_dir,
    "filtered",
    paste0("batch", batch_number)
  )

  dir.create(batch_filtered_dir, showWarnings = FALSE, recursive = TRUE)

  filt_Fs <- file.path(
    batch_filtered_dir,
    paste0(sample_names, "_1_filt.fastq")
  )

  filt_Rs <- file.path(
    batch_filtered_dir,
    paste0(sample_names, "_2_filt.fastq")
  )

  names(filt_Fs) <- sample_names
  names(filt_Rs) <- sample_names

  filter_output <- filterAndTrim(
    seq_Fs,
    filt_Fs,
    seq_Rs,
    filt_Rs,
    truncLen = trunc_len,
    maxEE = max_ee,
    trimLeft = trim_left,
    verbose = TRUE,
    multithread = FALSE
  )

  # -------------------------------------------------------------------------- #
  # Error learning
  # -------------------------------------------------------------------------- #

  errF <- learnErrors(
    filt_Fs,
    nbases = 1e8,
    multithread = FALSE
  )

  errR <- learnErrors(
    filt_Rs,
    nbases = 1e8,
    multithread = FALSE
  )

  plot_errors_F <- plotErrors(errF)
  plot_errors_R <- plotErrors(errR)

  ggsave(
    filename = file.path(
      output_dir, "qc",
      paste0("batch", batch_number, "_forward_error_model.png")
    ),
    plot = plot_errors_F,
    width = 10,
    height = 10,
    dpi = 300
  )

  ggsave(
    filename = file.path(
      output_dir, "qc",
      paste0("batch", batch_number, "_reverse_error_model.png")
    ),
    plot = plot_errors_R,
    width = 10,
    height = 10,
    dpi = 300
  )

  # -------------------------------------------------------------------------- #
  # ASV inference and paired-read merging
  # -------------------------------------------------------------------------- #

  dadaFs <- dada(
    filt_Fs,
    err = errF,
    multithread = FALSE,
    pool = TRUE
  )

  dadaRs <- dada(
    filt_Rs,
    err = errR,
    multithread = FALSE,
    pool = TRUE
  )

  mergers <- mergePairs(
    dadaFs,
    filt_Fs,
    dadaRs,
    filt_Rs,
    verbose = TRUE,
    maxMismatch = max_mismatch,
    minOverlap = min_overlap
  )

  seqtab <- makeSequenceTable(mergers)

  # -------------------------------------------------------------------------- #
  # Read tracking
  # -------------------------------------------------------------------------- #

  tracking <- data.frame(
    row.names = sample_names,
    input = filter_output[, 1],
    filtered = filter_output[, 2],
    dadaF = sapply(dadaFs, get_n),
    dadaR = sapply(dadaRs, get_n),
    merged = sapply(mergers, get_n),
    finalPercReadsKept = round(
      rowSums(seqtab) / filter_output[, 1] * 100,
      1
    )
  )

  # -------------------------------------------------------------------------- #
  # Save batch outputs
  # -------------------------------------------------------------------------- #

  saveRDS(
    seqtab,
    file.path(
      output_dir,
      "batch_objects",
      paste0("seqtab_batch", batch_number, ".rds")
    )
  )

  saveRDS(
    filter_output,
    file.path(
      output_dir,
      "batch_objects",
      paste0("filter_output_batch", batch_number, ".rds")
    )
  )

  saveRDS(
    dadaFs,
    file.path(
      output_dir,
      "batch_objects",
      paste0("dadaFs_batch", batch_number, ".rds")
    )
  )

  saveRDS(
    dadaRs,
    file.path(
      output_dir,
      "batch_objects",
      paste0("dadaRs_batch", batch_number, ".rds")
    )
  )

  saveRDS(
    mergers,
    file.path(
      output_dir,
      "batch_objects",
      paste0("mergers_batch", batch_number, ".rds")
    )
  )

  saveRDS(
    tracking,
    file.path(
      output_dir,
      "batch_objects",
      paste0("tracking_batch", batch_number, ".rds")
    )
  )

  write.table(
    tracking,
    file = file.path(
      output_dir,
      paste0("tracking_batch", batch_number, ".tsv")
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  message("Finished sequencing batch ", batch_number, ".")

  seqtab
}


# ============================================================================ #
# Process the four sequencing batches
# ============================================================================ #

seqtab_batches <- vector("list", length(batch_dirs))

for (i in seq_along(batch_dirs)) {
  seqtab_batches[[i]] <- process_dada_batch(
    path = batch_dirs[[i]],
    batch_number = i
  )
}

names(seqtab_batches) <- names(batch_dirs)


# ============================================================================ #
# Merge sequence tables
# ============================================================================ #

seqtab_all <- do.call(
  mergeSequenceTables,
  seqtab_batches
)

saveRDS(
  seqtab_all,
  file.path(output_dir, "seqtab_all.rds")
)


# ============================================================================ #
# Chimera removal
# ============================================================================ #

seqtab_nochim <- removeBimeraDenovo(
  seqtab_all,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

saveRDS(
  seqtab_nochim,
  file.path(output_dir, "seqtab_nochim.rds")
)

write.table(
  seqtab_nochim,
  file = file.path(output_dir, "seqtab_nochim.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


# ============================================================================ #
# Taxonomic assignment
# ============================================================================ #

taxonomy <- assignTaxonomy(
  seqtab_nochim,
  silva_trainset,
  multithread = TRUE
)

taxonomy <- addSpecies(
  taxonomy,
  silva_species
)

saveRDS(
  taxonomy,
  file.path(output_dir, "taxonomy.rds")
)

write.table(
  taxonomy,
  file = file.path(output_dir, "taxonomy.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


# ============================================================================ #
# Export representative ASV sequences
# ============================================================================ #

asv_sequences <- getSequences(seqtab_nochim)

writeLines(
  asv_sequences,
  con = file.path(output_dir, "ASV_sequences.txt")
)


# ============================================================================ #
# Final checks
# ============================================================================ #

cat("\nSequence-processing summary\n")
cat("---------------------------\n")
cat("Samples after chimera removal:", nrow(seqtab_nochim), "\n")
cat("ASVs after chimera removal:", ncol(seqtab_nochim), "\n")
cat("Total reads after chimera removal:", sum(seqtab_nochim), "\n")
cat("Taxonomy rows:", nrow(taxonomy), "\n")

stopifnot(
  ncol(seqtab_nochim) == nrow(taxonomy)
)

cat("\nDADA2 processing and taxonomy assignment completed.\n")
