# ============================================================================ #
# 02_phyloseq_filter_microdecon.R
# ============================================================================ #
#
# Purpose:
# Create the initial phyloseq object from the DADA2 outputs, apply the
# taxonomic filters used in the study, and run microDecon separately for the
# five processing groups.
#
# This script reproduces the confirmed portion of the original workflow up to
# the merged microDecon table.
#
# IMPORTANT:
# - Metadata are read using the FINAL SampleID names and are not renamed.
# - Singleton ASVs are NOT removed, matching the original analysis.
# - The "all samples + five controls" microDecon block from the original script
#   was explicitly labelled "test only" and is not included here.
# - The original workflow next used a manually prepared file named
#   decont_asvtab_taxa_filt.xlsx. That step is NOT reconstructed here because
#   the original script does not document how that file was produced.
#
# Inputs:
# - output_sequence_processing/seqtab_nochim.rds
# - output_sequence_processing/taxonomy.rds
# - study metadata workbook with final SampleID names
#
# Outputs:
# - ps_raw.rds
# - ps_filt.rds
# - asvtab.rds
# - taxonomy_filtered.rds
# - decontaminated_G1.rds ... decontaminated_G5.rds
# - decont_asvtab_taxa_G1.rds ... decont_asvtab_taxa_G5.rds
# - decont_asvtab_taxa.txt
#
# ============================================================================ #


# ============================================================================ #
# Packages
# ============================================================================ #

library(phyloseq)
library(Biostrings)
library(openxlsx)
library(janitor)
library(microDecon)
library(dplyr)
library(tidyr)
library(purrr)


# ============================================================================ #
# User paths
# ============================================================================ #

dada2_dir <- "output_sequence_processing"

# EDIT THIS PATH LOCALLY.
metadata_path <- "metadata.xlsx"
metadata_sheet <- "Sheet1"

output_dir <- "output_phyloseq_microdecon"
microdecon_dir <- file.path(output_dir, "microdecon")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(microdecon_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================ #
# Load DADA2 outputs
# ============================================================================ #

seqtab_nochim <- readRDS(
  file.path(dada2_dir, "seqtab_nochim.rds")
)

taxonomy <- readRDS(
  file.path(dada2_dir, "taxonomy.rds")
)

stopifnot(
  ncol(seqtab_nochim) == nrow(taxonomy)
)


# ============================================================================ #
# Assign compact ASV identifiers
# ============================================================================ #

# The original DADA2 table used nucleotide sequences as taxon names.
# Compact ASV identifiers were assigned before constructing phyloseq.

asv_sequences <- colnames(seqtab_nochim)

ndigits <- nchar(length(asv_sequences))
asv_ids <- sprintf(
  paste0("ASV_%0", ndigits, "d"),
  seq_along(asv_sequences)
)

colnames(seqtab_nochim) <- asv_ids
rownames(taxonomy) <- asv_ids
names(asv_sequences) <- asv_ids

asv_seq <- Biostrings::DNAStringSet(asv_sequences)
names(asv_seq) <- asv_ids

Biostrings::writeXStringSet(
  asv_seq,
  filepath = file.path(output_dir, "asv_seq.fasta"),
  format = "fasta"
)


# ============================================================================ #
# Load final metadata
# ============================================================================ #

metadata <- openxlsx::read.xlsx(
  metadata_path,
  sheet = metadata_sheet,
  rowNames = FALSE
)

if (!"SampleID" %in% colnames(metadata)) {
  stop("Metadata must contain a column named 'SampleID'.")
}

if (anyDuplicated(metadata$SampleID)) {
  stop("Duplicated SampleID values were found in metadata.")
}

if (any(is.na(metadata$SampleID))) {
  stop("Missing SampleID values were found in metadata.")
}

# Convert Excel serial dates when needed.
if ("Date" %in% colnames(metadata) && is.numeric(metadata$Date)) {
  metadata$Date <- janitor::excel_numeric_to_date(
    metadata$Date,
    date_system = "modern"
  )
}

# The metadata may contain samples without usable sequence data.
# Every sequenced sample, however, must occur exactly once in metadata.
missing_metadata <- setdiff(
  rownames(seqtab_nochim),
  metadata$SampleID
)

if (length(missing_metadata) > 0) {
  stop(
    "These sequenced samples are absent from metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}

metadata <- metadata[
  match(rownames(seqtab_nochim), metadata$SampleID),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    metadata$SampleID,
    rownames(seqtab_nochim)
  )
)

rownames(metadata) <- metadata$SampleID


# ============================================================================ #
# Create raw phyloseq object
# ============================================================================ #

ps_raw <- phyloseq(
  otu_table(seqtab_nochim, taxa_are_rows = FALSE),
  tax_table(taxonomy),
  sample_data(metadata),
  refseq(asv_seq)
)

saveRDS(
  ps_raw,
  file.path(output_dir, "ps_raw.rds")
)


# ============================================================================ #
# Taxonomic filtering
# ============================================================================ #

# Taxa removed in the original analysis.
filterlist <- c(
  "Chloroplast",
  "Mitochondria",
  "Eukaryota",
  "Escherichia-Shigella",
  "Staphylococcus",
  "Lactobacillus",
  "Streptococcus",
  "Corynebacterium"
)

ps_filt <- ps_raw

for (rank in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) {

  if (rank %in% rank_names(ps_filt)) {

    tax_values <- as.character(tax_table(ps_filt)[, rank])

    keep <- is.na(tax_values) | !(tax_values %in% filterlist)

    ps_filt <- prune_taxa(
      taxa_names(ps_filt)[keep],
      ps_filt
    )
  }
}

# Remove ASVs lacking Kingdom- or Phylum-level assignments.
tax_df <- as.data.frame(
  tax_table(ps_filt),
  stringsAsFactors = FALSE
)

keep_assigned <- (
  !is.na(tax_df$Kingdom) &
  trimws(tax_df$Kingdom) != "" &
  !is.na(tax_df$Phylum) &
  trimws(tax_df$Phylum) != ""
)

ps_filt <- prune_taxa(
  rownames(tax_df)[keep_assigned],
  ps_filt
)

# IMPORTANT:
# Singletons were intentionally retained in the original workflow.
# No singleton filter is applied here.

saveRDS(
  ps_filt,
  file.path(output_dir, "ps_filt.rds")
)


# ============================================================================ #
# Export filtered ASV and taxonomy tables for microDecon
# ============================================================================ #

asvtab <- as(
  otu_table(ps_filt),
  "matrix"
)

if (!taxa_are_rows(ps_filt)) {
  asvtab <- t(asvtab)
}

# microDecon input below expects ASVs in rows and samples in columns.
stopifnot(
  identical(
    rownames(asvtab),
    taxa_names(ps_filt)
  )
)

taxonomy_filtered <- as.data.frame(
  tax_table(ps_filt),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

saveRDS(
  asvtab,
  file.path(output_dir, "asvtab.rds")
)

saveRDS(
  taxonomy_filtered,
  file.path(output_dir, "taxonomy_filtered.rds")
)


# ============================================================================ #
# Define the five processing groups
# ============================================================================ #

group_samples <- list(

  G1 = c(
    "YB1_T0",
    "CHF1_T1", "CHF2_T1", "CHF3_T1",
    "HCP1_T1", "HCP2_T1", "HCP3_T1",
    "YSB1_T2", "YSB2_T2", "YSB3_T2",
    "CHF1_T2", "CHF2_T2", "CHF3_T2",
    "HCP1_T2", "HCP2_T2", "HCP3_T2",
    "YSB1_T3", "YSB2_T3", "YSB3_T3",
    "CHF1_T3", "CHF2_T3", "CHF3_T3",
    "YSB1_T4", "YSB2_T4", "YSB3_T4",
    "CHF1_T4", "CHF2_T4", "CHF3_T4",
    "YSB1_T5", "YSB2_T5", "YSB3_T5",
    "CHF1_T5", "CHF2_T5", "CHF3_T5",
    "YSB1_T6", "YSB2_T6", "YSB3_T6",
    "CHF1_T6", "CHF2_T6", "CHF3_T6",
    "YSB1_T7", "YSB2_T7", "YSB3_T7",
    "CHF1_T7", "CHF2_T7", "CHF3_T7",
    "YSB1_T8", "YSB2_T8", "YSB3_T8",
    "CHF1_T8", "CHF2_T8", "CHF3_T8",
    "YSB1_T9", "YSB2_T9", "YSB3_T9",
    "CHF1_T9", "CHF2_T9"
  ),

  G2 = c(
    "YB2_T0",
    "HCP1_T3", "HCP2_T3", "HCP3_T3",
    "HCP1_T4", "HCP2_T4", "HCP3_T4",
    "HCP1_T5", "HCP2_T5", "HCP3_T5",
    "HCP1_T6", "HCP2_T6", "HCP3_T6",
    "HCP1_T7", "HCP2_T7", "HCP3_T7",
    "HCP1_T8", "HCP2_T8", "HCP3_T8"
  ),

  G3 = c(
    "YB3_T0",
    "CHF3_T9",
    "HCP1_T9", "HCP2_T9", "HCP3_T9",
    "YSB1_T10", "YSB2_T10", "YSB3_T10",
    "CHF1_T10", "CHF2_T10", "CHF3_T10",
    "HCP1_T10", "HCP2_T10", "HCP3_T10",
    "YSB1_T11", "YSB2_T11", "YSB3_T11",
    "CHF1_T11", "CHF2_T11", "CHF3_T11",
    "HCP1_T11", "HCP2_T11", "HCP3_T11",
    "YSB1_T12", "YSB2_T12", "YSB3_T12",
    "CHF1_T12", "CHF2_T12", "CHF3_T12",
    "HCP1_T12", "HCP2_T12", "HCP3_T12"
  ),

  G4 = c(
    "YB4_T0",
    "YSB1_T13", "YSB2_T13", "YSB3_T13",
    "CHF1_T13", "CHF2_T13", "CHF3_T13",
    "HCP1_T13", "HCP2_T13", "HCP3_T13",
    "YSB1_T14", "YSB2_T14", "YSB3_T14",
    "CHF1_T14", "CHF2_T14", "CHF3_T14",
    "HCP1_T14", "HCP2_T14", "HCP3_T14",
    "YSB1_T15", "YSB2_T15", "YSB3_T15",
    "CHF1_T15", "CHF2_T15", "CHF3_T15",
    "HCP1_T15", "HCP2_T15", "HCP3_T15",
    "YSB1_T16", "YSB2_T16", "YSB3_T16",
    "CHF1_T16", "CHF2_T16",
    "HCP1_T16", "HCP2_T16", "HCP3_T16"
  ),

  G5 = c(
    "YB5_T0",
    "YSB1_T17", "YSB2_T17", "YSB3_T17",
    "CHF1_T17", "CHF2_T17", "CHF3_T17",
    "HCP1_T17", "HCP2_T17", "HCP3_T17",
    "YSB1_T18", "YSB2_T18", "YSB3_T18",
    "CHF2_T18", "CHF3_T18",
    "HCP1_T18", "HCP2_T18", "HCP3_T18",
    "YSB1_T19", "YSB2_T19", "YSB3_T19",
    "CHF1_T19", "CHF2_T19", "CHF3_T19",
    "HCP1_T19", "HCP2_T19", "HCP3_T19",
    "YSB1_T20", "YSB2_T20", "YSB3_T20",
    "CHF1_T20", "CHF2_T20", "CHF3_T20",
    "HCP1_T20", "HCP2_T20", "HCP3_T20"
  )
)

# Values passed to numb.ind in the original analysis.
group_numb_ind <- list(
  G1 = c(26, 24, 6),
  G2 = c(18),
  G3 = c(10, 9, 12),
  G4 = c(11, 12, 12),
  G5 = c(11, 12, 12)
)


# ============================================================================ #
# Check processing-group sample assignments
# ============================================================================ #

all_group_samples <- unlist(
  group_samples,
  use.names = FALSE
)

duplicated_group_samples <- unique(
  all_group_samples[duplicated(all_group_samples)]
)

if (length(duplicated_group_samples) > 0) {
  stop(
    "Samples occur in more than one microDecon group: ",
    paste(duplicated_group_samples, collapse = ", ")
  )
}

missing_group_samples <- setdiff(
  all_group_samples,
  colnames(asvtab)
)

if (length(missing_group_samples) > 0) {
  stop(
    "Samples required by the original microDecon groups are absent from ",
    "the filtered ASV table: ",
    paste(missing_group_samples, collapse = ", ")
  )
}


# ============================================================================ #
# Prepare taxonomy string for microDecon
# ============================================================================ #

taxonomy_for_decon <- taxonomy_filtered

tax_ranks <- intersect(
  c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  ),
  colnames(taxonomy_for_decon)
)

taxonomy_for_decon$Taxa <- apply(
  taxonomy_for_decon[, tax_ranks, drop = FALSE],
  1,
  function(x) paste(x, collapse = ";")
)


# ============================================================================ #
# Run microDecon separately for the five processing groups
# ============================================================================ #

decon_tables <- vector(
  "list",
  length(group_samples)
)

names(decon_tables) <- names(group_samples)

for (group_name in names(group_samples)) {

  message("Running microDecon for ", group_name, "...")

  samples_this_group <- group_samples[[group_name]]

  group_counts <- as.data.frame(
    asvtab[, samples_this_group, drop = FALSE],
    check.names = FALSE
  )

  group_counts$ASV <- rownames(group_counts)

  final_table <- merge(
    group_counts,
    taxonomy_for_decon["Taxa"],
    by.x = "ASV",
    by.y = "row.names",
    sort = FALSE
  )

  final_table <- final_table[
    ,
    c(
      "ASV",
      samples_this_group,
      "Taxa"
    )
  ]

  decontaminated <- microDecon::decon(
    data = final_table,
    numb.blanks = 1,
    runs = 4,
    thresh = 0.9,
    numb.ind = group_numb_ind[[group_name]],
    taxa = TRUE
  )

  decon_table <- as.data.frame(
    decontaminated[["decon.table"]],
    check.names = FALSE
  )

  saveRDS(
    decontaminated,
    file.path(
      microdecon_dir,
      paste0("decontaminated_", group_name, ".rds")
    )
  )

  saveRDS(
    decon_table,
    file.path(
      microdecon_dir,
      paste0("decont_asvtab_taxa_", group_name, ".rds")
    )
  )

  decon_tables[[group_name]] <- decon_table
}


# ============================================================================ #
# Merge the five microDecon outputs
# ============================================================================ #

decont_asvtab_taxa <- decon_tables |>
  purrr::reduce(
    dplyr::full_join,
    by = c("ASV", "Taxa")
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ tidyr::replace_na(.x, 0)
    )
  ) |>
  dplyr::select(
    ASV,
    Taxa,
    dplyr::everything()
  )

saveRDS(
  decont_asvtab_taxa,
  file.path(
    microdecon_dir,
    "decont_asvtab_taxa.rds"
  )
)

write.table(
  decont_asvtab_taxa,
  file = file.path(
    microdecon_dir,
    "decont_asvtab_taxa.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================================ #
# Unresolved historical step
# ============================================================================ #

# In the original analysis, the merged table above was followed by an
# externally prepared workbook:
#
#   decont_asvtab_taxa_filt.xlsx
#
# with separate sheets named "asvtab" and "taxonomy".
#
# The original master script does not contain the code that generated that
# workbook. Therefore this curated script intentionally stops here rather than
# inventing or altering a scientific filtering step.
#
# Once the historical workbook and/or the exact filtering procedure is
# confirmed, the next script can reconstruct:
#
#   decont_asvtab_taxa
#       -> ps_mdecon_v0
#       -> remove zero-count ASVs
#       -> ps_mdecon_v1
#       -> remove controls
#       -> ps_mdecon_v2
#       -> DECIPHER alignment
#       -> FastTree
#       -> midpoint-rooted final phyloseq object
#
# ============================================================================ #
