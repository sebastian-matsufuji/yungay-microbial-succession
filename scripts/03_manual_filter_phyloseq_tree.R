# ============================================================================ #
# 03_manual_filter_phyloseq_tree.R
# ============================================================================ #
#
# Purpose:
# Rebuild the post-decontamination phyloseq object from the manually curated
# microDecon workbook used in the original analysis, remove negative controls,
# align the retained ASV sequences with DECIPHER, and construct the FastTree
# phylogeny used for downstream analyses.
#
# Historical workflow:
#
#   merged microDecon table
#       -> manual control-informed ASV filtering
#       -> decont_asvtab_taxa_filt.xlsx
#       -> phyloseq
#       -> remove zero-count ASVs
#       -> remove negative controls
#       -> DECIPHER alignment
#       -> FastTree
#       -> midpoint rooting
#
# IMPORTANT:
# - The manual filtering step is NOT re-created here. The curated workbook is
#   treated as an input because the original filtering was performed manually
#   after inspection of the negative control associated with each processing
#   group.
# - Metadata use the final SampleID names and are not renamed in this script.
# - Both the original FastTree phyloseq object and a midpoint-rooted version are
#   saved, matching the historical analysis files.
#
# Inputs:
# - decont_asvtab_taxa_filt.xlsx
#     sheet "asvtab"
#     sheet "taxonomy"
# - final study metadata
# - asv_seq.fasta produced by 02_phyloseq_filter_microdecon.R
#
# Outputs:
# - ps_mdecon_v1.rds
# - ps_mdecon_v2.rds
# - alignment.rds
# - alignment.fasta
# - ps_mdecon_v2_FTree.rds
# - ps_mdecon_v2_FTree_rooted.rds
#
# ============================================================================ #


# ============================================================================ #
# Packages
# ============================================================================ #

library(phyloseq)
library(Biostrings)
library(DECIPHER)
library(openxlsx)
library(janitor)
library(ape)
library(phytools)
library(picante)


# ============================================================================ #
# User paths
# ============================================================================ #

# EDIT THESE PATHS LOCALLY.

manual_filter_path <- "decont_asvtab_taxa_filt.xlsx"
metadata_path <- "metadata.xlsx"
metadata_sheet <- "Sheet1"

sequence_fasta <- file.path(
  "output_phyloseq_microdecon",
  "asv_seq.fasta"
)

output_dir <- "output_phyloseq_tree"
fasttree_dir <- file.path(output_dir, "FastTree")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fasttree_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================ #
# Load manually curated ASV and taxonomy tables
# ============================================================================ #

asvtab <- openxlsx::read.xlsx(
  manual_filter_path,
  sheet = "asvtab",
  rowNames = TRUE
)

taxonomy <- openxlsx::read.xlsx(
  manual_filter_path,
  sheet = "taxonomy",
  rowNames = TRUE
)

asvtab <- as.matrix(asvtab)
taxonomy <- as.matrix(taxonomy)

storage.mode(asvtab) <- "numeric"

if (anyDuplicated(rownames(asvtab))) {
  stop("Duplicated ASV IDs were found in the manually curated ASV table.")
}

if (anyDuplicated(rownames(taxonomy))) {
  stop("Duplicated ASV IDs were found in the manually curated taxonomy table.")
}

if (!setequal(rownames(asvtab), rownames(taxonomy))) {
  stop(
    "ASV IDs in the manually curated count and taxonomy tables do not match."
  )
}

taxonomy <- taxonomy[
  rownames(asvtab),
  ,
  drop = FALSE
]


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

if ("Date" %in% colnames(metadata) && is.numeric(metadata$Date)) {
  metadata$Date <- janitor::excel_numeric_to_date(
    metadata$Date,
    date_system = "modern"
  )
}

missing_metadata <- setdiff(
  colnames(asvtab),
  metadata$SampleID
)

if (length(missing_metadata) > 0) {
  stop(
    "These samples from the curated ASV table are absent from metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}

metadata <- metadata[
  match(colnames(asvtab), metadata$SampleID),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    metadata$SampleID,
    colnames(asvtab)
  )
)

rownames(metadata) <- metadata$SampleID


# ============================================================================ #
# Load representative ASV sequences
# ============================================================================ #

asv_seq <- Biostrings::readDNAStringSet(
  filepath = sequence_fasta,
  format = "fasta"
)

missing_sequences <- setdiff(
  rownames(asvtab),
  names(asv_seq)
)

if (length(missing_sequences) > 0) {
  stop(
    "Representative sequences are missing for these ASVs: ",
    paste(missing_sequences, collapse = ", ")
  )
}

asv_seq <- asv_seq[
  rownames(asvtab)
]

stopifnot(
  identical(
    names(asv_seq),
    rownames(asvtab)
  )
)


# ============================================================================ #
# Create post-microDecon phyloseq object
# ============================================================================ #

ps_mdecon_v0 <- phyloseq(
  otu_table(asvtab, taxa_are_rows = TRUE),
  tax_table(taxonomy),
  sample_data(metadata),
  refseq(asv_seq)
)


# ============================================================================ #
# Remove zero-count ASVs
# ============================================================================ #

ps_mdecon_v1 <- filter_taxa(
  ps_mdecon_v0,
  function(x) sum(x) > 0,
  prune = TRUE
)

saveRDS(
  ps_mdecon_v1,
  file.path(output_dir, "ps_mdecon_v1.rds")
)


# ============================================================================ #
# Remove negative extraction controls
# ============================================================================ #

control_ids <- c(
  "YB1_T0",
  "YB2_T0",
  "YB3_T0",
  "YB4_T0",
  "YB5_T0"
)

controls_present <- intersect(
  control_ids,
  sample_names(ps_mdecon_v1)
)

ps_mdecon_v2 <- prune_samples(
  !sample_names(ps_mdecon_v1) %in% control_ids,
  ps_mdecon_v1
)

saveRDS(
  ps_mdecon_v2,
  file.path(output_dir, "ps_mdecon_v2.rds")
)

cat(
  "Negative controls removed:",
  paste(controls_present, collapse = ", "),
  "\n"
)

cat(
  "Environmental samples retained:",
  nsamples(ps_mdecon_v2),
  "\n"
)

cat(
  "ASVs retained:",
  ntaxa(ps_mdecon_v2),
  "\n"
)


# ============================================================================ #
# DECIPHER multiple-sequence alignment
# ============================================================================ #

seqs_ps <- refseq(ps_mdecon_v2)

Biostrings::writeXStringSet(
  seqs_ps,
  filepath = file.path(
    fasttree_dir,
    "seqs_ps_mdecon_v2.fasta"
  ),
  append = FALSE,
  compress = FALSE,
  format = "fasta"
)

alignment <- DECIPHER::AlignSeqs(
  Biostrings::DNAStringSet(seqs_ps),
  anchor = NA
)

saveRDS(
  alignment,
  file.path(fasttree_dir, "alignment.rds")
)

Biostrings::writeXStringSet(
  alignment,
  filepath = file.path(
    fasttree_dir,
    "alignment.fasta"
  ),
  format = "fasta"
)


# ============================================================================ #
# FastTree phylogeny
# ============================================================================ #

alignment_file <- file.path(
  fasttree_dir,
  "alignment.fasta"
)

tree_file <- file.path(
  fasttree_dir,
  "tree_file"
)

# The original analysis ran FastTree outside R with the following command:
#
#   FastTree -gtr -nt alignment.fasta > tree_file
#
# Run that command from inside the FastTree output directory after this script
# has generated alignment.fasta.
#
# Example:
#
#   cd output_phyloseq_tree/FastTree
#   FastTree -gtr -nt alignment.fasta > tree_file
#
# The script can be run again after tree_file has been created.

if (!file.exists(tree_file)) {

  message(
    "\nDECIPHER alignment completed.\n",
    "FastTree has not yet been run.\n\n",
    "Run the following command inside:\n",
    fasttree_dir,
    "\n\n",
    "FastTree -gtr -nt alignment.fasta > tree_file\n\n",
    "Then run this script again to import, root, and save the tree.\n"
  )

} else {

  # -------------------------------------------------------------------------- #
  # Import unrooted FastTree tree
  # -------------------------------------------------------------------------- #

  tree <- ape::read.tree(
    tree_file
  )

  if (!setequal(tree$tip.label, taxa_names(ps_mdecon_v2))) {
    stop(
      "FastTree tip labels do not match the ASVs in ps_mdecon_v2."
    )
  }

  ps_mdecon_v2_FTree <- merge_phyloseq(
    ps_mdecon_v2,
    tree
  )

  saveRDS(
    ps_mdecon_v2_FTree,
    file.path(
      output_dir,
      "ps_mdecon_v2_FTree.rds"
    )
  )


  # ========================================================================== #
  # Midpoint root and reconcile tree with phyloseq taxa
  # ========================================================================== #

  tree_rooted <- phytools::midpoint.root(
    tree
  )

  taxa_reference <- data.frame(
    present = rep(1, ntaxa(ps_mdecon_v2)),
    row.names = taxa_names(ps_mdecon_v2)
  )

  matched <- picante::match.phylo.data(
    tree_rooted,
    taxa_reference
  )

  tree_rooted <- matched$phy

  stopifnot(
    setequal(
      tree_rooted$tip.label,
      taxa_names(ps_mdecon_v2)
    )
  )

  ps_mdecon_v2_FTree_rooted <- merge_phyloseq(
    ps_mdecon_v2,
    tree_rooted
  )

  saveRDS(
    ps_mdecon_v2_FTree_rooted,
    file.path(
      output_dir,
      "ps_mdecon_v2_FTree_rooted.rds"
    )
  )


  # ========================================================================== #
  # Final summary
  # ========================================================================== #

  cat("\nPhylogenetic-tree processing completed.\n")
  cat(
    "Samples:",
    nsamples(ps_mdecon_v2_FTree_rooted),
    "\n"
  )
  cat(
    "ASVs:",
    ntaxa(ps_mdecon_v2_FTree_rooted),
    "\n"
  )
  cat(
    "Tree tips:",
    length(phy_tree(ps_mdecon_v2_FTree_rooted)$tip.label),
    "\n"
  )
  cat(
    "Tree rooted:",
    ape::is.rooted(phy_tree(ps_mdecon_v2_FTree_rooted)),
    "\n"
  )
}
