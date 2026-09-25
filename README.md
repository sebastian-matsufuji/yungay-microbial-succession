# Yungay microbial succession after rare rainfall in the hyperarid Atacama Desert

This repository contains the analysis code associated with the manuscript:

**From Rain to Brine to Crust: Divergent Prokaryotic Community Assembly Trajectories during Desiccation of Ephemeral (Hyper)saline Lakes in the Hyperarid Atacama Desert**

The study examines prokaryotic community succession and assembly across three ephemeral saline basins in the Yungay region of the hyperarid Atacama Desert following an exceptional rainfall event in June 2017.

## Data availability

Raw 16S rRNA gene amplicon reads are deposited in the **NCBI Sequence Read Archive (SRA)** under BioProject:

**PRJNA1527592**

Raw sequencing files are therefore not duplicated in this repository.

## Repository structure

```text
scripts/
├── 01_DADA2_taxonomy.R
├── 02_phyloseq_filter_microdecon.R
├── 03_manual_filter_phyloseq_tree.R
├── 04_alpha_diversity.R
├── 05_taxonomic_composition_LEfSe.R
├── 06_temporal_beta_diversity.R
├── 07_CCA_community_environment.R
├── 08_eLSA_family_environment.R
├── 09_SPIEC_EASI_networks.R
├── 10_iCAMP_community_assembly.R
├── 11_betaNTI_RCbray_community_assembly.R
├── 12_FAPROTAX_functional_inference.R
├── 13_PICRUSt2_MetaCyc_LinDA.R
├── 14_physicochemical_environment.R
├── 15_environmental_beta_diversity.R
└── 16_Mantel_community_environment.R
```

## Analysis workflow

### Sequence processing

**01_DADA2_taxonomy.R**  
DADA2 processing of 16S rRNA gene amplicon sequences, including filtering, dereplication, ASV inference, paired-end merging, chimera removal, and taxonomic assignment using SILVA.

**02_phyloseq_filter_microdecon.R**  
Construction and initial filtering of the phyloseq dataset, including removal of non-target taxa and control-informed decontamination with microDecon.

**03_manual_filter_phyloseq_tree.R**  
Integration of the manually reviewed ASV table, removal of controls and zero-abundance taxa, sequence alignment, and incorporation of the FastTree phylogeny.

### Alpha diversity

**04_alpha_diversity.R**  
Repeated-rarefaction estimates of observed ASV richness, Shannon diversity, and Faith's phylogenetic diversity, including temporal and Wet–Dry comparisons.

### Taxonomic composition

**05_taxonomic_composition_LEfSe.R**  
Order-level taxonomic composition and family-level LEfSe analyses.

### Temporal beta diversity

**06_temporal_beta_diversity.R**  
Temporal changes in Bray–Curtis and weighted UniFrac community dissimilarity, including similarity to the initial community and temporal distance-decay analyses.

### Community–environment relationships

**07_CCA_community_environment.R**  
Canonical correspondence analyses relating community composition to physicochemical and geochemical gradients.

**08_eLSA_family_environment.R**  
Preparation and analysis of family-level temporal associations with physicochemical variables using Extended Local Similarity Analysis (eLSA).

### Association networks

**09_SPIEC_EASI_networks.R**  
SPIEC-EASI association networks, network topology, shared-node changes, and Wet–Dry edge rewiring.

### Community assembly

**10_iCAMP_community_assembly.R**  
Basin-specific iCAMP analysis of ecological assembly processes and their temporal Wet–Dry variation.

**11_betaNTI_RCbray_community_assembly.R**  
Abundance-weighted betaMNTD, betaNTI, RCBray, and ecological-process classification of temporal community turnover.

### Functional inference

**12_FAPROTAX_functional_inference.R**  
FAPROTAX functional annotation workflow and statistical analysis of Wet–Dry differences in predicted functional groups.

**13_PICRUSt2_MetaCyc_LinDA.R**  
PICRUSt2 MetaCyc pathway prediction workflow and LinDA differential-abundance analyses.

### Environmental analyses

**14_physicochemical_environment.R**  
Temporal physicochemical variation, Wet–Dry basin comparisons, total organic carbon summaries, and environmental correlation matrices.

**15_environmental_beta_diversity.R**  
Relationships between pairwise differences in water content, salinity, pH, and ORP and Bray–Curtis or weighted UniFrac community dissimilarity.

**16_Mantel_community_environment.R**  
Basin-specific Mantel analyses linking community dissimilarity with physicochemical and geochemical variation, including BH-FDR correction and temporal sensitivity analyses.

## External software and manual steps

Some parts of the workflow require software or procedures outside R:

- **FastTree v2.2.0** for phylogenetic tree inference.
- **microDecon** for control-informed decontamination.
- A manually reviewed filtering step following decontamination.
- **eLSA** for temporal family–environment association analysis.
- **FAPROTAX** for functional annotation.
- **PICRUSt2 v2.6.3** for MetaCyc pathway prediction.

The corresponding R scripts document the relevant inputs and outputs for these steps.

## Reproducibility notes

The scripts were curated from the analysis workflow used to generate the manuscript figures, tables, and statistical results.

They preserve the principal analytical parameters and calculations used in the study while removing obsolete exploratory code, duplicated blocks, and machine-specific absolute paths.

This repository should be interpreted as the analysis code associated with the manuscript rather than as a fully containerized one-command workflow.

## Software environment

The final analyses were based primarily on:

- R 4.5.3
- Bioconductor 3.22
- Windows 10 x64

Individual R package requirements are specified through `library()` calls within the scripts.

## Citation

A permanent software citation and DOI will be added after creation of the archived repository release.

## License

This repository is distributed under the **MIT License**.
