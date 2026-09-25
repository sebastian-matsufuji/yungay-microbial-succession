Yungay microbial succession after rare rainfall in the hyperarid Atacama Desert
This repository contains the analysis code associated with the manuscript:
From Rain to Brine to Crust: Divergent Prokaryotic Community Assembly Trajectories during Desiccation of Ephemeral (Hyper)saline Lakes in the Hyperarid Atacama Desert
Raw 16S rRNA gene amplicon reads are deposited in the NCBI Sequence Read Archive under BioProject PRJNA1527592.
Repository structure
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
Analysis workflow
The scripts cover sequence processing, alpha and beta diversity, taxonomic composition, CCA and eLSA community–environment analyses, SPIEC-EASI association networks, iCAMP and betaNTI/RCBray community assembly analyses, FAPROTAX and PICRUSt2 functional inference, and physicochemical/Mantel analyses.
External software and manual steps
Some analyses include steps outside R, including FastTree v2.2.0, eLSA, FAPROTAX, and PICRUSt2 v2.6.3. A manually reviewed control-informed filtering step is also part of the decontamination workflow. The corresponding scripts document the inputs and outputs used for these steps.
Reproducibility notes
The scripts were curated from the analysis workflow used to generate the manuscript figures, tables, and statistical results. They preserve the analytical parameters and calculations used in the study while removing obsolete exploratory blocks, duplicated code, and machine-specific absolute paths.
The repository should be interpreted as the analysis code associated with the manuscript rather than as a fully containerized one-command workflow.
Software environment
R 4.5.3
Bioconductor 3.22
Windows 10 x64
Package requirements are listed through library() calls in the individual scripts.
Citation
A permanent software citation and DOI will be added after creation of the archived repository release.
License
MIT License.
