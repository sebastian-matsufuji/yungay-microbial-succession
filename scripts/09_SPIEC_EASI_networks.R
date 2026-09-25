# ============================================================================ #
# 09_SPIEC_EASI_networks.R
# ============================================================================ #
#
# Purpose:
# Reproduce the SPIEC-EASI association-network analyses supporting Figure 7
# and Supplementary Table 7.
#
# Historical detail preserved from the original analysis:
# - HCP Wet = T1-T9 and HCP Dry = T10-T20.
# - CHF and YSB use the metadata Phase variable.
#
# The exact final code used to select the five hub nodes reported in
# Supplementary Table 8 is not present in the master script. Therefore this
# curated script exports all node-centrality metrics but does not invent a new
# retrospective hub-ranking rule.
#
# Input:
# - output_phyloseq_tree/ps_mdecon_v2_FTree.rds
#
# Main outputs:
# - six fitted SPIEC-EASI models
# - network/node/edge summaries
# - Figure 7A-D
# - Table_Supplementary7_network_topology.xlsx
# - Supplementary_Table8_all_node_centralities.xlsx
#
# ============================================================================ #


# ============================================================================ #
# Packages
# ============================================================================ #

library(phyloseq)
library(microbiomeutilities)
library(SpiecEasi)
library(pulsar)
library(genefilter)
library(igraph)
library(tidygraph)
library(ggraph)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(openxlsx)


# ============================================================================ #
# Reproducibility and paths
# ============================================================================ #

set.seed(2026)

phyloseq_path <- file.path(
  "output_phyloseq_tree",
  "ps_mdecon_v2_FTree.rds"
)

output_dir <- "output_SPIEC_EASI"
model_dir <- file.path(output_dir, "models")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")

dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================ #
# Shared settings
# ============================================================================ #

basin_levels <- c("HCP", "CHF", "YSB")
phase_levels <- c("Wet", "Dry")

basin_labels <- c(
  "HCP" = "Herradura Clay Pan",
  "CHF" = "Ckoirama Halite Field",
  "YSB" = "Yungay Station Basin"
)

phase_colors <- c(
  "Wet" = "#4C78A8",
  "Dry" = "#D59A3A"
)

edge_colors <- c(
  "Positive" = "#007A5E",
  "Negative" = "#B83A00"
)

order_levels_rank <- c(
  "Other", "Unclassified", "Nanosalinales", "Longimicrobiales",
  "Bradymonadales", "Oscillospirales", "Lachnospirales", "Deinococcales",
  "Paenibacillales", "Bacillales", "Mycobacteriales", "Propionibacteriales",
  "Micrococcales", "Bacteroidales", "Cytophagales", "Sphingobacteriales",
  "Flavobacteriales", "Azospirillales", "Rhodobacterales", "Hyphomicrobiales",
  "Caulobacterales", "Sphingomonadales", "Cardiobacteriales", "Enterobacterales",
  "Lysobacterales", "Burkholderiales", "Pseudomonadales", "Halobacterales"
)

selected_orders <- setdiff(
  order_levels_rank,
  c("Other", "Unclassified")
)

order_colors <- c(
  "Other" = "#9E9E9E",
  "Unclassified" = "#343A40",
  "Nanosalinales" = "gold",
  "Longimicrobiales" = "chartreuse4",
  "Bradymonadales" = "lightpink",
  "Oscillospirales" = "peachpuff4",
  "Lachnospirales" = "lightgoldenrod",
  "Deinococcales" = "#ae017e",
  "Paenibacillales" = "#F4A261",
  "Bacillales" = "#e6550d",
  "Mycobacteriales" = "darkred",
  "Propionibacteriales" = "#fc9272",
  "Micrococcales" = "firebrick2",
  "Bacteroidales" = "olivedrab1",
  "Cytophagales" = "limegreen",
  "Sphingobacteriales" = "#66C2A5",
  "Flavobacteriales" = "#d9f0a3",
  "Azospirillales" = "#c7e9b4",
  "Rhodobacterales" = "#2c7fb8",
  "Hyphomicrobiales" = "#41b6c4",
  "Caulobacterales" = "#a1dab4",
  "Sphingomonadales" = "#253494",
  "Cardiobacteriales" = "#AAB0E6",
  "Enterobacterales" = "#D6BCEB",
  "Lysobacterales" = "#810f7c",
  "Burkholderiales" = "#dadaeb",
  "Pseudomonadales" = "#7A1FA2",
  "Halobacterales" = "#48D1CC"
)


# ============================================================================ #
# Load and preprocess phyloseq object
# ============================================================================ #

ps <- readRDS(phyloseq_path)

if (!taxa_are_rows(ps)) {
  otu_table(ps) <- t(otu_table(ps))
}

ps_besthit <- microbiomeutilities::format_to_besthit(ps)

if (!all(taxa_names(ps) == taxa_names(ps_besthit))) {
  taxa_names(ps_besthit) <- taxa_names(ps)
}

tax_table(ps_besthit) <- tax_table(ps_besthit)[, 1:7]

colnames(tax_table(ps_besthit)) <- c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

# Species-rank agglomeration with unresolved assignments retained.
ps_species <- tax_glom(
  ps_besthit,
  taxrank = "Species",
  NArm = FALSE
)

saveRDS(
  ps_species,
  file.path(model_dir, "ps_species_agglomerated.rds")
)


# ============================================================================ #
# Historical basin-phase subsets
# ============================================================================ #

make_network_subset <- function(ps_object, basin_code, phase_name) {

  basin_name <- basin_labels[[basin_code]]

  ps_basin <- subset_samples(
    ps_object,
    Ephemeral_Basin == basin_name
  )

  if (basin_code == "HCP") {

    if (phase_name == "Wet") {
      ps_sub <- subset_samples(
        ps_basin,
        Time %in% paste0("T", 1:9)
      )
    } else {
      ps_sub <- subset_samples(
        ps_basin,
        Time %in% paste0("T", 10:20)
      )
    }

  } else {

    ps_sub <- subset_samples(
      ps_basin,
      Phase == phase_name
    )
  }

  ps_sub <- prune_taxa(
    taxa_sums(ps_sub) > 0,
    ps_sub
  )

  ps_sub
}


# ============================================================================ #
# Node filtering
# ============================================================================ #

filter_network_nodes <- function(ps_sub) {

  # Relative abundance is used only for node filtering.
  ps_ra <- transform_sample_counts(
    ps_sub,
    function(x) 100 * x / sum(x)
  )

  otu_table(ps_ra)[is.na(otu_table(ps_ra))] <- 0

  # Historical filter: >0.01% relative abundance in at least one sample.
  keep_taxa <- filter_taxa(
    ps_ra,
    genefilter::kOverA(
      k = 1,
      A = 0.01
    ),
    prune = FALSE
  )

  list(
    ps_ra = prune_taxa(keep_taxa, ps_ra),
    ps_counts = prune_taxa(keep_taxa, ps_sub)
  )
}


# ============================================================================ #
# Infer one SPIEC-EASI network
# ============================================================================ #

infer_network <- function(ps_sub, basin_code, phase_name) {

  phase_code <- ifelse(
    phase_name == "Wet",
    "w",
    "d"
  )

  network_id <- paste0(
    basin_code,
    "_",
    phase_code
  )

  message("Running SPIEC-EASI: ", network_id)

  filtered <- filter_network_nodes(ps_sub)
  ps_ra_filt <- filtered$ps_ra
  ps_counts_filt <- filtered$ps_counts

  set.seed(2026)

  se_mb <- SpiecEasi::spiec.easi(
    ps_counts_filt,
    method = "mb",
    sel.criterion = "stars",
    verbose = TRUE,
    pulsar.select = TRUE,
    lambda.min.ratio = 1e-3,
    nlambda = 50,
    pulsar.params = list(
      rep.num = 100,
      thresh = 0.05,
      subsample.ratio = 0.8,
      seed = 2026
    )
  )

  saveRDS(
    se_mb,
    file.path(
      model_dir,
      paste0("se.mb_", network_id, ".rds")
    )
  )

  graph <- SpiecEasi::adj2igraph(
    SpiecEasi::getRefit(se_mb),
    vertex.attr = list(
      name = taxa_names(ps_ra_filt)
    )
  )

  beta_mat <- as.matrix(
    SpiecEasi::symBeta(
      SpiecEasi::getOptBeta(se_mb)
    )
  )

  edges <- as.data.frame(
    igraph::as_edgelist(graph),
    stringsAsFactors = FALSE
  )

  colnames(edges) <- c("from", "to")

  edges$weight <- apply(
    edges,
    1,
    function(z) {
      i <- match(z["from"], taxa_names(ps_ra_filt))
      j <- match(z["to"], taxa_names(ps_ra_filt))
      beta_mat[i, j]
    }
  )

  edges <- edges |>
    mutate(
      Basin = basin_code,
      Phase = phase_name,
      Network_ID = network_id,
      Sign = case_when(
        weight > 0 ~ "Positive",
        weight < 0 ~ "Negative",
        TRUE ~ NA_character_
      ),
      abs_weight = abs(weight)
    )

  tax_df <- as.data.frame(
    tax_table(ps_ra_filt)
  ) |>
    rownames_to_column("ASV")

  tax_df$Order <- as.character(tax_df$Order)
  tax_df$Order <- sub("^o__", "", tax_df$Order)
  tax_df$Order[is.na(tax_df$Order) | tax_df$Order == ""] <- "Unclassified"

  tax_df$Order <- ifelse(
    tax_df$Order %in% selected_orders |
      tax_df$Order == "Unclassified",
    tax_df$Order,
    "Other"
  )

  otu_ra <- as(otu_table(ps_ra_filt), "matrix")

  if (taxa_are_rows(ps_ra_filt)) {
    otu_ra <- t(otu_ra)
  }

  mean_ra <- apply(
    otu_ra,
    2,
    mean
  )

  node_abundance <- log2(mean_ra)

  graph_tbl <- graph |>
    as_tbl_graph() |>
    activate(edges) |>
    mutate(
      Sign = edges$Sign,
      weight = edges$weight
    ) |>
    activate(nodes) |>
    left_join(
      tax_df,
      by = c("name" = "ASV")
    ) |>
    mutate(
      abundance = node_abundance[
        match(name, names(node_abundance))
      ],
      degree = centrality_degree(),
      betweenness = centrality_betweenness(),
      eigenvector = centrality_eigen(),
      closeness = centrality_closeness()
    )

  nodes <- graph_tbl |>
    activate(nodes) |>
    as_tibble() |>
    transmute(
      Basin = basin_code,
      Phase = phase_name,
      Network_ID = network_id,
      ASV = name,
      Kingdom,
      Phylum,
      Class,
      Order,
      Family,
      Genus,
      Species,
      degree,
      betweenness,
      eigenvector,
      closeness,
      abundance
    )

  comp <- igraph::components(graph)
  clusters <- igraph::cluster_louvain(graph)

  n_nodes <- igraph::vcount(graph)
  n_edges <- igraph::ecount(graph)
  n_positive <- sum(edges$Sign == "Positive", na.rm = TRUE)
  n_negative <- sum(edges$Sign == "Negative", na.rm = TRUE)

  network_summary <- tibble(
    Basin = basin_code,
    Phase = phase_name,
    Network_ID = network_id,
    Samples = nsamples(ps_counts_filt),
    Nodes = n_nodes,
    Edges = n_edges,
    Positive_Interactions = n_positive,
    Negative_Interactions = n_negative,
    Positive_pct = 100 * n_positive / n_edges,
    Negative_pct = 100 * n_negative / n_edges,
    Density = igraph::edge_density(graph),
    Global_transitivity = igraph::transitivity(
      graph,
      type = "global"
    ),
    Diameter = igraph::diameter(
      graph,
      directed = FALSE,
      weights = NA
    ),
    Mean_path_length = igraph::mean_distance(
      graph,
      directed = FALSE,
      weights = NA
    ),
    Components = comp$no,
    Largest_component = max(comp$csize),
    Largest_component_pct = 100 * max(comp$csize) / n_nodes,
    Assortativity_degree = igraph::assortativity_degree(
      graph,
      directed = FALSE
    ),
    Mean_degree = mean(nodes$degree),
    Mean_betweenness = mean(nodes$betweenness),
    Mean_eigenvector = mean(nodes$eigenvector),
    Mean_closeness = mean(nodes$closeness, na.rm = TRUE),
    Modules = length(unique(igraph::membership(clusters))),
    Modularity = igraph::modularity(clusters),
    Mean_abs_weight = mean(edges$abs_weight, na.rm = TRUE),
    Min_abs_weight = min(edges$abs_weight, na.rm = TRUE),
    Max_abs_weight = max(edges$abs_weight, na.rm = TRUE)
  )

  opt_idx <- se_mb$select$stars$opt.index

  stars_summary <- tibble(
    Basin = basin_code,
    Phase = phase_name,
    Network_ID = network_id,
    Method = "mb",
    Optimal_lambda_position = opt_idx,
    Stability = se_mb$select$stars$summary[opt_idx]
  )

  openxlsx::write.xlsx(
    list(
      Network_summary = network_summary,
      Node_metrics = nodes,
      Edge_table = edges,
      STARS = stars_summary
    ),
    file.path(
      table_dir,
      paste0("SPIEC_EASI_", network_id, ".xlsx")
    ),
    overwrite = TRUE
  )

  list(
    id = network_id,
    basin = basin_code,
    phase = phase_name,
    model = se_mb,
    graph = graph,
    graph_tbl = graph_tbl,
    nodes = nodes,
    edges = edges,
    summary = network_summary,
    stars = stars_summary
  )
}


# ============================================================================ #
# Run the six networks
# ============================================================================ #

network_design <- tidyr::crossing(
  Basin = basin_levels,
  Phase = phase_levels
)

networks <- list()

for (i in seq_len(nrow(network_design))) {

  basin_i <- network_design$Basin[i]
  phase_i <- network_design$Phase[i]

  ps_sub <- make_network_subset(
    ps_species,
    basin_i,
    phase_i
  )

  result_i <- infer_network(
    ps_sub,
    basin_i,
    phase_i
  )

  networks[[result_i$id]] <- result_i
}


# ============================================================================ #
# Compile network outputs
# ============================================================================ #

network_summary_all <- map_dfr(
  networks,
  "summary"
) |>
  mutate(
    Basin = factor(Basin, levels = basin_levels),
    Phase = factor(Phase, levels = phase_levels)
  ) |>
  arrange(Basin, Phase)

node_metrics_all <- map_dfr(
  networks,
  "nodes"
) |>
  mutate(
    Basin = factor(Basin, levels = basin_levels),
    Phase = factor(Phase, levels = phase_levels)
  )

edge_table_all <- map_dfr(
  networks,
  "edges"
)

stars_all <- map_dfr(
  networks,
  "stars"
)

openxlsx::write.xlsx(
  list(
    Network_summary = network_summary_all,
    Node_metrics = node_metrics_all,
    Edge_table = edge_table_all,
    STARS = stars_all
  ),
  file.path(
    table_dir,
    "SPIEC_EASI_network_data_compiled.xlsx"
  ),
  overwrite = TRUE
)


# ============================================================================ #
# Figure 7A | Six association networks
# ============================================================================ #

plot_network <- function(net, show_legend = FALSE) {

  set.seed(2026)

  layout_fr <- create_layout(
    net$graph_tbl,
    layout = "fr",
    weights = abs(
      net$graph_tbl |>
        activate(edges) |>
        pull(weight)
    )
  )

  ggraph(layout_fr) +
    geom_edge_link(
      aes(
        edge_colour = Sign,
        edge_width = abs(weight)
      ),
      alpha = 0.8
    ) +
    scale_edge_colour_manual(
      values = edge_colors,
      name = "Association"
    ) +
    scale_edge_width(
      range = c(0.3, 1.5),
      limits = c(0, 0.55),
      breaks = seq(0.1, 0.5, 0.1),
      name = "Weight"
    ) +
    geom_node_point(
      aes(
        fill = Order,
        size = abundance,
        shape = Kingdom
      ),
      color = "black",
      stroke = 0.7
    ) +
    scale_size_continuous(
      range = c(1.5, 7),
      guide = "none"
    ) +
    scale_fill_manual(
      values = order_colors,
      breaks = order_levels_rank,
      name = "Order"
    ) +
    scale_shape_manual(
      values = c(
        "Bacteria" = 21,
        "Archaea" = 24
      ),
      name = "Domain"
    ) +
    labs(
      title = paste0(
        basin_labels[[net$basin]],
        "\n",
        net$phase
      )
    ) +
    theme_void() +
    theme(
      plot.title = element_text(
        size = 11,
        face = "bold",
        hjust = 0.5
      ),
      legend.position = if (show_legend) "right" else "none"
    )
}

p_HCP_w <- plot_network(networks$HCP_w)
p_CHF_w <- plot_network(networks$CHF_w)
p_YSB_w <- plot_network(networks$YSB_w)
p_HCP_d <- plot_network(networks$HCP_d)
p_CHF_d <- plot_network(networks$CHF_d, show_legend = TRUE)
p_YSB_d <- plot_network(networks$YSB_d)

figure7a <- (
  p_HCP_w + p_CHF_w + p_YSB_w
) /
  (
    p_HCP_d + p_CHF_d + p_YSB_d
  ) +
  patchwork::plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(figure_dir, "Figure7A_SPIEC_EASI_networks.tiff"),
  plot = figure7a,
  width = 18,
  height = 10,
  units = "in",
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  file.path(figure_dir, "Figure7A_SPIEC_EASI_networks.svg"),
  plot = figure7a,
  width = 18,
  height = 10,
  units = "in",
  device = "svg"
)


# ============================================================================ #
# Figure 7B | Degree of nodes shared between Wet and Dry networks
# ============================================================================ #

shared_degree <- inner_join(
  node_metrics_all |>
    filter(Phase == "Wet") |>
    select(
      Basin,
      ASV,
      degree_wet = degree
    ),
  node_metrics_all |>
    filter(Phase == "Dry") |>
    select(
      Basin,
      ASV,
      degree_dry = degree
    ),
  by = c("Basin", "ASV")
)

degree_tests <- shared_degree |>
  group_by(Basin) |>
  group_modify(
    ~ {
      test <- wilcox.test(
        .x$degree_dry,
        .x$degree_wet,
        paired = TRUE,
        exact = FALSE
      )

      tibble(
        Shared_ASVs = nrow(.x),
        Median_degree_wet = median(.x$degree_wet),
        Median_degree_dry = median(.x$degree_dry),
        P_value = test$p.value
      )
    }
  ) |>
  ungroup() |>
  mutate(
    P_adjust_BH = p.adjust(
      P_value,
      method = "BH"
    ),
    Significance = case_when(
      P_adjust_BH < 0.001 ~ "***",
      P_adjust_BH < 0.01 ~ "**",
      P_adjust_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

shared_degree_long <- shared_degree |>
  pivot_longer(
    cols = c(degree_wet, degree_dry),
    names_to = "Phase",
    values_to = "Degree"
  ) |>
  mutate(
    Phase = recode(
      Phase,
      degree_wet = "Wet",
      degree_dry = "Dry"
    ),
    Phase = factor(
      Phase,
      levels = phase_levels
    ),
    Basin = factor(
      Basin,
      levels = basin_levels
    )
  )

figure7b <- ggplot(
  shared_degree_long,
  aes(
    x = Phase,
    y = Degree,
    fill = Phase
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_point(
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    size = 1.2,
    alpha = 0.4
  ) +
  facet_wrap(
    ~ Basin,
    nrow = 1,
    scales = "free_y",
    labeller = as_labeller(basin_labels)
  ) +
  scale_fill_manual(
    values = phase_colors
  ) +
  labs(
    x = NULL,
    y = "Degree"
  ) +
  theme_cowplot() +
  theme(
    legend.position = "none",
    strip.background = element_blank()
  )

ggsave(
  file.path(figure_dir, "Figure7B_shared_ASV_degree.tiff"),
  plot = figure7b,
  width = 7.5,
  height = 3.8,
  units = "in",
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  file.path(figure_dir, "Figure7B_shared_ASV_degree.svg"),
  plot = figure7b,
  width = 7.5,
  height = 3.8,
  units = "in",
  device = "svg"
)


# ============================================================================ #
# Figure 7C | Edge rewiring among shared nodes
# ============================================================================ #

add_edge_key <- function(df) {

  df |>
    mutate(
      node1 = pmin(from, to),
      node2 = pmax(from, to),
      Edge_key = paste(
        node1,
        node2,
        sep = "||"
      )
    )
}

calculate_rewiring <- function(basin_code) {

  wet <- networks[[paste0(basin_code, "_w")]]
  dry <- networks[[paste0(basin_code, "_d")]]

  shared_nodes <- intersect(
    wet$nodes$ASV,
    dry$nodes$ASV
  )

  wet_edges <- wet$edges |>
    filter(
      from %in% shared_nodes,
      to %in% shared_nodes
    ) |>
    add_edge_key()

  dry_edges <- dry$edges |>
    filter(
      from %in% shared_nodes,
      to %in% shared_nodes
    ) |>
    add_edge_key()

  wet_keys <- unique(wet_edges$Edge_key)
  dry_keys <- unique(dry_edges$Edge_key)

  conserved <- intersect(wet_keys, dry_keys)
  lost <- setdiff(wet_keys, dry_keys)
  gained <- setdiff(dry_keys, wet_keys)
  edge_union <- union(wet_keys, dry_keys)

  n_union <- length(edge_union)

  tibble(
    Basin = basin_code,
    Shared_ASVs = length(shared_nodes),
    Wet_edges_shared_ASVs = length(wet_keys),
    Dry_edges_shared_ASVs = length(dry_keys),
    Conserved = length(conserved),
    Lost_in_Dry = length(lost),
    Gained_in_Dry = length(gained),
    Edge_union = n_union,
    Jaccard = length(conserved) / n_union,
    Rewiring = 1 - Jaccard,
    Conserved_pct = 100 * length(conserved) / n_union,
    Lost_in_Dry_pct = 100 * length(lost) / n_union,
    Gained_in_Dry_pct = 100 * length(gained) / n_union
  )
}

rewiring_summary <- map_dfr(
  basin_levels,
  calculate_rewiring
)

rewiring_plot_data <- rewiring_summary |>
  select(
    Basin,
    Conserved_pct,
    Lost_in_Dry_pct,
    Gained_in_Dry_pct
  ) |>
  pivot_longer(
    cols = -Basin,
    names_to = "Edge_status",
    values_to = "Percentage"
  ) |>
  mutate(
    Edge_status = recode(
      Edge_status,
      Conserved_pct = "Conserved",
      Lost_in_Dry_pct = "Lost in Dry",
      Gained_in_Dry_pct = "Gained in Dry"
    ),
    Edge_status = factor(
      Edge_status,
      levels = c(
        "Lost in Dry",
        "Conserved",
        "Gained in Dry"
      )
    ),
    Basin = factor(
      Basin,
      levels = basin_levels
    )
  )

rewiring_colors <- c(
  "Lost in Dry" = unname(phase_colors["Wet"]),
  "Conserved" = "#BDBDBD",
  "Gained in Dry" = unname(phase_colors["Dry"])
)

figure7c <- ggplot(
  rewiring_plot_data,
  aes(
    x = Basin,
    y = Percentage,
    fill = Edge_status
  )
) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f%%",
        Percentage
      )
    ),
    position = position_stack(
      vjust = 0.5
    ),
    size = 3.4
  ) +
  scale_fill_manual(
    values = rewiring_colors
  ) +
  scale_x_discrete(
    labels = basin_labels
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    x = NULL,
    y = "Edges in wet-dry union (%)",
    fill = NULL
  ) +
  theme_cowplot()

ggsave(
  file.path(figure_dir, "Figure7C_edge_rewiring.tiff"),
  plot = figure7c,
  width = 7.2,
  height = 4.4,
  units = "in",
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  file.path(figure_dir, "Figure7C_edge_rewiring.svg"),
  plot = figure7c,
  width = 7.2,
  height = 4.4,
  units = "in",
  device = "svg"
)


# ============================================================================ #
# Figure 7D | Positive inferred associations
# ============================================================================ #

positive_plot_data <- network_summary_all |>
  mutate(
    Basin = factor(
      Basin,
      levels = basin_levels
    ),
    Phase = factor(
      Phase,
      levels = phase_levels
    )
  )

figure7d <- ggplot(
  positive_plot_data,
  aes(
    x = Positive_pct,
    y = Basin,
    color = Phase
  )
) +
  geom_line(
    aes(group = Basin),
    color = "grey65",
    linewidth = 0.8
  ) +
  geom_point(
    size = 3.5
  ) +
  scale_color_manual(
    values = phase_colors
  ) +
  scale_y_discrete(
    labels = basin_labels
  ) +
  labs(
    x = "Positive inferred associations (%)",
    y = NULL,
    color = NULL
  ) +
  theme_cowplot()

ggsave(
  file.path(figure_dir, "Figure7D_positive_associations.tiff"),
  plot = figure7d,
  width = 7.2,
  height = 4.0,
  units = "in",
  device = "tiff",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  file.path(figure_dir, "Figure7D_positive_associations.svg"),
  plot = figure7d,
  width = 7.2,
  height = 4.0,
  units = "in",
  device = "svg"
)


# ============================================================================ #
# Supplementary Table 7
# ============================================================================ #

network_wide <- network_summary_all |>
  select(
    Basin,
    Phase,
    Nodes,
    Edges,
    Mean_degree,
    Global_transitivity,
    Mean_path_length,
    Largest_component_pct,
    Positive_pct
  ) |>
  pivot_wider(
    names_from = Phase,
    values_from = c(
      Nodes,
      Edges,
      Mean_degree,
      Global_transitivity,
      Mean_path_length,
      Largest_component_pct,
      Positive_pct
    )
  )

supp_table7 <- network_wide |>
  left_join(
    degree_tests,
    by = "Basin"
  ) |>
  left_join(
    rewiring_summary,
    by = "Basin"
  )

openxlsx::write.xlsx(
  supp_table7,
  file.path(
    table_dir,
    "Table_Supplementary7_network_topology.xlsx"
  ),
  overwrite = TRUE
)


# ============================================================================ #
# Supplementary Table 8 support
# ============================================================================ #

# The historical master contains an older exploratory hub-selection rule:
# top 20 nodes by degree, followed by top 10 by betweenness.
#
# The final manuscript instead reports five hubs per network selected using
# degree, betweenness, and eigenvector centrality. The exact final ranking code
# is not present in the master script. For transparency, export every node and
# all four centrality metrics without recreating an undocumented ranking rule.

openxlsx::write.xlsx(
  node_metrics_all,
  file.path(
    table_dir,
    "Supplementary_Table8_all_node_centralities.xlsx"
  ),
  overwrite = TRUE
)


# ============================================================================ #
# Supporting exports
# ============================================================================ #

openxlsx::write.xlsx(
  list(
    Shared_degree = shared_degree,
    Degree_tests = degree_tests
  ),
  file.path(
    table_dir,
    "Figure7B_shared_ASV_degree_data.xlsx"
  ),
  overwrite = TRUE
)

openxlsx::write.xlsx(
  list(
    Rewiring_summary = rewiring_summary,
    Plot_data = rewiring_plot_data
  ),
  file.path(
    table_dir,
    "Figure7C_edge_rewiring_data.xlsx"
  ),
  overwrite = TRUE
)

openxlsx::write.xlsx(
  positive_plot_data,
  file.path(
    table_dir,
    "Figure7D_positive_associations_data.xlsx"
  ),
  overwrite = TRUE
)


# ============================================================================ #
# Final summary
# ============================================================================ #

cat("\nSPIEC-EASI analyses completed.\n\n")

print(
  network_summary_all |>
    select(
      Basin,
      Phase,
      Nodes,
      Edges,
      Mean_degree,
      Global_transitivity,
      Mean_path_length,
      Largest_component_pct,
      Positive_pct
    )
)

cat("\nShared-node paired degree tests:\n")
print(degree_tests)

cat("\nEdge rewiring among shared nodes:\n")
print(rewiring_summary)
