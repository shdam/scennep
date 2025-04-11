## code to prepare `example_counts` dataset
## This script simulates a small scRNA-seq count matrix
## for use in package examples and tests.

# --- Configuration ---
set.seed(123) 

n_genes <- 200    # Number of genes (rows)
n_cells <- 90    # Number of cells (columns)
n_groups <- 3    # Number of distinct cell groups/types

# --- Simulate Cell Groups ---
group_sizes <- diff(c(0, sort(sample(1:(n_cells-1), n_groups - 1)), n_cells))
cell_group_ids <- rep(paste0("Group", 1:n_groups), times = group_sizes)
cell_order <- sample(1:n_cells)
cell_group_ids <- cell_group_ids[order(cell_order)]
cell_names <- paste0("Cell_", 1:n_cells)
names(cell_group_ids) <- cell_names

# --- Simulate Gene Properties ---
gene_names <- paste0("Gene", 1:n_genes)

# Baseline expression
gene_baseline_meanlog <- rnorm(n_genes, mean = -2, sd = 1)
gene_baseline_means <- exp(gene_baseline_meanlog)

# Dispersion parameter
gene_dispersions_log <- rnorm(n_genes, mean = log(10), sd = 0.5)
gene_dispersions <- exp(gene_dispersions_log)

# Marker Genes
n_marker_genes_per_group <- 10
marker_gene_indices <- sample(1:n_genes, n_groups * n_marker_genes_per_group, replace = FALSE)
marker_list <- split(marker_gene_indices, rep(1:n_groups, each = n_marker_genes_per_group))
names(marker_list) <- paste0("Group", 1:n_groups)
marker_effect_multiplier <- 20

# --- Simulate Counts ---
example_counts <- matrix(0, nrow = n_genes, ncol = n_cells)
rownames(example_counts) <- gene_names
colnames(example_counts) <- cell_names

# Cell size factors
cell_size_factors_log <- rnorm(n_cells, mean = 0, sd = 0.5)
cell_size_factors <- exp(cell_size_factors_log)
names(cell_size_factors) <- cell_names

# Generate expected means (lambda) and sample counts
for (j in 1:n_cells) {
    cell_name <- cell_names[j]
    group_id <- cell_group_ids[cell_name]
    
    expected_means <- gene_baseline_means * cell_size_factors[cell_name]
    
    if (group_id %in% names(marker_list)) {
        markers_for_this_group <- marker_list[[group_id]]
        expected_means[markers_for_this_group] <- expected_means[markers_for_this_group] * marker_effect_multiplier
    }
    
    expected_means <- expected_means + 1e-6 
    
    example_counts[, j] <- rnbinom(n_genes, size = gene_dispersions, mu = expected_means)
}

# --- Create Metadata DataFrames (not used) ---
cell_metadata <- data.frame(
    cell_id = cell_names,
    group = factor(cell_group_ids), # Store groups as factor
    simulated_size_factor = cell_size_factors,
    total_counts = colSums(example_counts),
    row.names = cell_names
)

gene_metadata <- data.frame(
    gene_id = gene_names,
    simulated_baseline_mean = gene_baseline_means,
    simulated_dispersion = gene_dispersions,
    is_marker = gene_names %in% gene_names[unlist(marker_list)],
    marker_group = NA_character_, # Add column for which group it marks (if any)
    total_counts = rowSums(example_counts),
    row.names = gene_names
)

# Populate the marker_group column
for (grp_name in names(marker_list)) {
    gene_indices <- marker_list[[grp_name]]
    gene_metadata$marker_group[gene_indices] <- grp_name
}
gene_metadata$marker_group <- factor(gene_metadata$marker_group) # Store as factor

# --- Store objects ---


# Save counts and metadata in a list
# small_scrna_data <- list(
#     counts = counts_matrix,
#     cell_meta = cell_metadata,
#     gene_meta = gene_metadata
# )

usethis::use_data(example_counts, overwrite = TRUE)
