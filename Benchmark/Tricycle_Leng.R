# Load required libraries
library(anndata)
library(SingleCellExperiment)
library(scran)
library(tricycle)

# Read the h5ad file containing single-cell data
res <- anndata::read_h5ad('/remote-home/share/data2/wbl/wbl_cellcycle/ground-truth/ref_dataset_Leng.h5ad')

# Extract and transpose expression matrix
expr <- as.matrix(res$X)  # Get expression matrix
expr <- t(expr)  # Transpose to genes x cells format
rownames(expr) <- res$var_names  # Set gene names as row names
colnames(expr) <- rownames(res$obs)  # Set cell names as column names

# Check dimensions for consistency
print(dim(expr))  # Dimensions of expression matrix (genes x cells)
print(length(rownames(res$var)))  # Number of genes in var
print(length(rownames(res$obs)))  # Number of cells in obs

# Update rownames with gene symbols from var metadata
rownames(expr) <- rownames(res$var)

# Create SingleCellExperiment object
# - counts: expression matrix
# - colData: cell metadata (using obs rownames as cell identifiers)
# - rowData: gene metadata (gene symbols)
sce <- SingleCellExperiment(
  assays = list(counts = expr),
  colData = data.frame(cell_id = rownames(res$obs)),  # Store cell IDs
  rowData = data.frame(gene_name = rownames(res$var))  # Store gene symbols
)

# Normalize counts using scran's logNormCounts
sce <- logNormCounts(sce)

# Project cells into cell cycle space using tricycle
# - Uses gene symbols (SYMBOL) and human genome annotation
sce <- project_cycle_space(sce, gname.type = "SYMBOL", species = "human")

# Estimate cell cycle stages using Schwabe model
sce <- estimate_Schwabe_stage(sce, gname.type = "SYMBOL", species = 'human')

# Extract cell cycle stage assignments
ccStage <- as.data.frame(sce$CCStage)

# Get projection coordinates from tricycle embedding
reducedDimNames(sce)  # Check available reduced dimensions
projection_vectors <- reducedDim(sce, "tricycleEmbedding")  # Get coordinates
head(projection_vectors)

# Create plotting dataframe with projection coordinates and original stages
plot_df <- as.data.frame(projection_vectors)
plot_df$stage <- res$obs$stage  # Add original stage annotations

# Preview and save results
print(head(plot_df))
write.csv(plot_df, '/root/Cycle/HarmoCycle_v1/Result/Tricycle/Projection_Leng.csv')
