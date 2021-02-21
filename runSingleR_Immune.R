# Script to run singleR for automated cell type annotation

# use R/4.0
library(SingleR)
library(celldex)
library(Seurat)

# read in seurat object
colon <- readRDS("./clustered_full_colon_proj_seurat.rds")
analysis_parent_folder <- "./immune_normalize_and_scale/"
setwd(analysis_parent_folder)

# set up reference
ref <- HumanPrimaryCellAtlasData()
types_to_use <- c("DC","Epithelial_cells","B_cell","Neutrophils","T_cells","Monocyte",
	"Endothelial_cells","Neurons","Macrophage","NK_cell",
	"BM","Platelets","Fibroblasts","Astrocyte","Myelocyte","Pre-B_cell_CD34-","Pro-B_cell_CD34+","Pro-Myelocyte")
ref <- ref[,(colData(ref)$label.main %in% types_to_use)]

# Run singleR
singler.pred <- SingleR(test = as.SingleCellExperiment(colon), ref = ref, labels = ref$label.fine)

# Add to seurat object the you can plot
colon <- AddMetaData(colon, metadata = singler.pred$labels, col.name = "SingleR.labels")
