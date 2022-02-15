# Script to analyze stromal scRNA data from scHTAN and scHuBMAP Projects
# WRB 2020

##############################################################################################################################
# Import libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ArchR)
library(viridis)
library(DoubletFinder)
library(Rcpp)
library(harmony)
library(future)
library(Matrix)
set.seed(1)

##############################################################################################################################
# Analysis steps
# 1) QC
# 2) Normalize and scale data, plot UMAP,find markers
# 3) Remove bad clusters and redo analysis
# 4) ID Cluster Markers
# 5) Plot Markers
# 6) Initial cluster identification
# 7) Find markers

execute_steps <- c(1,2,3,4,5,6,7)

# Define variables
sample_name <- "all_samples"
analysis_parent_folder <- "./stromal_analysis/"

# Create directory
if (!dir.exists(paste0(analysis_parent_folder))){
	dir.create(paste0(analysis_parent_folder))
}
setwd(analysis_parent_folder)

colon <- readRDS("./initial_clustering/initial_clustering_stromal.rds")

###############################################################################################################################
# Define Functions
vln_plot <- function(features, save_name){
	pdf(save_name, width = 20, onefile=F)
	print(VlnPlot(colon, features = features, group.by = "orig.ident", pt.size = 0, cols = c(rep("#D51F26",100)))+
	geom_boxplot(outlier.shape = NA, alpha = 0.6)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))+
	scale_x_discrete(labels=paste0(data.frame(table(colon@meta.data$orig.ident))$Var1, "\n n = ",  data.frame(table(colon@meta.data$orig.ident))$Freq)))
	dev.off()
}

normalize_and_dim_reduce <- function (colon, sample_name){
	colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
	colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(colon)
	colon <- ScaleData(colon, features = all.genes)
	colon <- RunPCA(colon, features = VariableFeatures(object = colon))
	print(pdf(paste0("./pca", sample_name, ".pdf")))
	DimPlot(colon, reduction = "pca")
	dev.off()
	colon <- FindNeighbors(colon, dims = 1:20)
	colon <- FindClusters(colon, resolution = 0.5)
	colon <- RunUMAP(colon, dims = 1:20)
	return(colon)
}

plotUMAP <- function(colon){
	pdf(paste0("./UMAPclustering" , ".pdf"))
	print(DimPlot(colon, reduction = "umap", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAP_samples.pdf"), width = 12)
	print(DimPlot(colon, reduction = "umap", group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))#, cols = (ArchRPalettes$stallion))
	dev.off()

	pdf(paste0("./UMAP_disease_state.pdf"), width = 6.5, onefile=F)
	print(DimPlot(colon, reduction = "umap", group.by = "DiseaseState",
		cols = c("#D51F26", "#272E6A", "#208A42", "#89288F" ,"#F47D2B" ,"#FEE500" ,"#8A9FD1", "#C06CAB","#90D5E4", "#89C75F")) + theme_ArchR())
	dev.off()
}

seurat_feature_plot <- function(sample_name, reduction, cell_type, markers){
	p1 <- FeaturePlot(colon, features = markers, reduction = reduction, sort.cell = TRUE, combine = FALSE, pt.size = 5)
	fix.sc <- scale_colour_gradientn(colours = ArchRPalettes$blueYellow)
	if (length(p1)==1){
		width <- 4
		height <- 4
	} else if (length(p1)==2){
		width <- 8
		height <- 4
	} else if (length(p1)<5){
		width <- 8
		height <- 8
	} else if (length(p1)<7){
		width <- 12
		height <- 8
	} else if (length(p1)<10){
		width <- 12
		height <- 12
	} else if (length(p1)<13){
		width <- 16
		height <- 12
	} else if (length(p1)<17){
		width <- 16
		height <- 16
	}
	pdf(paste0("./", reduction, "_feature_plot_", sample_name, "_", cell_type ,".pdf"), width = width, height = height)
	print(CombinePlots(lapply(p1, function (x) AugmentPlot(x + fix.sc))))
	dev.off()
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 1) QC
if (1 %in% execute_steps){
	if (!dir.exists(paste0(analysis_parent_folder, "all_samples_qc_plots"))){
		dir.create(paste0(analysis_parent_folder, "all_samples_qc_plots"))
	}
	setwd(paste0(analysis_parent_folder, "all_samples_qc_plots"))

	colon[["percent.mt"]] <- PercentageFeatureSet(colon, pattern = "^MT-")

	# Now subset the project and make some nice qc plots
	colon <- subset(colon, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 10000)
	vln_plot(c("nFeature_RNA"), paste0("./n_genes_violin_", sample_name, ".pdf"))
	vln_plot(c("nCount_RNA"), paste0("./n_counts_violin_", sample_name, ".pdf"))
	vln_plot(c("percent.mt"), paste0("./pMT_violin_", sample_name, ".pdf"))

	setwd(analysis_parent_folder)
}


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 2) Normalize and scale data, plot UMAP, find markers
if (2 %in% execute_steps){
	colon <- normalize_and_dim_reduce(colon, sample_name)
	plotUMAP(colon)
	colon.markers <- FindAllMarkers(colon, only.pos = TRUE, min.pct = 0.1, max.cells.per.ident = 250)
}


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 3) Remove bad clusters and redo analysis
if (3 %in% execute_steps){
	# Save ids of cells to remove
	bad_clusters <- bad_clusters <- c(3,12,16)

	# Subset
	colon_new <- DietSeurat(subset(colon, subset = seurat_clusters %ni% bad_clusters))

	# Redo the initial analysis steps
	makeQCplots(colon_new, analysis_parent_folder, sample_name)
	colon_new <- normalize_and_dim_reduce(colon_new, sample_name)
	colon_new <- FindClusters(colon_new, resolution = 1.0)
	plotUMAP(colon_new)
	colon <- colon_new

	# Save ids of cells to remove
	bad_clusters <- c(5)

	# Subset
	colon_new <- DietSeurat(subset(colon, subset = seurat_clusters %ni% bad_clusters))

	# Redo the initial analysis steps
	makeQCplots(colon_new, analysis_parent_folder, sample_name)
	colon_new <- normalize_and_dim_reduce(colon_new, sample_name)
	colon_new <- FindClusters(colon_new, resolution = 1.0)
	plotUMAP(colon_new)
	colon <- colon_new
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 4) ID Cluster Markers
if (4 %in% execute_steps){
	colon.markers <- FindAllMarkers(colon, only.pos = TRUE, min.pct = 0.05, max.cells.per.ident = 250)

	pdf(paste0("./Cell_marker_heatmap_", sample_name ,".pdf"), width = 12, height = 20)
	DoHeatmap(colon, features = top10$gene, label=FALSE, group.colors = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)) + scale_fill_viridis()
	dev.off()
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 5) Plot Markers
if (5 %in% execute_steps){
	reductions_to_plot <- c("umap")
	for (reduction in reductions_to_plot){
		seurat_feature_plot(sample_name, reduction, "Fibroblasts", c("COL1A1", "COL1A2", "COL6A1", "COL6A2", "FAP", "CBLN2", "SPOCK1", "ACSS3"))
		seurat_feature_plot(sample_name, reduction, "Fibroblasts_cell_subtypes", c("RSPO3", "CCL11", "WNT5B", "BMP4", "CHI3L1", "ACTA2", "WNT2B"))
		seurat_feature_plot(sample_name, reduction, "Myofibroblasts", c("SYT10", "SOSTDC1", "DES", "MYH11", "TAGLN", "ACTA2", "TPM4"))
		seurat_feature_plot(sample_name, reduction, "stromal_other", c("FAM110D", "INHBB", "NPR1", "NOVA2", "GPIHBP1", "SOX17", "VWF", "PLVAP", "CDH5", "S100B"))
		seurat_feature_plot(sample_name, reduction, "Pericytes", c("MCAM", "COX4I2", "KCNJ8", "HIGD1B", "RGS5", "NOTCH3", "HEYL", "FAM162B"))
		seurat_feature_plot(sample_name, reduction, "Microvascular", c("PLVAP","CD36","DYSF","NRP1","SH3BP5","EXOC3L2","FABP5","VWA1","BAALC","PRSS23","RAPGEF4","APLN","HTRA1"))
		seurat_feature_plot(sample_name, reduction, "SchwanCell", c("S100A1", "SOX10", "EGR2", "MBP", "MPZ", "GAP43", "NCAM", "P75NTR"))
		seurat_feature_plot(sample_name, reduction, "Nerve", c("MAP2", "RBFOX3", "DLG4", "SYP"))
	}
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 6) Initial cluster identification
if (6 %in% execute_steps){
	colon <- FindClusters(colon, resolution = 1.0)
	new.cluster.ids <- c(
		"Myofibroblasts/Smooth Muscle 1",#0
		"Myofibroblasts/Smooth Muscle 1", #1
		"Crypt Fibroblasts 1", #2
		"Myofibroblasts/Smooth Muscle 1", #3
		"Villus Fibroblasts WNT5B+", #4
		"Endothelial", #5
		"Crypt Fibroblasts 3", #6
		"Myofibroblasts/Smooth Muscle 1", #7
		"Myofibroblasts/Smooth Muscle 1", #8
		"Myofibroblasts/Smooth Muscle 3", #9 GREM1+
		"Crypt Fibroblasts 2", #10
		"Crypt Fibroblasts 4", #11 RSPO3+
		"Myofibroblasts/Smooth Muscle 2", #12
		"Cancer Associated Fibroblasts", #13
		"Pericytes", #14
		"Glia", #15
		"Lymphatic Endothelial Cells", #16
		"Endothelial",#17 Post-capillary venules
		"Endothelial", #18
		"Crypt Fibroblasts 4",#19 RSPO3+
		"Unknown", #20 #neurons
		"Adipocytes", #21
		"Myofibroblasts/Smooth Muscle 1", #22
		"Neurons")#23

	identities <- as.character(colon@meta.data$seurat_clusters)
	for (i in 0:length(new.cluster.ids)){
		identities[identities==as.character(i)] <- new.cluster.ids[i+1]
	}
	colon <- AddMetaData(colon, identities, col.name = "CellType")
	pdf(paste0("./UMAP_cell_type.pdf"), width = 6.2)
	AugmentPlot(DimPlot(colon, reduction = "umap", group.by = "CellType", cols = paletteDiscrete(values = unique(colon@meta.data$CellType), set = "stallion", reverse = FALSE))+theme_ArchR(),dpi = 500)#, cols = (ArchRPalettes$stallion))
	dev.off()

	pdf(paste0("./UMAP_cell_type_vector.pdf"), width = 6.2, onefile = F)
	DimPlot(colon, reduction = "umap", group.by = "CellType", cols = paletteDiscrete(values = unique(colon@meta.data$CellType), set = "stallion", reverse = FALSE))+theme_ArchR()
	dev.off()

	saveRDS(colon, file = "colon_stromal_all_samples_filtered.rds")
	saveRDS(DietSeurat(colon), "diet_clustered_full_colon_proj_seurat.rds")
}


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 7) Find markers
if (7 %in% execute_steps){
	new.cluster.ids <- c(
		"Myofibroblasts/Smooth Muscle 1", #0
		"Myofibroblasts/Smooth Muscle 1", #1
		"Crypt Fibroblasts 1", #2
		"Myofibroblasts/Smooth Muscle 1", #3
		"Villus Fibroblasts WNT5B+", #4
		"Endothelial", #5
		"Crypt Fibroblasts 3", #6
		"Myofibroblasts/Smooth Muscle 1", #7
		"Myofibroblasts/Smooth Muscle 1", #8
		"Myofibroblasts/Smooth Muscle 3", #9 GREM1+
		"Crypt Fibroblasts 2", #10
		"Crypt Fibroblasts 4", #11 RSPO3+
		"Myofibroblasts/Smooth Muscle 2", #12
		"Cancer Associated Fibroblasts", #13
		"Pericytes", #14
		"Glia", #15
		"Lymphatic Endothelial Cells", #16
		"Endothelial",#17 Post-capillary venules
		"Endothelial", #18
		"Crypt Fibroblasts 4",#19 RSPO3+
		"Unknown", #20 #neurons
		"Adipocytes", #21
		"Myofibroblasts/Smooth Muscle 1", #22
		"Neurons")#23
	names(new.cluster.ids) <- levels(colon)
	colon <- RenameIdents(colon, new.cluster.ids)

	# ID markers
	colon.markers <- FindAllMarkers(colon_named, only.pos = TRUE, min.pct = 0.05, max.cells.per.ident = 250, test.use = "MAST")
	colon.markers.CAF <- colon.markers[colon.markers$cluster == "Cancer Associated Fibroblasts",]
	colon.up <- colon.markers.CAF[colon.markers.CAF$avg_logFC>0,]
}


