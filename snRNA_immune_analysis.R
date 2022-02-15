# Script to analyze immune subset of scRNA data from scHTAN and scHuBMAP Projects
# WRB 2020-2021

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
library(Rmagic)
library(ggpubr)

# Analysis steps
# 1) QC and filtering (filtering should be done previously)
# 2) Normalize and scale data
# 3) Plot UMAPs
# 5) ID Cluster Markers and compare to published data
# 6) Remove bad clusters and redo analysis
# 7) Plot Markers
# 8) Cluster identification 
# 9) Cell type fraction plots
# 10) Module scores for IFN-gamma signaling

# Set steps to run
execute_steps <- c(1,2,3,5,6,7,9,10)

# Load previously defined seurat object that contains only the cells classified as immune cells
colon <- readRDS("./colon_immune_all_samples_initial.rds") # this is the one with all the cells from the initial clustering post doublet finder doublet removal

# Define variables
sample_name <- "all_samples" # used as label for saving plots
analysis_parent_folder <- "./immune_results/"
setwd(analysis_parent_folder)

###############################################################################################################################
# Define Functions #
###############################################################################################################################
vln_plot <- function(features, save_name){
	pdf(save_name, width = 20, onefile=F)
	print(VlnPlot(colon, features = features, group.by = "orig.ident", pt.size = 0, cols = c(rep("#D51F26",100)))+
	geom_boxplot(outlier.shape = NA, alpha = 0.6)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))+
	scale_x_discrete(labels=paste0(data.frame(table(colon@meta.data$orig.ident))$Var1, "\n n = ",  data.frame(table(colon@meta.data$orig.ident))$Freq)))
	dev.off()
}

normalize_and_dim_reduce <- function (colon, sample_name){
	# stanfard seurat pipeline
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

plotUMAPandRunHarmony <- function(colon, run_harmony, version){
	# plot umaps
	pdf(paste0("./UMAPclustering", version, ".pdf"))
	print(DimPlot(colon, reduction = "umap", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAP_samples", version, ".pdf"), width = 12)
	print(DimPlot(colon, reduction = "umap", group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAP_disease_state", version, ".pdf"), width = 6.5, onefile=F)
	print(DimPlot(colon, reduction = "umap", group.by = "DiseaseState",
		cols = c("#D51F26", "#272E6A", "#208A42", "#89288F" ,"#F47D2B" ,"#FEE500" ,"#8A9FD1", "#C06CAB","#90D5E4", "#89C75F")) + theme_ArchR())
	dev.off()

	# run harmony
	if (run_harmony){
		library(harmony)
		colon <- RunHarmony(colon, "orig.ident") # will use PCA
		colon <- RunUMAP(colon, dims = 1:20, reduction = "harmony", reduction.name = "umapharmony")
		pdf(paste0("./UMAPharmony_samples.pdf"), width = 12)
		print(DimPlot(colon, reduction = "umapharmony", group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))#, cols = (ArchRPalettes$stallion))
		dev.off()

		pdf(paste0("./UMAPharmony_diseasestate.pdf"), width = 12)
		print(DimPlot(colon, reduction = "umapharmony", group.by = "DiseaseState", cols = paletteDiscrete(values = unique(colon@meta.data$DiseaseState), set = "stallion", reverse = FALSE)))#, cols = (ArchRPalettes$stallion))
		dev.off()

		pdf(paste0("./UMAPharmony_clustering.pdf"), width = 12)
		print(DimPlot(colon, reduction = "umapharmony", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
		dev.off()

		colon <- FindNeighbors(colon, reduction = "harmony", dims = 1:20)
		colon <- FindClusters(colon, resolution = 2.0)

		pdf(paste0("./UMAPharmony_clustering_new.pdf"), width = 12)
		print(DimPlot(colon, reduction = "umapharmony", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
		dev.off()
	}
	return(colon)
}

seurat_feature_plot <- function(colon, sample_name, reduction, cell_type, markers){
	# make seurat feature plots for multiple markers and arrange into a grid
	# colon: seurat_object
	# sample_name: included in the plot save name
	# reduction: reduction to plot e.g. UMAP
	# cell_type: the name of the cell type the markers correspond to, will be added to plot name
	# markers: list of arkers to plot
	p1 <- FeaturePlot(colon, features = markers, reduction = reduction, sort.cell = TRUE, combine = FALSE, pt.size = 4)
	fix.sc <- scale_colour_gradientn(colours = ArchRPalettes$greyMagma)
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
	# colon <- subset(colon, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 10000) # Already done
	vln_plot(c("nFeature_RNA"), paste0("./n_genes_violin_", sample_name, ".pdf"))
	vln_plot(c("nCount_RNA"), paste0("./n_counts_violin_", sample_name, ".pdf"))
	vln_plot(c("percent.mt"), paste0("./pMT_violin_", sample_name, ".pdf"))

	setwd(analysis_parent_folder)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 2) Normalize and scale data
if (2 %in% execute_steps){
	# This is here to reproduce the analysis exactly becuase the name for this changed, but the name for this has to be changed 
	# back at the end to match the metadata, you can remove these lines but the result will be slightly different
	colon@meta.data$orig.ident[colon@meta.data$orig.ident == "A002-C-010-R0"]<- "A002-C-025-R0"
	colon <- RenameCells(colon, new.names = gsub("A002-C-010-R0", "A002-C-025-R0", colnames(colon)))

	colon <- normalize_and_dim_reduce(colon, sample_name)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 3) Plot UMAPs and run harmony
if (3 %in% execute_steps){
	colon <- plotUMAPandRunHarmony(colon, TRUE, version = "initial")
	saveRDS(colon, "clustered_full_colon_immune_proj_seurat_pre_doublet_filter.rds")
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 5) ID Cluster Markers
if (5 %in% execute_steps){
	colon.markers <- FindAllMarkers(colon, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, max.cells.per.ident = 250)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 6) Remove bad clusters and redo analysis
if (6 %in% execute_steps){
	# Save ids of cells to remove
	bad_clusters <- c(7,21,23,24)

	# Subset seurat object
	colon_new <- DietSeurat(subset(colon, subset = seurat_clusters %ni% bad_clusters))

	# Rerun find clusters w/ desired resolution
	colon_new <- normalize_and_dim_reduce(colon_new, sample_name)
	colon_new <- plotUMAPandRunHarmony(colon_new, FALSE, version = "bad_clusters_removed")
	colon <- colon_new
	colon <- FindClusters(colon, resolution = 4.0)

	# Plot the clustering of the filtered project
	pdf(paste0("./UMAP_clustering_bad_clusters_removed.pdf"), width = 12)
	plot = (DimPlot(colon, reduction = "umap", group.by = "RNA_snn_res.4", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
	plot = LabelClusters(plot = plot, id = "RNA_snn_res.4")
	print(plot)
	dev.off()

	saveRDS(colon, "clustered_full_colon_immune_proj_seurat_post_doublet_filter.rds")
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 7) Plot Markers
if (7 %in% execute_steps){
	reductions_to_plot <- c("umap")
	for (reduction in reductions_to_plot){
		seurat_feature_plot(colon, sample_name, reduction, "B_cells", c("PAX5", "MS4A1", "CD19", "IGLL5", "VPREB3"))
		seurat_feature_plot(colon, sample_name, reduction, "GC_B_cells", c("SERPINA9", "HRK", "HTR3A", "TCL6", "CD180", "FCRLA"))
		seurat_feature_plot(colon, sample_name, reduction, "Plasma_B_cells", c("SSR4", "IGLL5", "IGLL1", "AMPD1"))
		seurat_feature_plot(colon, sample_name, reduction, "Mast_cells", c("TPSAB1", "HDC", "CTSG", "CMA1", "KRT1", "IL1RAPL1", "GATA2"))
		seurat_feature_plot(colon, sample_name, reduction, "CD69pos_Mast", c("CMA1", "IL1RAPL1", "CD69"))
		seurat_feature_plot(colon, sample_name, reduction, "NK", c("KLRF1", "SH2D1B", "SH2D1B", "NCAM1", "FCGR3A")) # NCAM1 and FCGR3A from Blish lab paper
		seurat_feature_plot(colon, sample_name, reduction, "Monocytes_macrophages", c("CD14", "CLEC9A", "FCGR1A", "LILRB2", "CD209", "CD1E", "FOLR2","FABP3","PLA2G2D"))
		seurat_feature_plot(colon, sample_name, reduction, "T_cells", c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "TBX21", "IL7R", "CD4", "CD2"))
		seurat_feature_plot(colon, sample_name, reduction, "Tregs", c("BATF","TNFRSF4", "FOXP3","CTLA4","LAIR2"))
		seurat_feature_plot(colon, sample_name, reduction, "T_memory", c("BACH2", "IFNG", "STIM2", "ID2", "IFNAR1", "IL12RB2", "PTPRC"))
		seurat_feature_plot(colon, sample_name, reduction, "T_naive", c("CCR7", "CD28", "ETS1"))
		seurat_feature_plot(colon, sample_name, reduction, "activated_CD4", c("IL4R", "STAT1", "MAL", "SOCS1", "IL2", "ODC1", "WARS"))
		seurat_feature_plot(colon, sample_name, reduction, "T_activated", c("TNF", "IFNG", "JUN", "FOS", "CD69", "REL"))
		seurat_feature_plot(colon, sample_name, reduction, "Th17_CD4", c("IL17A", "CTSH", "KLRB1", "IL26"))
		seurat_feature_plot(colon, sample_name, reduction, "T_exhauseted", c("PDCD1", "HAVCR2", "LAG3","CD101", "CD38", "CXCR6", "TIGIT"))
		seurat_feature_plot(colon, sample_name, reduction, "T_term_exhausted", c("TOX", "GZMB", "ENTPD1", "ITGAE"))
	}
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 8) Cluster identification - overclustered and combined clusters with similar markers
if (8 %in% execute_steps){
	# change back sample names if changed before
	colon@meta.data$orig.ident[colon@meta.data$orig.ident == "A002-C-025-R0"]<- "A002-C-010-R0"
	colon <- RenameCells(colon, new.names = gsub("A002-C-025-R0", "A002-C-010-R0", colnames(colon)))

	new.cluster.ids <- c(
		"Memory B", #0-
		"Naive B", #1-
		"Naive B", #2-
		"Plasma", #3-
		"Memory B", #4-
		"Plasma", #5-
		"Naive B", #6-
		"Macrophages", #7-
		"Memory B", #8-
		"Memory B", #9-
		"Naive T", #10-
		"CD4+", #11-
		"Naive B", #12-
		"CD4+", #13-
		"CD8+", #14-
		"Naive T", #15-
		"CD8+", #16-
		"Tregs", #17-
		"Macrophages", #18-
		"CD4+", #19-
		"CD4+", #20-
		"CD8+", #21-
		"GC", #22-
		"CD4+", #23-
		"Macrophages", #24-
		"CD4+", #25-
		"Macrophages", #26-
		"CD4+", #27-
		"GC", #28-
		"Memory B", #29-
		"Macrophages", #30-
		"GC", #31-
		"NK", #32-
		"GC", #33-
		"Macrophages", #34-
		"Macrophages", #35-
		"GC", #36-
		"GC", #37-
		"Plasma", #38-
		"CD4+", #39-
		"Macrophages", #40-
		"Macrophages", #41-
		"Mast", #42-
		"ILCs", #43-
		"DC", #44-
		"CD8+", #45-
		"CD4+", #46-
		"GC", #47-
		"Macrophages" #48-
		)

	pal <- c("#D51F26")
	names(pal) <- "CD4+"
	pal["CD8+"] <- "#89288F"
	pal["DC"] <- "#F47D2B"
	pal["GC"] <- "#8A9FD1"
	pal["ILCs"] <- "#C06CAB"
	pal["Macrophages"] <- "#D8A767"
	pal["Mast"] <- "#89C75F"
	pal["Memory B"] <- "#F37B7D"
	pal["Naive B"] <- "#9983BD"
	pal["NK"] <- "#D24B27"
	pal["Plasma"] <- "#3BBCA8"
	pal["Tregs"] <- "#6E4B9E"
	pal["Naive T"] <- "#0C727C" 

	identities <- as.character(colon@meta.data$seurat_clusters)
	for (i in 0:length(new.cluster.ids)){
		identities[identities==as.character(i)] <- new.cluster.ids[i+1]
	}
	colon <- AddMetaData(colon, identities, col.name = "CellType")
	pdf(paste0("./UMAP_cell_type_initial.pdf"), width = 7, onefile = F)
	DimPlot(colon, reduction = "umap", pt.size = 0.4, group.by = "CellType", cols = pal) + theme_ArchR()#, cols = (ArchRPalettes$stallion))
	dev.off()

	##############################################################################################################################
	##############################################################################################################################
	# save
	saveRDS(colon, "clustered_full_colon_immune_proj_seurat_post_doublet_filter.rds")
	saveRDS(DietSeurat(colon), "diet_clustered_full_colon_immune_proj_seurat_post_doublet_filter.rds") 
}
				  
##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 9) Cell type fraction plots
if (9 %in% execute_steps){
	# cell type colored by disease state fraction plot
	cM <- confusionMatrix(paste0(colon@meta.data$CellType), paste0(colon@meta.data$orig.ident))
	cM <- data.frame(cM / Matrix::rowSums(cM))
	samplenames <- c()
	cellnames <- c()
	values <- c()
	diseaseStates <- c()
	for (sample in colnames(cM)){
	    diseaseState <- unique(colon@meta.data[colon@meta.data$orig.ident == gsub("\\.", "-", sample),]$DiseaseState)
	for (celltype in rownames(cM)){
	    samplenames <- c(samplenames, sample)
	    cellnames <- c(cellnames, celltype)
	    diseaseStates <- c(diseaseStates, diseaseState)
	    values <- c(values, cM[celltype,sample])
	}
	}
	data <- data.frame(samplenames,cellnames,values, diseaseStates)
	data <- data[order(sapply(data$diseaseStates, function(x) which(x == c("Normal", "Unaffected", "Polyp", "Adenocarcinoma")))), ]

	colOrder <-  c("CD4+","Tregs","CD8+","Naive T","NK","ILCs","DC","Macrophages","Mast","Memory B","Naive B","GC","Plasma")

	data <- data[order(sapply(data$cellnames, function(x) which(x == colOrder))), ]
	data$cellnames <- factor(data$cellnames, levels = colOrder)

	counts <- table(colon@meta.data$CellType)[colOrder]
	axis_labels <- paste0(levels(data$cellnames), " (N = ", counts, ")")

	#data <- data[order(data$diseaseStates),]
	data$samplenames <- factor(data$samplenames, levels = unique(data$samplenames))
	pdf(paste0(paste("Sample-CellType-Bar-Fraction_Samplev2", sample_name, sep = "-"), ".pdf"), width = 8, height = 12, onefile=F)
	p <- ggplot(data, aes(fill=samplenames, y=values, x=cellnames)) + 
	geom_bar(position="stack", stat="identity") + theme_ArchR() + 
	ylab("Fraction Cells") + scale_x_discrete(labels= axis_labels)+
	xlab("Sample") + theme(axis.text.x = element_text(angle = 90)) + theme(axis.text=element_text(size=20,hjust=0.95,vjust=0.2),
	    axis.title=element_text(size=20)) +
	scale_fill_manual("legend", values = c(
			    colorRampPalette(c("#007849","#1E392A"))(8), 
			    colorRampPalette(c("#008F95", "#0375B4","#062F4F"))(17),
			    colorRampPalette(c("#94618E", "#6E3667","#813772"))(36),
			    colorRampPalette(c("#E24E42", "#fe3401","#B82601"))(5)))
	print(p)
	dev.off()

	##############################################################################################################################
	##############################################################################################################################
	# Make celltype fraction plot with samples on x axis and stacked bars colored by cell type
	cM <- confusionMatrix(paste0(colon@meta.data$CellType), paste0(colon@meta.data$orig.ident))
	cM <- data.frame(t(cM) / Matrix::rowSums(t(cM)))
	samplenames <- c()
	cellnames <- c()
	values <- c()
	diseaseStates <- c()
	for (sample in rownames(cM)){
	    diseaseState <- unique(colon@meta.data[colon@meta.data$orig.ident == sample,]$DiseaseState)
	    for (celltype in colnames(cM)){
	        samplenames <- c(samplenames, sample)
	        cellnames <- c(cellnames, celltype)
	        diseaseStates <- c(diseaseStates, diseaseState)
	        values <- c(values, cM[sample, celltype])
	    }
	}
	data <- data.frame(samplenames,cellnames,values, diseaseStates)
	data <- data[order(samplenames), ]
	data <- data[order(sapply(data$diseaseStates, function(x) which(x == c("Normal", "Unaffected", "Polyp", "Adenocarcinoma")))), ]
	data$samplenames <- factor(data$samplenames, levels = unique(data$samplenames))

	test <- data.frame(table(colon@meta.data$orig.ident))
	rownames(test) <- test$Var1

	pdf(paste0(paste("Sample-CellType-Bar-Fraction_Sample", sample_name, sep = "-"), ".pdf"), width = 20, height = 7, onefile=F)
	ggplot(data, aes(fill=cellnames, y=values, x=samplenames)) + 
	    geom_bar(position="stack", stat="identity") + theme_ArchR() + 
	  ylab("Fraction Cells") +
	  xlab("Sample") + theme(axis.text.x = element_text(angle = 90)) +
	  scale_x_discrete(labels=paste0(levels(data$samplenames), " (n = ",  test[levels(data$samplenames),]$Freq, ")")) +
	  scale_fill_manual("legend", values = c("#D51F26", "#89288F", "#F47D2B", "#8A9FD1", "#C06CAB", "#D8A767", "#89C75F", "#F37B7D", "#9983BD", "#0C727C", "#D24B27", "#3BBCA8", "#6E4B9E"))
	dev.off()

	##############################################################################################################################
	##############################################################################################################################
	# Make celltype fraction plot with celltypes on x axis and stacked bars colored by sample/donor
	cM <- confusionMatrix(paste0(colon@meta.data$CellType), paste0(colon@meta.data$orig.ident))
	cM <- data.frame(cM / Matrix::rowSums(cM))
	samplenames <- c()
	cellnames <- c()
	values <- c()
	donors <- c()
	for (sample in colnames(cM)){
		Donor <- unique(colon@meta.data[colon@meta.data$orig.ident == gsub("\\.", "-", sample),]$Donor)
		for (celltype in rownames(cM)){
		    samplenames <- c(samplenames, sample)
		    cellnames <- c(cellnames, celltype)
		    donors <- c(donors, Donor)
		    values <- c(values, cM[celltype,sample])
		}
	}
	data <- data.frame(samplenames,cellnames,values, donors)
	data <- data[order(sapply(data$donors, function(x) which(x == c("B001", "B004", "F", "A001", "A002", "A008", "A010", "A014", "A015", "A018", "CRC1", "CRC2", "CRC3")))), ]

	colOrder <-  c("CD4+","Tregs","CD8+","Naive T","NK","ILCs","DC","Macrophages","Mast","Memory B","Naive B","GC","Plasma")

	data <- data[order(sapply(data$cellnames, function(x) which(x == colOrder))), ]
	data$cellnames <- factor(data$cellnames, levels = colOrder)
	counts <- table(colon@meta.data$CellType)[colOrder]
	axis_labels <- paste0(levels(data$cellnames), " (N = ", counts, ")")

	#data <- data[order(data$diseaseStates),]
	data$samplenames <- factor(data$samplenames, levels = unique(data$samplenames))
	pdf(paste0(paste("Sample-CellType-Bar-Fraction_Sample-ByDonor", sample_name, sep = "-"), ".pdf"), width = 8, height = 12, onefile=F)
	p <- ggplot(data, aes(fill=samplenames, y=values, x=cellnames)) + 
	 geom_bar(position="stack", stat="identity") + theme_ArchR() + 
	ylab("Fraction Cells") + scale_x_discrete(labels= axis_labels)+
	xlab("Sample") + theme(axis.text.x = element_text(angle = 90)) + theme(axis.text=element_text(size=20,hjust=0.95,vjust=0.2),
	    axis.title=element_text(size=20)) +
	scale_fill_manual("legend", values = c(
			    colorRampPalette(c("#007849","#1E392A"))(4), 
			    colorRampPalette(c("#008F95", "#0375B4","#062F4F"))(4),
			    colorRampPalette(c("#94618E", "#6E3667","#813772"))(2),
			    colorRampPalette(c("#E24E42", "#fe3401","#B82601"))(11),
			    colorRampPalette(c("#fbd2b6", "#F47D2B","#913f08"))(13),
			    colorRampPalette(c("#fff7b3", "#FEE500"))(2),
			    colorRampPalette(c("#dae1f1", "#8A9FD1"))(2),
			    colorRampPalette(c("#e0b8d6", "#C06CAB","#6b2e5c"))(10),
			    colorRampPalette(c("#b1e7df", "#3BBCA8","#1f6157"))(13),
			    colorRampPalette(c("#c1e8f0", "#90D5E4"))(2),
			    colorRampPalette(c("#F37B7D"))(1),
			    colorRampPalette(c("#9983BD"))(1),
			    colorRampPalette(c("#D24B27"))(1)))
	print(p)
	dev.off()

	##############################################################################################################################
	##############################################################################################################################
	# Boxplots and Wilcoxin tests for differential abundance
	cM <- t(as.matrix(confusionMatrix(paste0(colon@meta.data$CellType), paste0(colon@meta.data$orig.ident))))
	#combine and remove replicate samples
	cM["A002-C-010",] <- cM["A002-C-010-R0",]+cM["A002-C-010",]
	cM["A002-C-121",] <- cM["A002-C-121-R0",]+cM["A002-C-121",]
	cM <- cM[rownames(cM) %ni% c("A002-C-010-R0", "A002-C-121-R0"),]

	cM <- cM[rowSums(cM)>50,] #must have >50 cells
	cM <- cM/rowSums(cM)

	for (cellTypeCompared in c("MemoryB","CD4pos","GC","NaiveB","CD8pos","NaiveT","Macrophages","Mast","Plasma","Tregs","NK","DC","ILCs")){
	  dataTemp <- data[,c(cellTypeCompared, "DiseaseState")]
	  colnames(dataTemp) <- c("cellTypeCompared", "DiseaseState")
	  wilcox <- compare_means(cellTypeCompared ~ DiseaseState, ref.group = "Unaffected", p.adjust.method = "bonferroni", method='wilcox.test', data = dataTemp)
	  pdf(paste0(paste("cell-type-boxplot-padj-vs-unaffected-", cellTypeCompared, sep = "-"), ".pdf"), width = 4, height = 4, onefile=F)
	  p <- ggplot(data, aes_string(x = "DiseaseState", y = cellTypeCompared)) +
	      geom_boxplot(position = position_dodge(), outlier.shape = NA) + geom_jitter(color="black", size=2, alpha=0.9) + stat_pvalue_manual(wilcox, label = "p.adj", y.position = max(data[,cellTypeCompared])*c(1.1, 1.2,1.3))+ theme_ArchR()
	  print(p)
	  dev.off()
	}
}


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 10) Module scores for IFN-gamma signaling
# https://www.gsea-msigdb.org/gsea/msigdb/cards/REACTOME_INTERFERON_GAMMA_SIGNALING
if (10 %in% execute_steps){
	REACTOME_INTERFERON_GAMMA_SIGNALING <- c("B2M",
		"CAMK2A",
		"CAMK2B",
		"CAMK2D",
		"CAMK2G",
		"CD44",
		"CIITA",
		"FCGR1A",
		"FCGR1B",
		"GBP1",
		"GBP2",
		"GBP3",
		"GBP4",
		"GBP5",
		"GBP6",
		"GBP7",
		"HLA-A",
		"HLA-B",
		"HLA-C",
		"HLA-DPA1",
		"HLA-DPB1",
		"HLA-DQA1",
		"HLA-DQA2",
		"HLA-DQB1",
		"HLA-DQB2",
		"HLA-DRA",
		"HLA-DRB1",
		"HLA-DRB3",
		"HLA-DRB4",
		"HLA-DRB5",
		"HLA-E",
		"HLA-F",
		"HLA-G",
		"HLA-H",
		"ICAM1",
		"IFI30",
		"IFNG",
		"IFNGR1",
		"IFNGR2",
		"IRF1",
		"IRF2",
		"IRF3",
		"IRF4",
		"IRF5",
		"IRF6",
		"IRF7",
		"IRF8",
		"IRF9",
		"JAK1",
		"JAK2",
		"MID1",
		"MT2A",
		"NCAM1",
		"OAS1",
		"OAS2",
		"OAS3",
		"OASL",
		"PIAS1",
		"PML",
		"PRKCD",
		"PTAFR",
		"PTPN1",
		"PTPN11",
		"PTPN2",
		"PTPN6",
		"SOCS1",
		"SOCS3",
		"SP100",
		"STAT1",
		"SUMO1",
		"TRIM10",
		"TRIM14",
		"TRIM17",
		"TRIM2",
		"TRIM21",
		"TRIM22",
		"TRIM25",
		"TRIM26",
		"TRIM29",
		"TRIM3",
		"TRIM31",
		"TRIM34",
		"TRIM35",
		"TRIM38",
		"TRIM45",
		"TRIM46",
		"TRIM48",
		"TRIM5",
		"TRIM6",
		"TRIM62",
		"TRIM68",
		"TRIM8",
		"VCAM1")
	colon <- AddModuleScore(
	object = colon,
	features = list(REACTOME_INTERFERON_GAMMA_SIGNALING),
	ctrl = 50,
	name = 'REACTOME_INTERFERON_GAMMA_SIGNALING'
	)

	#https://www.gsea-msigdb.org/gsea/msigdb/cards/ST_INTERFERON_GAMMA_PATHWAY
	ST_INTERFERON_GAMMA_PATHWAY <- c("CISH",
			"ELP2",
			"IFNG",
			"IFNGR1",
			"JAK1",
			"JAK2",
			"PLA2G2A",
			"PTPRU",
			"REG1A",
			"STAT1")
	colon <- AddModuleScore(
	object = colon,
	features = list(ST_INTERFERON_GAMMA_PATHWAY),
	ctrl = 50,
	name = 'ST_INTERFERON_GAMMA_PATHWAY'
	)

	pal <- c("#D51F26")
	names(pal) <- "CD4+"
	pal["CD8+"] <- "#89288F"
	pal["DC"] <- "#F47D2B"
	pal["GC"] <- "#8A9FD1"
	pal["ILCs"] <- "#C06CAB"
	pal["Macrophages"] <- "#D8A767"
	pal["Mast"] <- "#89C75F"
	pal["Memory B"] <- "#F37B7D"
	pal["Naive B"] <- "#9983BD"
	pal["NK"] <- "#D24B27"
	pal["Plasma"] <- "#3BBCA8"
	pal["Tregs"] <- "#6E4B9E"
	pal["Naive T"] <- "#0C727C" 

	pdf('REACTOME_INTERFERON_GAMMA_SIGNALING.pdf', width = 20, onefile=F)
	print(VlnPlot(colon, features = 'REACTOME_INTERFERON_GAMMA_SIGNALING1', group.by = "CellType", pt.size = 0, cols = pal)+theme_ArchR())
	dev.off()

	pdf('REACTOME_INTERFERON_GAMMA_SIGNALING_sample.pdf', width = 20, onefile=F)
	print(VlnPlot(colon, features = 'REACTOME_INTERFERON_GAMMA_SIGNALING1', group.by = "orig.ident", pt.size = 0, cols = c(rep("#D51F26",100)))+theme_ArchR())
	dev.off()

	pdf('REACTOME_INTERFERON_GAMMA_SIGNALING_diseasestate.pdf', width = 20, onefile=F)
	print(VlnPlot(colon, features = 'REACTOME_INTERFERON_GAMMA_SIGNALING1', group.by = "DiseaseState", pt.size = 0, cols = c("#fe3401", "#0375B4", "#6E3667", "#1E392A"))+theme_ArchR())
	dev.off()


	pdf('ST_INTERFERON_GAMMA_PATHWAY_celltype.pdf', width = 20, onefile=F)
	print(VlnPlot(colon, features = 'ST_INTERFERON_GAMMA_PATHWAY1', group.by = "CellType", pt.size = 0, cols = pal)+theme_ArchR())
	dev.off()

	pdf('ST_INTERFERON_GAMMA_PATHWAY_sample.pdf', width = 20, onefile=F)
	print(VlnPlot(colon, features = 'ST_INTERFERON_GAMMA_PATHWAY1', group.by = "orig.ident", pt.size = 0, cols = c(rep("#D51F26",100)))+theme_ArchR())
	dev.off()

	pdf('ST_INTERFERON_GAMMA_PATHWAY_diseasestate.pdf', width = 20, onefile=F)
	print(VlnPlot(colon, features = 'ST_INTERFERON_GAMMA_PATHWAY1', group.by = "DiseaseState", pt.size = 0, cols = c("#fe3401", "#0375B4", "#6E3667", "#1E392A"))+theme_ArchR())
	dev.off()
}
