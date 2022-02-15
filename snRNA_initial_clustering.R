# Script to analyze scRNA data from scHTAN and scHuBMAP Projects
# For this script we will just make some qc plots and do simple clustering to divide the samples into groups
# WRB 2020

# Import libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ArchR)
library(viridis)
library(DoubletFinder)
set.seed(1)
# These are the analysis steps tht will be performed below, set execute steps to reflect which ones you want to do, note some are dependent on previous steps
# 1) Create seurat objects for individual 10x runs, run doblet finder and filter most likely doublets, merge into a seurat object containing all samples
# 2) QC and filtering
# 3) Add metadata
# 4) Normalize and scale data
# 5) Plot UMAPs
# 6) ID Cluster Markers
# 7) Initial cluster identification and seperation into compartments

execute_steps <- c(1,2,3,4,5,6,7)

# Define variables
sample_name <- "all_samples"
individual_qc_and_dublet_plot_location <- "./initial_clustering/doublet_analysis/"
analysis_parent_folder <- "./initial_clustering/"
setwd(analysis_parent_folder)
path_to_metadata <- "./hubmap_htan_metadata_atac_and_rna_final.csv"


# Define sets and locations of files for initial processing
scRNA_data_path <- "./data/"
scRNA_set_names <- c("all_samples_s12a", "all_samples_s12b", "all_samples_s5", "all_samples_s6", "all_samples_s7", "all_samples_s8", "all_samples_s3", "all_samples_s4", "all_samples_s9", "all_samples_s10", "all_samples_s11")
scRNA_sets <- list(c("EP007","A014-C-201","A002-C-010","A001-C-207","A001-C-124"), 
			c("B004-A-204", "B004-A-104", "A002-C-106", "A001-C-223", "A002-C-204", "A014-C-111"),
			c("A015-C-203","A015-C-204","A014-C-040","A002-C-201","A002-C-203"), 
			c("B001-A-301","B001-A-401","A015-C-008","A015-C-208","A014-C-114","A014-C-043","A001-C-202", "A001-C-119"),
			c("A001-C-108","A002-C-021","A002-C-212","A002-C-205","A014-C-101","A014-C-108", "A001-C-104"),
			c("A015-C-005","A015-C-006","A015-C-106","A002-C-114","A015-C-104","A015-C-202"), 
			c("A001-C-014","A001-C-023","A002-C-016","A002-C-024","A014-C-001","A014-C-054","A015-C-002","A015-C-010"),
			c("A015-C-109","B004-A-004","A002-C-121-R0", "A002-C-010-R0"), 
			c("A001-C-007", "A001-C-203", "A002-C-116", "A002-C-121", "A008-E-008", "A008-E-015", "A010-E-018", "A010-E-023", "A014-C-008", "A014-C-052"),
			c("A015-C-001", "A018-E-013", "A018-E-020", "B001-A-406", "B001-A-501", "B004-A-008", "EP034", "EP072B", "EP091", "A022-E-022"),
			c("CRC1_8810", "CRC2_15564", "CRC3_11773")) 

# Define functions
seurat_standard_normalize_and_scale <- function(colon, cluster, cluster_resolution){
	# colon is seurat object, 
	colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
	colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(colon)
	colon <- ScaleData(colon, features = all.genes)
	colon <- RunPCA(colon, features = VariableFeatures(object = colon))
	if (cluster){
		colon <- FindNeighbors(colon, dims = 1:20)
		colon <- FindClusters(colon, resolution = cluster_resolution)
	}
	colon <- RunUMAP(colon, dims = 1:20)
	return(colon)
}

make_seurat_object_and_doublet_removal <- function(data_directory, project_name){
	# function for basic seurat based qc and doubletfinder based doublet removal
	colon.data <- Read10X(data.dir = data_directory)
	currentSample <- CreateSeuratObject(counts = colon.data, project = project_name, min.cells = 3, min.features = 40)
	currentSample[["percent.mt"]] <- PercentageFeatureSet(currentSample, pattern = "^MT-")

	# qc plot-pre filtering
	pdf(paste0("./qc_plots_", project_name, "_prefiltered.pdf"))
	print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05))
	dev.off()
	pdf(paste0("./qc_plots_", project_name, "_prefiltered_no_points.pdf"))
	print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
	dev.off()

	# filter everything to 400 unique genes/cell
	currentSample <- subset(currentSample, subset = nFeature_RNA > 400 & nFeature_RNA < 4000)
	
	# Normalize and make UMAP
	currentSample <- seurat_standard_normalize_and_scale(currentSample, FALSE)

	# Run doublet finder
	nExp_poi <- round(0.08*length(currentSample@meta.data$orig.ident)*length(currentSample@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
	seu_colon <- doubletFinder_v3(currentSample, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	print(head(seu_colon@meta.data))
	
	# rename columns
	seu_colon$doublet.class <- seu_colon[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]]
	seu_colon[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]] <- NULL
	pann <- grep(pattern="^pANN", x=names(seu_colon@meta.data), value=TRUE)
	seu_colon$pANN <- seu_colon[[pann]]
	seu_colon[[pann]] <- NULL

	# plot pre and post doublet finder results
	pdf(paste0("./UMAP_pre_double_removal", project_name, ".pdf"))
	print(DimPlot(seu_colon, reduction = "umap", group.by = "doublet.class", cols = c("#D51F26", "#272E6A")))
	dev.off()
	seu_colon <- subset(seu_colon, subset = doublet.class != "Doublet")
	pdf(paste0("./UMAP_post_double_removal", project_name, ".pdf"))
	print(DimPlot(seu_colon, reduction = "umap", cols = c("#D51F26")))
	dev.off()

	# Remove extra stuff and return filtered Seurat object
	seu_colon <- DietSeurat(seu_colon, counts=TRUE, data=TRUE, scale.data=FALSE, assays="RNA")
	return(seu_colon)
}

seurat_qc_plots <- function(colon, sample_name){
	# Make some basic qc plots
	pdf(paste0("./seurat_nFeature_plots_", sample_name, ".pdf"), width = 40, height = 15)
	print(VlnPlot(colon, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.2))
	dev.off()

	pdf(paste0("./seurat_nCount_plots_", sample_name, ".pdf"), width = 40, height = 15)
	print(VlnPlot(colon, features = c("nCount_RNA"), ncol = 1, pt.size = 0.2))
	dev.off()

	pdf(paste0("./seurat_pMT_plots_", sample_name, ".pdf"), width = 40, height = 15)
	print(VlnPlot(colon, features = c("percent.mt"), ncol = 1, pt.size = 0.2))
	dev.off()
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 1) Create seurat objects for individual 10x runs, run doblet finder and filter most likely doublets, merge into a seurat object containing all samples
if (1 %in% execute_steps){
	setwd(individual_qc_and_dublet_plot_location)
	for (j in 1:length(scRNA_set_names)){
		samples <- scRNA_sets[[j]]
		print(paste0(scRNA_data_path, samples[1], "/outs/filtered_feature_bc_matrix/"))
		data_directory <- paste0(scRNA_data_path, samples[1], "/outs/filtered_feature_bc_matrix/")
		sample1 <- make_seurat_object_and_doublet_removal(data_directory, samples[1])
		seu_list <- c()
		for (i in 2:length(samples)){
			data_directory <- paste0(scRNA_data_path, samples[i], "/outs/filtered_feature_bc_matrix/")
			seu_list <- c(seu_list, make_seurat_object_and_doublet_removal(data_directory, samples[i]))
		}
		current_merge <- merge(sample1, y = seu_list, add.cell.ids = samples, project = scRNA_set_names[j])
		if (j==1){
			colon <- current_merge
		} else if (j>1){
			colon <- merge(colon, y = current_merge, project = "full_colon_project")
		}
	}
	setwd(analysis_parent_folder)
	colon[["percent.mt"]] <- PercentageFeatureSet(colon, pattern = "^MT-")
	saveRDS(colon, "initial_full_colon_proj_seurat.rds")
} elif (2 %in% execute_steps) {
	# if 1 is not a step, but 2 is then you need to load the initial object
	setwd(analysis_parent_folder)
	colon <- readRDS("initial_full_colon_proj_seurat.rds")
	# clean up if necessary:
	pann <- grep(pattern="^pANN", x=names(colon@meta.data), value=TRUE)
	for (i in pann){
	  colon[[i]] <- NULL
	}
	pann <- grep(pattern="^DF.classification", x=names(colon@meta.data), value=TRUE)
	for (i in pann){
	  colon[[i]] <- NULL
	}
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 2) QC
if (2 %in% execute_steps){
	# create and set working directory to save qc plots
	if (!dir.exists(paste0(analysis_parent_folder, "all_samples_qc_plots"))){
		dir.create(paste0(analysis_parent_folder, "all_samples_qc_plots"))
	}
	setwd(paste0(analysis_parent_folder, "all_samples_qc_plots"))

	# make the standard seurat qc plots
	seurat_qc_plots(colon, sample_name)

	# Now subset the project (if not done already)
	colon <- subset(colon, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 10000)

	# leave the qc directory
	setwd(analysis_parent_folder)
}


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 3) Add metadata
if (3 %in% execute_steps){
	metadata <- read.table(path_to_metadata, header = TRUE, sep = ",", stringsAsFactors=FALSE)
	
	# remove atac column
	metadata <- metadata[,colnames(metadata)[2:28]]
	colnames(metadata) <- c("Sample", colnames(metadata)[2:27])
	metadata <- metadata[metadata$Sample != "",]

	meta_data_types <- colnames(metadata)
	for (i in 2:length(meta_data_types)){
		identities <- colon[['orig.ident']]
		for (j in 1:length(metadata$Sample)){
			identities[identities==metadata$Sample[j]] <- metadata[j,meta_data_types[i]]
		}
		colon <- AddMetaData(colon, identities$orig.ident, col.name = meta_data_types[i])
	}
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 4) Normalize and scale data
if (4 %in% execute_steps){
	colon <- seurat_standard_normalize_and_scale(colon, TRUE, 1.0)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 5) Plot UMAPs
if (5 %in% execute_steps){
	# Plot by clustering, sample, and disease state
	colon <- FindClusters(colon, resolution = 0.5)
	pdf(paste0("./UMAPclustering" , ".pdf"), onefile=F)
	print(DimPlot(colon, reduction = "umap", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAP_samples.pdf"), width = 12, onefile=F)
	print(DimPlot(colon, reduction = "umap", group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAP_donor.pdf"), width = 12, onefile=F)
	print(DimPlot(colon, reduction = "umap", group.by = "Donor", cols = paletteDiscrete(values = unique(colon@meta.data$Donor), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAP_disease_state.pdf"), width = 6.5, onefile=F)
	print(DimPlot(colon, reduction = "umap", group.by = "DiseaseState",
		cols = c("#D51F26", "#D7CEC7", "#89288F", "#D7CEC7")) + theme_ArchR())
	dev.off()

	saveRDS(colon, "clustered_full_colon_proj_seurat.rds")
}


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 6) ID Cluster Markers
if (6 %in% execute_steps){
	# Find cluster specific markers
	colon.markers <- FindAllMarkers(colon, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, max.cells.per.ident = 250)
	top10 <- colon.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
}


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 7) Initial cluster identification
if (7 %in% execute_steps){
	# Rough initial cluster identification
	new.cluster.ids <- c("Epithelial", #0
		"Epithelial",
		"Epithelial",
		"Epithelial",
		"Epithelial",
		"Epithelial",#5
		"Epithelial",
		"Epithelial",
		"Immune",
		"Stromal",
		"Immune",#10
		"Epithelial",
		"Stromal",
		"Epithelial",
		"Immune",
		"Epithelial",#15
		"Stromal",
		"Stromal",
		"Immune",
		"Tuft",
		"Epithelial",#20
		"Immune",
		"Epithelial",
		"Enteroendocrine",
		"Stromal",
		"Immune"#25
		)

	identities <- as.character(colon@meta.data$seurat_clusters)
	for (i in 0:length(new.cluster.ids)){
		identities[identities==as.character(i)] <- new.cluster.ids[i+1]
	}
	colon <- AddMetaData(colon, identities, col.name = "CellTypeInitial")
	pdf(paste0("./UMAP_cell_type_initial.pdf"), width = 12)
	AugmentPlot(DimPlot(colon, reduction = "umap", group.by = "CellTypeInitial"))
	dev.off()

	pdf(paste0("./UMAP_cell_type_initial_v2.pdf"), width = 12)
	plot = DimPlot(colon, reduction = "umap", group.by = "CellTypeInitial")
	plot = LabelClusters(plot = plot, id = "CellTypeInitial")
	AugmentPlot(plot)
	dev.off()

	# Save ids of immune, epithelial, and stromal cells
	immune_clusters <- c(8,10,14,18,21,25)
	stromal_clusters <- c(9,12,16,17,24)
	epithelial_clusters <- c(0:7,11,13,15,19,20,22,23)
	epithelial_specialized_clusters <- c(19,23)
	epithelial_nonspecialized_clusters <- c(0:7,12,14,20,21,25)

	# Subset to make immune, stromal, and epithelial projects
	colon_immune <- DietSeurat(subset(colon, subset = seurat_clusters %in% immune_clusters))
	colon_stromal <- DietSeurat(subset(colon, subset = seurat_clusters %in% stromal_clusters))
	colon_epitehlial <- DietSeurat(subset(colon, subset = seurat_clusters %in% epithelial_clusters))
	colon_specialized_epitehlial <- DietSeurat(subset(colon, subset = seurat_clusters %in% epithelial_specialized_clusters))
	colon_nonspecialized_epitehlial <- DietSeurat(subset(colon, subset = seurat_clusters %in% epithelial_nonspecialized_clusters))

	saveRDS(colon_immune, file = "colon_immune_all_samples_initial.rds")
	saveRDS(colon_stromal, file = "colon_stromal_all_samples_initial.rds")
	saveRDS(colon_epitehlial, file = "colon_epithelial_all_samples_initial.rds")
	saveRDS(colon_specialized_epitehlial, file = "colon_specialized_epitehlial_all_samples_initial.rds")
	saveRDS(colon_nonspecialized_epitehlial, file = "colon_nonspecialized_epitehlial_all_samples_initial.rds")
}



