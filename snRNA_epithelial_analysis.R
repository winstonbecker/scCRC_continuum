# Script to analyze scRNA data for epithelail cells from scHTAN colon project
# See publication "Single-cell analyses reveal a continuum of cell state and composition changes 
# in the malignant transformation of polyps to colorectal cancer."
# WRB 2020-2021
# Substantial portions of the LSI dimensionality reduction and projections were adapted or copied from Granja et al 2019

# Analysis steps
# 1) Load epithelial project and subset
# 2) QC and filtering
# 3) LSI Dim reduction
# 4) ID Cluster Markers
# 5) Plot Markers
# 6) Define clusters
# 7) Project all samples into normal manifold
# 8) Compute differential tests
# 9) Compute malignancy continuum
# 10) Cluster genes for heatmap

###############################################################################################################################
###############################################################################################################################
# Import libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(ArchR)
library(viridis)
library(DoubletFinder)
library(Rcpp)
library(harmony)
library(future)
library(Matrix)
set.seed(1)

###############################################################################################################################
###############################################################################################################################
# Define variables
execute_steps <- 1:10

# folder to save the results
analysis_parent_folder <- "./epithelial_results/"

# load seurat object containing all epithelial cells
colon_full <- readRDS("./initial_clustering/colon_epithelial_all_samples.rds")

# load metadata
metadata <- read.table("./hubmap_htan_metadata_atac_and_rna_final.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)

# LSI Params
resolution <- c(0.1,0.2,0.4,0.8)
umapNeighbors <- 50
umapMinDist <- 0.45
umapDistMetric <- "cosine"
nTop <- 1600
initial_nPCs <- 1:8
final_nPCs <- c(1:4,6:8)
nPCs <- initial_nPCs

# Create directory
if (!dir.exists(paste0(analysis_parent_folder))){
	dir.create(paste0(analysis_parent_folder))
}
setwd(analysis_parent_folder)

###############################################################################################################################
###############################################################################################################################
# Define Functions--as noted above, multiple functions copied from Granja et al 2019
sourceCpp(code='
  #include <Rcpp.h>
  using namespace Rcpp;
  using namespace std;
  // [[Rcpp::export]]
  Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n) {
    const int nv = j.size();
    const int nm = rm.size();
    Rcpp::NumericVector rv(nm);
    Rcpp::NumericVector rit(nm);
    int current;
    // Calculate RowVars Initial
    for (int i = 0; i < nv; ++i) {
      current = j(i) - 1;
      rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
      rit(current) = rit(current) + 1;
    }
    // Calculate Remainder Variance
    for (int i = 0; i < nm; ++i) {
      rv(i) = rv(i) + (n - rit(i))*rm(i)*rm(i);
    }
    rv = rv / (n - 1);
    return(rv);
  }'
)

sparseRowVariances <- function (m){
    #Compute Fast Sparse Row Variances--From Granja et al 2019
    rM <- Matrix::rowMeans(m)
    rV <- computeSparseRowVariances(m@i + 1, m@x, rM, ncol(m))
    return(rV)
}

getVarGenesFilterBlacklist <- function(mat, nTopGenes = 2000, blacklist = NULL){
  # Get the top nTopGenes variable genes present in mat (a gene x sample/cell matrix)
  # If blacklist is present, do not return any genes in the blacklist
  if(!is.null(blacklist)){
    mat <- mat[!rownames(mat) %in% blacklist,]
  }
  # compute row variances and return the most variable genes for eith matrix or sparse matrix
  if (class(mat) == "matrix"){
    rownames(mat)[head(order(rowVars(mat), decreasing = TRUE), nTopGenes)]
  } else {
    rownames(mat)[head(order(sparseRowVariances(mat), decreasing = TRUE), nTopGenes)]
  }
}

calcLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL){
	#From Granja et al 2019
    set.seed(1)

    #TF IDF LSI adapted from flyATAC
    if(binarize){
        message(paste0("Binarizing matrix..."))
        mat@x[mat@x > 0] <- 1
    }

    #Calc RowSums and ColSums
    colSm <- Matrix::colSums(mat)
    rowSm <- Matrix::rowSums(mat)

    #Calc TF IDF
    message("Computing Term Frequency IDF...")
    freqs <- t(t(mat)/colSm)
    idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs

    #Calc SVD then LSI
    message("Computing SVD using irlba...")
    svd <- irlba::irlba(tfidf, nComponents, nComponents)
    svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
    diag(svdDiag) <- svd$d
    matSVD <- t(svdDiag %*% t(svd$v))
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))

    #Return Object
    out <- list(
        matSVD = matSVD, 
        rowSm = rowSm, 
        colSm = colSm, 
        svd = svd, 
        binarize = binarize, 
        nComponents = nComponents,
        seed = 1)

    out

}

#Helper function for summing sparse matrix groups
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

projectLSI <- function(mat, lsi, binarize){   
    
    #Get Same Features
    mat <- mat[lsi$varFeatures,]
    if(binarize){
        message(paste0("Binarizing matrix..."))
        mat@x[mat@x > 0] <- 1       
    }
    
    #Calc TF IDF
    rowsToZero <- which(lsi$rowSm == 0)
    setToZero <- which((mat@i + 1) %in% rowsToZero)
    if(length(setToZero) > 0){
        mat@x[setToZero] <- 0
    }

    message("Computing Term Frequency IDF...")
    freqs <- t(t(mat)/Matrix::colSums(mat))
    idf   <- as(log(1 + length(lsi$colSm) / lsi$rowSm), "sparseVector")
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
    if(length(Matrix::which(is.na(tfidf),arr.ind=TRUE)) > 0){
        tfidf[Matrix::which(is.na(tfidf),arr.ind=TRUE)] <- 0 #weird Inf * 0
    }

    #Calc V
    V <- t(tfidf) %*% lsi$svd$u %*% diag(1/lsi$svd$d)

    #Calc SVD then LSI
    message("Computing SVD using irlba...")
    svdDiag <- matrix(0, nrow=lsi$nComponents, ncol=lsi$nComponents)
    diag(svdDiag) <- lsi$svd$d
    matSVD <- t(svdDiag %*% t(V))
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
    
    return(matSVD)
}

seurat_feature_plot <- function(sample_name, reduction, cell_type, markers){
	# function to make grid layout of seurat feature plots
	p1 <- FeaturePlot(colon, features = markers, reduction = reduction, sort.cell = TRUE, combine = FALSE, pt.size = 2)
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
	print(CombinePlots(lapply(p1, function (x) AugmentPlot(x))))
	dev.off()
}

vln_plot <- function(features, save_name){
	pdf(save_name, width = 20, onefile=F)
	print(VlnPlot(colon, features = features, group.by = "orig.ident", pt.size = 0, cols = c(rep("#D51F26",100)))+
	geom_boxplot(outlier.shape = NA, alpha = 0.6)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))+
	scale_x_discrete(labels=paste0(data.frame(table(colon@meta.data$orig.ident))$Var1, "\n n = ",  data.frame(table(colon@meta.data$orig.ident))$Freq)))
	dev.off()
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 1) Subset to create normal project

if (1 %in% execute_steps){
	# define normal project containing all normal samples except B001-A-104, which was much lower quality and was not used in reference as a result
	colon <- DietSeurat(subset(colon_full, subset = DiseaseState == "Normal"))
	colon <- DietSeurat(subset(colon, subset = orig.ident != "B004-A-104"))
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

	# already done, but can add here if you got here a different way
	colon[["percent.mt"]] <- PercentageFeatureSet(colon, pattern = "^MT-")

	# Now subset the project and make some nice qc plots for just the epithelial cells--this was already done in the first step
	colon <- subset(colon, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 10000)
	vln_plot(c("nFeature_RNA"), paste0("./n_genes_violin_", sample_name, ".pdf"))
	vln_plot(c("nCount_RNA"), paste0("./n_counts_violin_", sample_name, ".pdf"))
	vln_plot(c("percent.mt"), paste0("./pMT_violin_", sample_name, ".pdf"))

	setwd(analysis_parent_folder)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 3) Dim reduction with LSI
if (3 %in% execute_steps){
	#Initialize list for storing iterative LSI output
	lsiOut <- list()
	nPCs <- initial_nPCs

	# Get the raw counts
	rawCounts <- GetAssayData(object = colon, slot = "counts")

	# Identify genes (mitochondrial, ribosomal, and HLA) to blacklist for dimensionality reduction
	mt_genes <- grep("^MT-", rownames(colon), value = TRUE)
	rps_genes <- grep("^RPS", rownames(colon), value = TRUE)
	rpl_genes <- grep("^RPL", rownames(colon), value = TRUE)
	hla_genes <- grep("^HLA-", rownames(colon), value = TRUE)
	blacklist <- c(mt_genes, rps_genes, rpl_genes, hla_genes)

	# Here we will use the normalized data in the data slot from seurat, which is ln rather than log2 normalized
	log2CP10k <- GetAssayData(object = colon, slot = "data")
	# alternatively could do:
	# log2CP10k <- as(log2(t(t(rawCounts)/Matrix::colSums(rawCounts)) * 10000 + 1), "matrix")

	for(i in seq_along(resolution)){
	    if(i == 1){
	    	#Initial LSI uses variances that are across all single cells and will have larger batch relationships
	        message("Running initial LSI...")
	        varGenes <- getVarGenesFilterBlacklist(log2CP10k, nTopGenes = nTop, blacklist = blacklist)
	    }else{
	    	message(sprintf("Running LSI %s of %s...", i,  length(resolution)))
	        # Calculate variable genes using clusters defined in previous LSI iteration
	        clusterMat <- edgeR::cpm(groupSums(rawCounts, clusters, sparse = TRUE), log=TRUE, prior.count = 3)
	        varGenes <- getVarGenesFilterBlacklist(clusterMat, nTopGenes = nTop, blacklist = blacklist)
	    }
	    
	    # Now run LSI and find clusters
	    LSIi <- calcLSI(rawCounts[varGenes,], nComponents = max(nPCs), binarize = FALSE)
	    colon[[paste0("LSI_iter",i)]] <- CreateDimReducObject(embeddings = LSIi$matSVD, key = sprintf("LSI%s_", i), assay = "RNA")
	    colon <- FindNeighbors(object = colon, reduction = paste0("LSI_iter",i), dims = nPCs, force.recalc = TRUE)
	    colon <- FindClusters(object = colon, resolution = resolution[i])
	    clusters <- Idents(colon)
	    
	    # Save LSI iteration
	    lsiOut[[paste0("LSI_iter",i)]] <- list(
	        lsiMat = LSIi$matSVD, 
	        varFeatures = varGenes, 
	        clusters = clusters,
	        colSm = LSIi$colSm,
	        rowSm = LSIi$rowSm,
	        svd = LSIi$svd, 
	        binarize = LSIi$binarize, 
	        nComponents = LSIi$nComponents
	    )
	}

	# Run uwot to compute the UMAP
	nPCs <- final_nPCs
	uwotUmap <- uwot::umap(
	    LSIi$matSVD[,nPCs], 
	    n_neighbors = umapNeighbors, 
	    min_dist = umapMinDist, 
	    metric = umapDistMetric, 
	    n_threads = 1, 
	    verbose = TRUE, 
	    ret_nn = TRUE,
	    ret_model = TRUE
	    )

	# Add to the seurat object and plot
	umap_vals <- uwotUmap[[1]][,1:2]
	rownames(umap_vals) <- rownames(LSIi$matSVD[,nPCs])
	colon[["uwot_UMAP"]] <- CreateDimReducObject(embeddings = umap_vals, key = "uwot_UMAP", assay = "RNA")

	pdf(paste0("./uwot_UMAP_samples_", nTop, "_vargenes_all_", max(nPCs), "PCs.pdf"), width = 8)
	print(DimPlot(colon, reduction = "uwot_UMAP", group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))#, cols = (ArchRPalettes$stallion))
	dev.off()

    # Recluster with final nPCs
    colon <- FindNeighbors(object = colon, reduction = paste0("LSI_iter",i), dims = nPCs, force.recalc = TRUE)
    colon <- FindClusters(object = colon, resolution = 1.0)

    pdf(paste0("./UMAPclustering" , ".pdf"))
    print(DimPlot(colon, reduction = "uwot_UMAP", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
    dev.off()
}


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 4) ID Cluster Markers to help with manual annotation
if (4 %in% execute_steps){
	colon.markers <- FindAllMarkers(colon, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, max.cells.per.ident = 250)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 5) Plot Some Markers on the UMAPs to help with annotation
if (5 %in% execute_steps){
	reductions_to_plot <- c("uwot_UMAP")
	for (reduction in reductions_to_plot){
		seurat_feature_plot(sample_name, reduction, "CyclingTA", c("TICRR","CDC25C"))
		seurat_feature_plot(sample_name, reduction, "ImmatureEnterocytes", c("SLC26A2","CA1"))
		seurat_feature_plot(sample_name, reduction, "Tuft", c("GNG13","SH2D7","SH2D6","TRPM5","AZGP1","KRT18","BMX","PSTPIP2","LRMP","PTGS1","IL17RB","HCK","PLCG2","ANXA13"))
		seurat_feature_plot(sample_name, reduction, "Best4posEnterocytes", c("BEST4", "CA7","OTOP2","OTOP3", "MYOM1","MT1G","MT1H"))
		seurat_feature_plot(sample_name, reduction, "General_Epithelial", c("EPCAM", "KRT8","KRT18"))
		seurat_feature_plot(sample_name, reduction, "Immature_Goblet", c("KLK1","ITLN1","WFDC2","CLCA1","LRRC26","RETNLB","SPINK4","AGR2"))
		seurat_feature_plot(sample_name, reduction, "Goblet", c("MUC2", "TFF1", "FCGBP","FFAR4","SYTL2","LGALS9B","BCAS1"))
		seurat_feature_plot(sample_name, reduction, "Stem", c("SMOC2", "RGMB", "LGR5", "ASCL2", "SOX9", "CD34"))
		seurat_feature_plot(sample_name, reduction, "Enteroendocrine", c("CRYBA2","SCGN","FEV","CHGA","GCG","SCG5","PCSK1N","PYY","NEUROD1","MS4A8","DDC"))
	}
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 6) Define clusters
if (6 %in% execute_steps){
	# cluster identification
	new.cluster.ids <- c("TA2", #0
        "TA1", #1
        "Enterocytes", #2
        "Enterocyte Progenitors", #3
        "Immature Goblet", #4
        "Immature Enterocytes", #5
        "TA2", #6
        "Immature Goblet", #7
        "Enterocytes", #8
        "TA1", #9
        "Stem", #10
        "Goblet", #11
        "CyclingTA", #12
        "Immature Enterocytes", #13
        "Best4+ Enterocytes", #14
        "Best4+ Enterocytes", #15
        "Goblet", #16
        "Enteroendocrine", #17
        "Tuft")#18

	# Set the cell types in the project
	identities <- data.frame(colon[['seurat_clusters']])
	identities$seurat_clusters <- as.character(colon[['seurat_clusters']]$seurat_clusters)
	for (i in 0:(length(new.cluster.ids)-1)){
		identities[identities$seurat_clusters==as.character(i),] <- new.cluster.ids[i+1]
	}

	# plot
	colon <- AddMetaData(colon, identities$seurat_clusters, col.name = "CellType")
	pdf(paste0("./UMAP_cell_type_id.pdf"), width = 6)
	DimPlot(colon, group.by = "CellType", reduction = "uwot_UMAP", cols = paletteDiscrete(values = unique(colon@meta.data$CellType), set = "stallion", reverse = FALSE)) + theme_ArchR()
	dev.off()
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 7) Project all samples into normal manifold
if (7 %in% execute_steps){
	# Subset data for a sample
	cell_type_df <- data.frame()
	for (sampleName in unique(colon_full@meta.data$orig.ident)){
	    # Subset seurat project to get sample only and get the normalized data
	    temp <- subset(colon_full, subset = orig.ident == sampleName)
		rawCounts <- GetAssayData(object = temp, slot = "counts")

	    # project in the previously defined lsi dimensionality reduction and add to temp seurat project--assumes lsiOut is still around, save it and load it if its not
	    lsiProjectionMat <- projectLSI(rawCounts,lsiOut[[paste0("LSI_iter",length(resolution))]], binarize = FALSE)
	    temp[["LSI_project"]] <- CreateDimReducObject(embeddings = as.matrix(lsiProjectionMat), key = sprintf("LSI_project"), assay = "RNA")

	    # project into umap and add to temp seurat project
	    umapProjection <- uwot::umap_transform(as.matrix(Embeddings(temp, reduction = "LSI_project"))[,nPCs], uwotUmap, verbose = TRUE)
	    proDF <- data.frame(X1 = umapProjection[,1], X2 = umapProjection[,2])
	    rownames(proDF) <- rownames(Embeddings(temp, reduction = "LSI_project"))
	    temp[["uwot_UMAP_projection"]] <- CreateDimReducObject(embeddings = as.matrix(proDF), key = "uwot_UMAP_projection", assay = "RNA")

	    # now plot the projection and the reference
	    refDF <- data.frame(X1 = uwotUmap$embedding[,1], X2 = uwotUmap$embedding[,2], Type = "Reference")
	    proDF <- data.frame(X1 = umapProjection[,1], X2 = umapProjection[,2], Type = "Sample")
	    projectionDF <- rbind(refDF, proDF)
	    p <- ggplot(projectionDF, aes(X1,X2,color=Type)) + 
	      geom_point(size = 0.1, alpha = 1) + 
	      xlab("UMAP Dimension 1") + 
	      ylab("UMAP Dimension 2") + 
	      theme_ArchR(baseSize = 10) + ggtitle("projection") + 
	      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
	      scale_color_manual(values=c("Reference"="#D5D5D5","Sample"="#008F95"))
	    pdf(paste0("./uwot_projection_", sampleName, "_", nTop, "_vargenes_all_", max(nPCs), "PCs.pdf"), width = 12)
	    print(p)
	    dev.off()

	    # Now id cell based on nearest neghbors
	    library(FNN)
	    input_knn <- 25
	    nPCs <- final_nPCs
	    svdReference <- as.data.frame(lsiOut[[paste0("LSI_iter",length(resolution))]]$lsiMat)
	    svdDisease <- as.data.frame(as.matrix(lsiProjectionMat))
	    knnDisease <- get.knnx(
	      data = svdReference[,nPCs],
	      query = svdDisease[,nPCs],
	      k = input_knn)

	    cellTypes <- c()
	    for (j in 1:length(knnDisease$nn.index[,1])){
	        types <- colon@meta.data$CellType[knnDisease$nn.index[j,]]
	        cellTypes <- c(cellTypes,labels(sort(table(types),decreasing=TRUE)[1]))
	    }
	    refDF$CellType <- "Reference"
	    proDF$CellType <- cellTypes
	    projectionDF2 <- rbind(refDF, proDF)

	    pal <- paletteDiscrete(values = unique(proDF$CellType), set = "stallion", reverse = FALSE)
	    pal["Reference"] <- "#D5D5D5"
	    p <- ggplot(projectionDF2, aes(X1,X2,color=CellType)) + 
	      geom_point(size = 0.2, alpha = 1) +
	      xlab("UMAP Dimension 1") + 
	      ylab("UMAP Dimension 2") +
	      theme_ArchR(baseSize = 10) + ggtitle(sampleName) + 
	      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
	      scale_color_manual(values=c("Reference"="#D5D5D5", "Best4+ Enterocytes"="#D51F26", "CyclingTA"="#272E6A", "Enterocyte Progenitors" = "#208A42", "Enterocytes" = "#89288F", 
	      	"Enteroendocrine" = "#F47D2B", "Goblet" = "#FEE500", "Immature Enterocytes" = "#8A9FD1", "Immature Goblet" = "#C06CAB", "Secretory TA" = "#D24B27", "Stem" = "#D8A767", "TA1" = "#90D5E4", "TA2" = "#89C75F", "Tuft" = "#F37B7D"))
	    pdf(paste0("./uwot_projection_", sampleName, "_", nTop, "_vargenes_all_", max(nPCs), "PCs.pdf"), width = 6)
	    print(p)
	    dev.off()

	    current_cell_types <- DataFrame("Cell" = rownames(svdDisease), "CellType" = cellTypes)
	    cell_type_df <- rbind(cell_type_df, current_cell_types)
	}
	write.csv(cell_type_df, "allEpithelialCellTypes.csv")
	colon_full <- AddMetaData(colon_full, cell_type_df$CellType, col.name = "CellType")

	# save seurat objects
	saveRDS(colon, "clustered_normal_colon_proj_seurat.rds")
	saveRDS(DietSeurat(colon), "diet_clustered_normal_colon_proj_seurat.rds")

	saveRDS(colon_full, "clustered_all_samples_epithelial_colon_proj_seurat.rds")
	saveRDS(DietSeurat(colon_full), "diet_clustered_all_samples_epithelial_colon_proj_seurat.rds")
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 8) Compute differential tests of stem cells relative to normal and unaffected stem cells
if (8 %in% execute_steps){
	colon_full <- readRDS("diet_clustered_all_samples_epithelial_colon_proj_seurat.rds")
	backgrounds <- c("Unaffected", "NormalColon")
	
	# first set up a colon test object where the idents include test and background groups
	for (background in backgrounds){
		if (background == "NormalColon"){
			sample_name_list <- paste0(unique(colon_full@meta.data$SimplifiedSampleName), "Stem")
			sample_name_list <- sample_name_list[sample_name_list != "NormalColonStem"]
			colon_full <- AddMetaData(colon_full, paste0(colon_full@meta.data$SimplifiedSampleName, colon_full@meta.data$CellType), col.name = "CellTypeSimplifiedSample")
			colon_test_object <- colon_full
			Idents(colon_test_object) <- "CellTypeSimplifiedSample"
			background_sample <- "NormalColonStem"
		}
		if (background == "Unaffected"){
			sample_name_list <- paste0(unique(colon_full@meta.data$DifferentialGroup), "Stem")
			sample_name_list <- sample_name_list[sample_name_list %ni% c("NormalColonStem", "UnaffectedStem")]

			simplified_sample_names <- colon_full@meta.data$DifferentialGroup
			colon_full <- AddMetaData(colon_full, paste0(simplified_sample_names, colon_full@meta.data$CellType), col.name = "CellTypeSimplifiedSample2")

			colon_test_object <- colon_full
			Idents(colon_test_object) <- "CellTypeSimplifiedSample2"
			background_sample <- "UnaffectedStem"
		}

		# identify samples with at least 100 cells
		new_sample_list <- c()
		for (i in 1:length(sample_name_list)){
			current_sample <- sample_name_list[i]
			if (sum(Idents(colon_test_object) == current_sample)>100){
				new_sample_list <- c(new_sample_list, current_sample)
			}
		}
		sample_name_list <- new_sample_list

		# compute differential genes using seurats findMarkers, merge into single object
		for (i in 1:length(sample_name_list)){
			current_sample <- sample_name_list[i]
			message(paste0("Computing differential genes for ", current_sample))
			# added use.test = "DESeq2" and got similar results
      			differential_test <- FindMarkers(colon_test_object, ident.1 = current_sample, ident.2 = background_sample, verbose = FALSE, min.pct = 0, logfc.threshold = 0, min.cells.feature = 0, max.cells.per.ident = 300, test.use = "MAST")
			colnames(differential_test) <- paste0(colnames(differential_test), current_sample)
			if (i == 1){
				colon_full_diff_test <- differential_test
			} else {
				message(paste0("Merging differential expression data for ", current_sample))
				colon_full_diff_test <- merge(colon_full_diff_test, differential_test, by=0, all=TRUE)
				rownames(colon_full_diff_test) <- colon_full_diff_test$Row.names
				colon_full_diff_test <- colon_full_diff_test[,colnames(colon_full_diff_test) %ni% c("Row.names")]
			}
		}
		saveRDS(colon_full_diff_test, paste0("./colon_full_diff_test_", background_sample, "_background.rds"))
	}
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 9) Compute malignancy continuum
if (9 %in% execute_steps){
	# load the differential test relative to normal stem cells
	colon_full_diff_test <- readRDS("./colon_full_diff_test_NormalColonStem_background.rds")
	background <- "Normal"
	typeA <- "Stem"

	# select just the logfc and p values
	colon_full_diff_test_avg_logFC <- colon_full_diff_test[,grepl("avg_logFC", colnames(colon_full_diff_test))]
	colon_full_diff_test_avg_logFC <- colon_full_diff_test_avg_logFC/0.69314718056
	colon_full_diff_test_p_val_adj <- colon_full_diff_test[,grepl("p_val_adj", colnames(colon_full_diff_test))]
	colon_full_diff_test_avg_logFC[is.na(colon_full_diff_test_avg_logFC)] <- 0
	colon_full_diff_test_p_val_adj[is.na(colon_full_diff_test_p_val_adj)] <- 1
	
	# select significant genes with the following cutoffs
	colon_full_diff_test_avg_logFC_significant <- colon_full_diff_test_avg_logFC[(rowSums(colon_full_diff_test_avg_logFC>0.5 & colon_full_diff_test_p_val_adj<0.05)>1 | rowSums(colon_full_diff_test_avg_logFC<(-0.5) & colon_full_diff_test_p_val_adj<0.05)>1), ]
	colon_full_diff_test_p_val_adj_significant <- colon_full_diff_test_p_val_adj[(rowSums(colon_full_diff_test_avg_logFC>0.5 & colon_full_diff_test_p_val_adj<0.05)>1 | rowSums(colon_full_diff_test_avg_logFC<(-0.5) & colon_full_diff_test_p_val_adj<0.05)>1), ]

	# compute the pcs on the logfold changes
	pcs <- prcomp(colon_full_diff_test_avg_logFC_significant)
	pc_df <- data.frame(pcs$rotation)
	pc_df$sample <- rownames(pc_df)

	# remove atac column from metadata
	metadata <- metadata[,colnames(metadata)[2:28]]
	colnames(metadata) <- c("Sample", colnames(metadata)[2:27])
	metadata <- metadata[metadata$Sample != "",]
	sample_names <- substr(rownames(pc_df),10,nchar(rownames(pc_df))-(nchar(typeA)))
	metadata <- metadata[metadata$SimplifiedSampleName %in% sample_names,]
	metadata <- metadata[!duplicated(metadata),]
	# remove specific partial duplicates "A002-C-010-R0", "A002-C-121-R0", "A014-C-114"
	metadata <- metadata[metadata$Sample %ni% c("A002-C-010-R0", "A002-C-121-R0", "A014-C-114"),]
	rownames(metadata) <- metadata$SimplifiedSampleName
	metadata <- metadata[sample_names,]

	# combine the pcs and the metadata
	pc_df <- cbind(pc_df, metadata)
	scalef <- 1
	pc_df$PC1 <- pc_df$PC1*scalef

	# plot the first two pcs
	p <- ggplot(pc_df, aes(x=PC2, y=PC1, color=DiseaseState)) +
	geom_point(size=2) + scale_color_manual(values=c("#D51F26", "#89288F", "#208A42")) + theme_ArchR()
	ggsave(plot=p,height=5,width=5, filename=paste0("All_Polyp_", typeA, "_vs_", background, "_", typeA, "_pca_on_diff_genes.pdf"), useDingbats=FALSE)

	# fit a spline (0.11 knot selected based on previous plot)
	require(splines)
	fit<-lm(PC1 ~ bs(PC2,knots = c(0.11)),data = pc_df )
	age.grid<-seq(from=-3000*scalef, to = 4500*scalef)/10000
	splinefit <- data.frame(age.grid,predict(fit,newdata = list(PC2=age.grid)))
	colnames(splinefit) <- c("PC2", "PC1")

	# plot the pcs with the spline fit
	p <- ggplot(pc_df, aes(x=PC2, y=PC1, color=DiseaseState)) +
	geom_point(size=2) + scale_color_manual(values=c("#D51F26", "#89288F", "#208A42")) + 
	xlab(paste0("PC2 (",100*summary(pcs)$importance["Proportion of Variance","PC2"], "%)")) + ylab(paste0("PC1 (",100*summary(pcs)$importance["Proportion of Variance","PC1"], "%)")) + theme_ArchR()
	p <- p+ geom_line(data=splinefit, colour="#CC0000") 
	ggsave(plot=p,height=5,width=5, filename=paste0("All_Polyp_", typeA, "_vs_", background, "_", typeA, "_pca_on_diff_genes.pdf"), useDingbats=FALSE)

	#get the nearest spline point for each sample
	x_vals <- c()
	for (i in 1:length(pc_df[,1])){
		x_vals <- c(x_vals,splinefit$PC2[which.min(sqrt((pc_df$PC2[i]-splinefit$PC2)**2+(pc_df$PC1[i]-splinefit$PC1)**2))])
	}
	pc_df$nearest_spline_x_vals <- x_vals
	pc_df <- pc_df[order(pc_df$nearest_spline_x_vals),]

	# save the results
	saveRDS(pc_df, paste0("pc_df_", background, "_background_", typeA, "_cellType.rds"))
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 10) Cluster genes for heatmap
if (10 %in% execute_steps){
	colon_full_diff_test <- readRDS("./colon_full_diff_test_UnaffectedStem_background.rds")#"./colon_full_diff_test.rds")
	background <- "Unaffected"
	typeA <- "Stem"

	colon_full_diff_test_avg_logFC <- colon_full_diff_test[,grepl("avg_logFC", colnames(colon_full_diff_test))]
	colon_full_diff_test_avg_logFC <- colon_full_diff_test_avg_logFC/0.69314718056
	colon_full_diff_test_p_val_adj <- colon_full_diff_test[,grepl("p_val_adj", colnames(colon_full_diff_test))]
	colon_full_diff_test_avg_logFC[is.na(colon_full_diff_test_avg_logFC)] <- 0
	colon_full_diff_test_p_val_adj[is.na(colon_full_diff_test_p_val_adj)] <- 1
	# 1/28/2020
	colon_full_diff_test_avg_logFC_significant <- colon_full_diff_test_avg_logFC[(rowSums(colon_full_diff_test_avg_logFC>=0.75 & colon_full_diff_test_p_val_adj<=0.05)>1 | rowSums(colon_full_diff_test_avg_logFC<=(-0.75) & colon_full_diff_test_p_val_adj<=0.05)>1), ]
	colon_full_diff_test_p_val_adj_significant <- colon_full_diff_test_p_val_adj[(rowSums(colon_full_diff_test_avg_logFC>=0.75 & colon_full_diff_test_p_val_adj<=0.05)>1 | rowSums(colon_full_diff_test_avg_logFC<=(-0.75) & colon_full_diff_test_p_val_adj<=0.05)>1), ]
	order_samples <- rownames(pc_df)[rownames(pc_df) %in% colnames(colon_full_diff_test_avg_logFC_significant)]

	set.seed(1)
	gene_clusters <- kmeans(cbind(colon_full_diff_test_avg_logFC_significant), 10, iter.max = 500, algorithm = "Lloyd")
}
