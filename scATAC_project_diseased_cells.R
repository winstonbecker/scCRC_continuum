# Script to project diseased cells into normal manifold
# Adapted from Granja et al 2019, Nature Biotech and Granja et al 2021, Nature Genetics

# Load Packages
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

#Load Genome Annotations
addArchRGenome("hg38")

#Set Threads to be used
addArchRThreads()
`%notin%` <- Negate(`%in%`)

############################################################
############################################################
# User defined things to set/load
setwd("./projects/")

# Load full project
proj_all_sample_peaks <- loadArchRProject(path = "./all_samples_HuBMAP_HTAN/")

# Define non-epithelial clusters in full project (these will be the ones you dont project)
nonEpithelialClusters <- paste0("C", 28:35)

# Load normal colon reference ArchR project
proj_hubmap_normals <- loadArchRProject(path = "./all_samples_HuBMAP_normal_colon_epithelial_cells_current_doublets_removed/")

# Load reference LSI
lsi <- getReducedDims(proj_hubmap_normals, reducedDims = "IterativeLSIhubmap_colon_epithelial", returnMatrix = FALSE)

# Load reference UMAP manifold
umap <- getEmbedding(proj_hubmap_normals, embedding = "UMAPhubmap_colon_epithelial", returnDF = FALSE)
umapManifold <- uwot::load_uwot(umap$params$uwotModel[1])

# Define path containing all arrow files
full_project_arrow_file_folder <- "./projects/all_samples_HuBMAP_HTAN/ArrowFiles/"

# Define samples to project
samplesToProject <- c("A002-C-106-D",
  "A014-C-002-D",
  "A014-C-043-D",
  "A002-C-203-D",
  "A002-C-204-D",
  "A015-C-006-D",
  "A001-C-002-D",
  "A014-C-013-D",
  "A002-C-201-D",
  "A002-C-212-D",
  "A014-C-111-D",
  "A015-C-204-D",
  "A002-C-114-D",
  "A018-E-020-D",
  "A015-C-104-D",
  "A014-C-040-D",
  "A001-C-119-D",
  "A002-C-025-S2",
  "A001-C-207-D",
  "A015-C-102-D",
  "A014-C-101-D",
  "A002-C-010-D",
  "A015-C-203-D",
  "A002-C-205-D",
  "B001-A-301-D",
  "EP072B-D",
  "EP034-D",
  "A015-C-202-D",
  "A015-C-106-D",
  "A001-C-202-D",
  "A022-E-002-D",
  "A001-C-108-D",
  "A001-C-218-D",
  "A001-C-104-D-R2",
  "A001-C-223-D",
  "A001-C-124-D-R2",
  "A015-C-005-D",
  "EP007-D",
  "A014-C-114-D",
  "A001-C-203-D",
  "A014-C-008-D",
  "A018-E-013-D",
  "A015-C-208-D",
  "A014-C-001-D",
  "A002-C-202-D-OCT",
  "B004-A-008-D",
  "A014-C-108-D",
  "A002-C-116-S2",
  "A014-C-052-D",
  "A014-C-201-D",
  "A002-C-021-D",
  "A015-C-008-D",
  "A002-C-016-D",
  "A015-C-010-D-R2",
  "A002-C-010-D-R1",
  "A010-E-018-D",
  "A001-C-124-D-R1",
  "A001-C-007-D",
  "A010-E-023-D",
  "A002-C-025-D",
  "A002-C-121-D",
  "B004-A-004-D-R2",
  "A015-C-109-D",
  "A008-E-015-D",
  "A002-C-116-D",
  "A015-C-206-D",
  "A015-C-002-D",
  "B001-A-401-D",
  "A002-C-024-D",
  "A002-C-010-S2",
  "A008-E-008-D",
  "A001-C-014-D",
  "B004-A-204-D",
  "A002-C-121-S2",
  "B001-A-406-D",
  "A014-C-054-D",
  "A001-C-123-D",
  "B001-A-302-D",
  "B001-A-501-D",
  "A001-C-023-D-R2",
  "A001-C-023-D-R1",
  "B004-A-004-D",
  "A001-C-104-D-R1",
  "A015-C-001-D",
  "EP091-D",
  "CRC-1-8810-D",
  "CRC-2-15564-D",
  "CRC-3-11773-D",
  "CRC-4-8456-D")

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Define Functions

projectLSI <- function(
  # Function from Granja et al.
  mat = NULL, 
  LSI = NULL, 
  returnModel = FALSE, 
  verbose = FALSE
  ){   
    
  out2 <- tryCatch({
   
    require(Matrix)
    set.seed(LSI$seed)
    
    #Get Same Features
    mat <- mat[LSI$idx,]

    #Binarize Matrix
    if(LSI$binarize){
        mat@x[mat@x > 0] <- 1       
    }
    
    #TF
    colSm <- Matrix::colSums(mat)
    if(any(colSm == 0)){
      exclude <- which(colSm==0)
      mat <- mat[,-exclude]
      colSm <- colSm[-exclude]
    }
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))

    if(LSI$LSIMethod == 1 | tolower(LSI$LSIMethod) == "tf-logidf"){

      #Adapted from Casanovich et al.

      #LogIDF
      idf   <- as(log(1 + LSI$nCol / LSI$rowSm), "sparseVector")

      #TF-LogIDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

    }else if(LSI$LSIMethod == 2 | tolower(LSI$LSIMethod) == "log(tf-idf)"){

      #Adapted from Stuart et al.

      #IDF
      idf   <- as(LSI$nCol / LSI$rowSm, "sparseVector")

      #TF-IDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

      #Log transform TF-IDF
      mat@x <- log(mat@x * LSI$scaleTo + 1)  

    }else if(LSI$LSIMethod == 3 | tolower(LSI$LSIMethod) == "logtf-logidf"){

      #LogTF
      mat@x <- log(mat@x + 1)

      #LogIDF
      idf   <- as(log(1 + LSI$nCol / LSI$rowSm), "sparseVector")

      #TF-IDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat


    }else{

      stop("LSIMethod unrecognized please select valid method!")

    }

    gc()

    #Clean Up Matrix
    idxNA <- Matrix::which(is.na(mat),arr.ind=TRUE)
    if(length(idxNA) > 0){
        mat[idxNA] <- 0
    }

    #Calc V
    V <- Matrix::t(mat) %*% LSI$svd$u %*% Matrix::diag(1/LSI$svd$d)

    #LSI Diagonal
    svdDiag <- matrix(0, nrow=LSI$nDimensions, ncol=LSI$nDimensions)
    diag(svdDiag) <- LSI$svd$d
    matSVD <- Matrix::t(svdDiag %*% Matrix::t(V))
    matSVD <- as.matrix(matSVD)
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))

    if(returnModel){
        X <- LSI$svd$u %*% diag(LSI$svd$d) %*% t(V)
        out <- list(matSVD = matSVD, V = V, X = X)
    }else{
        out <- matSVD
    }

    out

  }, error = function(e){

    errorList <- list(
      mat = mat,
      colSm = if(exists("colSm", inherits = FALSE)) colSm else "Error with colSm!",
      idf = if(exists("idf", inherits = FALSE)) idf else "Error with idf!",
      V = if(exists("V", inherits = FALSE)) V else "Error with V!",
      matSVD = if(exists("matSVD", inherits = FALSE)) matSVD else "Error with matSVD!"
    )
  })

  out2

}


############################################################
# Run script to project cells
############################################################
cell_type_df <- data.frame("Cell" = c(), "CellType" = c()) # initialize to save projected identities
gene_annotations <- getArchRGenome(geneAnnotation = TRUE) # load annotations

# iterate through samples, projecting cells
for (sampleName in samplesToProject){
  idxSample <- BiocGenerics::which(getCellColData(proj_all_sample_peaks, "Clusters") %notin% nonEpithelialClusters & getCellColData(proj_all_sample_peaks, "Sample") %in% c(sampleName))
  cellsSample <- proj_all_sample_peaks$cellNames[idxSample[["Clusters"]]]

  mat_se <- getMatrixFromArrow(
    ArrowFile = paste(paste(full_project_arrow_file_folder, sampleName, sep = ""), ".arrow", sep = ""),
    useMatrix = "TileMatrix",
    useSeqnames = NULL,
    cellNames = cellsSample[grep(sampleName,cellsSample)],
    ArchRProj = proj_all_sample_peaks,
    verbose = TRUE,
    binarize = TRUE
  )

  subset_rows <- paste(rowData(mat_se)$seqnames, rowData(mat_se)$start) %in% paste(lsi$LSIFeatures$seqnames, lsi$LSIFeatures$start)
  mat <- assay(mat_se)
  mat <- mat[subset_rows,]
  matSVD <- projectLSI(mat, lsi)
  lsiProjection <- projectLSI(mat, lsi)

  set.seed(1)
  umapProjection <- uwot::umap_transform(as.matrix(matSVD)[,1:25], umapManifold, verbose = TRUE)

  #Plot Projection
  refDF <- data.frame(X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2], Type = "Reference")
  proDF <- data.frame(X1 = umapProjection[,1], X2 = umapProjection[,2], Type = "Sample")
  projectionDF <- rbind(refDF, proDF)

  plotDir <- paste0("./")

  #Input Parameters
  input_knn <- 25
  scaleTo <- 10000
  nMax <- 500

  #LSI-SVD
  svdReference <- as.data.frame(lsi$matSVD)
  svdDisease <- as.data.frame(as.matrix(lsiProjection))

  #Differential Seed
  set.seed(1)

  #Cells that we are testing of disease
  idxDisease <- cellsSample

  #If the number of cells is greater than 5 continue
  stopifnot(length(idxDisease) > 5)

  #KNN Nearest Neighbor using FNN #find 25 nn cells
  library(FNN)
  knnDisease <- get.knnx(
      data = svdReference,
      query = svdDisease,#[idxDisease, ], #Subset by idxDisease 
      k = input_knn)

  #Determine the minimum KNN where reference cells are less than 1.25x disease cells
  i <- 0
  uniqueIdx <- unique(as.vector(knnDisease$nn.index))
  while(length(uniqueIdx) > 1.25 * length(idxDisease)){
      i <- i + 1
      uniqueIdx <- unique(as.vector(knnDisease$nn.index[,seq_len(input_knn-i)]))
  }

  cellTypes <- c()
  for (j in 1:length(knnDisease$nn.index[,1])){
    types <- getCellColData(proj_hubmap_normals, "CellType")[knnDisease$nn.index[j,],]
    cellTypes <- c(cellTypes,labels(sort(table(types),decreasing=TRUE)[1]))
  }

  current_cell_types <- DataFrame("Cell" = rownames(lsiProjection), "CellType" = cellTypes)
  write.csv(current_cell_types, paste0(sampleName,"-CellType.csv"))
  cell_type_df <- rbind(cell_type_df, current_cell_types)

  refDF$CellType <- "Reference"
  proDF$CellType <- cellTypes
  projectionDF2 <- rbind(refDF, proDF)

  p <- ggplot(projectionDF2, aes(X1,X2,color=CellType)) + 
      geom_point(size = 0.1, alpha = 1) +
      xlab("UMAP Dimension 1") + 
      ylab("UMAP Dimension 2") +
      theme_ArchR(baseSize = 10) + ggtitle(sampleName) + 
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
      scale_color_manual(values=c("Reference"="#D5D5D5", "Best4+ Enterocytes"="#D51F26", "Enterocyte Progenitors" = "#208A42", "Enterocytes" = "#89288F", 
        "Enteroendocrine" = "#F47D2B", "Goblet" = "#FEE500", "Immature Enterocytes" = "#8A9FD1", "Immature Goblet" = "#C06CAB", "Secretory TA" = "#D24B27", "Stem" = "#D8A767", "TA1" = "#90D5E4", "TA2" = "#89C75F"))
  plotPDF(p, name = paste0(sampleName,"-Projection-UMAP2.pdf"), ArchRProj = proj_all_sample_peaks, addDOC = FALSE, width = 5, height = 5)
}

write.csv(cell_type_df, paste0("all_samples_cell_types_projections.csv"))


