# scATAC analysis for HTAN FAP project
# WRB 2020-2021

# Script to run ArchR on scHTAN/HuBMAP data

# The following modules should be loaded:
# ml R/3.6.1

# The following is done in this script:
# 1) Make arrow files and project, bounds set to isolate population of cells
# 2) Create project
# 3) Add metadata to project
# 4) QC Plots
# 5) Initial Dimensionality reductions and clustering (Updated on 5/17/2020 to make nicer UMAP using 100K cells instead of 10K)
# 6) Marker Genes plotting and analysis

# Load libraries
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

# Set/Create Working Directory to Folder
setwd("./projects/")
`%notin%` <- Negate(`%in%`)

# Load Genome Annotations
addArchRGenome("hg38")

# Set Threads to be used
# threads <- 16
addArchRThreads()

# Which steps above you want to run
execute_steps <- c(1,2,3,4,5,6)


metadata <- read.table("./hubmap_htan_metadata_atac_and_rna_final.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)

pre_defined_doublets <- TRUE
cells_to_include <- read.table("./cells_in_initial_clustering.txt")


################################################################################################################
#..............................................................................................................#
################################################################################################################
# 1) Make arrow files

if (1 in execute_steps){
  # # Input Files
  # # switched 010-D and 025-D because there was clearly an issue
  #   "A002-C-024-S2" = "/oak/stanford/groups/wjg/wbecker/other/scATAC/20200226_nextSeq/A002_S2/outs/fragments.tsv.gz",
  inputFiles1p5K <- c(
    "A001-C-218-D" = "./scATAC/fragments_files/HTAN/A001-C-218-D_20200214_fragments.tsv.gz",
    "A002-C-025-S2" = "./scATAC/fragments_files/HTAN/A002-C-025-S2_20200310_fragments.tsv.gz")


  inputFiles2K <- c(
    "A002-C-016-D" = "./scATAC/fragments_files/HTAN/A002-C-016-D_20200715_fragments.tsv.gz",
    "A014-C-001-D" = "./scATAC/fragments_files/HTAN/A014-C-001-D_20200715_fragments.tsv.gz",
    "A015-C-002-D" = "./scATAC/fragments_files/HTAN/A015-C-002-D_20200715_fragments.tsv.gz",
    "A015-C-010-D-R2" = "./scATAC/fragments_files/HTAN/A015-C-010-D_20200715_fragments.tsv.gz",
    "A015-C-204-D" = "./scATAC/fragments_files/HTAN/A015-C-204-D_20200804_fragments.tsv.gz"
  )

  # # switched 010-D and 025-D because there was clearly an issue
  inputFiles3K <- c(
    "A001-C-014-D" = "./scATAC/fragments_files/HTAN/A001-C-014-D_20200715_fragments.tsv.gz",
    "A001-C-023-D-R1" = "./scATAC/fragments_files/HTAN/A001-C-023-D_20200214_fragments.tsv.gz",
    "A001-C-023-D-R2" = "./scATAC/fragments_files/HTAN/A001-C-023-D_20200715_fragments.tsv.gz",
    "A001-C-104-D-R1" = "./scATAC/fragments_files/HTAN/A001-C-104-D_20200214_fragments.tsv.gz",
    "A001-C-104-D-R2" = "./scATAC/fragments_files/HTAN/A001-C-104-D_20200811_fragments.tsv.gz",
    "A001-C-108-D" = "./scATAC/fragments_files/HTAN/A001-C-108-D_20200811_fragments.tsv.gz",
    "A001-C-119-D" = "./scATAC/fragments_files/HTAN/A001-C-119-D_20200811_fragments.tsv.gz",
    "A001-C-123-D" = "./scATAC/fragments_files/HTAN/A001-C-123-D_20200214_fragments.tsv.gz",
    "A001-C-203-D" = "./scATAC/fragments_files/HTAN/A001-C-203-D_20200804_fragments.tsv.gz",
    "A002-C-025-D" = "./scATAC/fragments_files/HTAN/A002-C-025-D_20200310_fragments.tsv.gz",
    "A002-C-010-S2" = "./scATAC/fragments_files/HTAN/A002-C-010-S2_20200310_fragments.tsv.gz",
    "A002-C-021-D" = "./scATAC/fragments_files/HTAN/A002-C-021-D_20200811_fragments.tsv.gz",
    "A002-C-024-D" = "./scATAC/fragments_files/HTAN/A002-C-024-D_20200715_fragments.tsv.gz",
    "A002-C-010-D-R1" = "./scATAC/fragments_files/HTAN/A002-C-010-D_20200310_fragments.tsv.gz",
    "A002-C-116-D" = "./scATAC/fragments_files/HTAN/A002-C-116-D_20200310_fragments.tsv.gz",
    "A002-C-116-S2" = "./scATAC/fragments_files/HTAN/A002-C-116-S2_20200310_fragments.tsv.gz",
    "A002-C-121-D" = "./scATAC/fragments_files/HTAN/A002-C-121-D_20200310_fragments.tsv.gz",
    "A002-C-121-S2" = "./scATAC/fragments_files/HTAN/A002-C-121-S2_20200310_fragments.tsv.gz",
    "A002-C-201-D" = "./scATAC/fragments_files/HTAN/A002-C-201-D_08042020_fragments.tsv.gz",
    "A002-C-202-D-OCT" = "./scATAC/fragments_files/HTAN/A002-C-202-D-OCT_20200214_fragments.tsv.gz",
    "A002-C-205-D" = "./scATAC/fragments_files/HTAN/A002-C-205-D_20200811_fragments.tsv.gz",
    "A002-C-212-D" = "./scATAC/fragments_files/HTAN/A002-C-212-D_20200811_fragments.tsv.gz",
    "A014-C-002-D" = "./scATAC/fragments_files/HTAN/A014-C-002-D_20200811_fragments.tsv.gz",
    "A014-C-013-D" = "./scATAC/fragments_files/HTAN/A014-C-013-D_20200804_fragments.tsv.gz",
    "A014-C-043-D" = "./scATAC/fragments_files/HTAN/A014-C-043-D_08042020_fragments.tsv.gz",
    "A014-C-054-D" = "./scATAC/fragments_files/HTAN/A014-C-054-D_20200715_fragments.tsv.gz",
    "A014-C-101-D" = "./scATAC/fragments_files/HTAN/A014-C-101-D_20200811_fragments.tsv.gz",
    "A014-C-108-D" = "./scATAC/fragments_files/HTAN/A014-C-108-D_20200811_fragments.tsv.gz",
    "A014-C-114-D" = "./scATAC/fragments_files/HTAN/A014-C-114-D_20200804_fragments.tsv.gz",
    "A015-C-005-D" = "./scATAC/fragments_files/HTAN/A015-C-005-D_20200811_fragments.tsv.gz",
    "A015-C-008-D" = "./scATAC/fragments_files/HTAN/A015-C-008-D_20200804_fragments.tsv.gz",
    "A015-C-102-D" = "./scATAC/fragments_files/HTAN/A015-C-102-D_20200811_fragments.tsv.gz",
    "A015-C-104-D" = "./scATAC/fragments_files/HTAN/A015-C-104-D_20200811_fragments.tsv.gz",
    "A015-C-106-D" = "./scATAC/fragments_files/HTAN/A015-C-106-D_20200811_fragments.tsv.gz",
    "A015-C-109-D" = "./scATAC/fragments_files/HTAN/A015-C-109-D_20200715_fragments.tsv.gz",
    "A015-C-203-D" = "./scATAC/fragments_files/HTAN/A015-C-203-D_20200804_fragments.tsv.gz",
    "A015-C-206-D" = "./scATAC/fragments_files/HTAN/A015-C-206-D_20200715_fragments.tsv.gz",
    "A015-C-208-D" = "./scATAC/fragments_files/HTAN/A015-C-208-D_20200804_fragments.tsv.gz",
    "F007-D" = "./scATAC/fragments_files/HTAN/F007-D_20200702_fragments.tsv.gz",
    "A001-C-007-D" = "./scATAC/fragments_files/HTAN/A001-C-007-D_20200817_fragments.tsv.gz",
    "A014-C-052-D" = "./scATAC/fragments_files/HTAN/A014-C-052-D_20200817_fragments.tsv.gz",
    "A014-C-008-D" = "./scATAC/fragments_files/HTAN/A014-C-008-D_20200817_fragments.tsv.gz",
    "A015-C-001-D" = "./scATAC/fragments_files/HTAN/A015-C-001-D_20200817_fragments.tsv.gz",
    "F034-D" = "./scATAC/fragments_files/HTAN/F034-D_20200817_fragments.tsv.gz",
    "F091-D" = "./scATAC/fragments_files/HTAN/F091-D_20200817_fragments.tsv.gz",
    "CRC-1-8810-D" = "./scATAC/fragments_files/HTAN/CRC-1-8810-D_20200917_fragments.tsv.gz",
    "CRC-2-15564-D" = "./scATAC/fragments_files/HTAN/CRC-2-15564-D_20200917_fragments.tsv.gz",
    "CRC-3-11773-D" = "./scATAC/fragments_files/HTAN/CRC-3-11773-D_20200917_fragments.tsv.gz",
    "CRC-4-8456-D" = "./scATAC/fragments_files/HTAN/CRC-4-8456-D_20200917_fragments.tsv.gz"
  )

  inputFiles4K <- c(
    "A001-C-207-D" = "./scATAC/fragments_files/HTAN/A001-C-207-D_20200702_fragments.tsv.gz",
    "A001-C-223-D" = "./scATAC/fragments_files/HTAN/A001-C-223-D_20200702_fragments.tsv.gz",
    "A002-C-010-D" = "./scATAC/fragments_files/HTAN/A002-C-010-D_20200702_fragments.tsv.gz",
    "A002-C-114-D" = "./scATAC/fragments_files/HTAN/A002-C-114-D_20200811_fragments.tsv.gz",
    "A002-C-203-D" = "./scATAC/fragments_files/HTAN/A002-C-203-D_08042020_fragments.tsv.gz",
    "A014-C-040-D" = "./scATAC/fragments_files/HTAN/A014-C-040-D_20200804_fragments.tsv.gz",
    "A014-C-111-D" = "./scATAC/fragments_files/HTAN/A014-C-111-D_20200702_fragments.tsv.gz",
    "A014-C-201-D" = "./scATAC/fragments_files/HTAN/A014-C-201-D_20200702_fragments.tsv.gz",
    "A015-C-006-D" = "./scATAC/fragments_files/HTAN/A015-C-006-D_20200811_fragments.tsv.gz",
    "F072B-D" = "./scATAC/fragments_files/HTAN/F072B-D_20200817_fragments.tsv.gz",
    "A008-E-008-D" = "./scATAC/fragments_files/HTAN/A008-E-008-D_20200817_fragments.tsv.gz",
    "A010-E-018-D" = "./scATAC/fragments_files/HTAN/A010-E-018-D_20200817_fragments.tsv.gz",
    "A018-E-013-D" = "./scATAC/fragments_files/HTAN/A018-E-013-D_20200817_fragments.tsv.gz",
    "A015-C-202-D" = "./scATAC/fragments_files/HTAN/A015-C-202-D_20200817_fragments.tsv.gz",
    "A010-E-023-D" = "./scATAC/fragments_files/HTAN/A010-E-023-D_20200817_fragments.tsv.gz",
    "A018-E-020-D" = "./scATAC/fragments_files/HTAN/A018-E-020-D_20200817_fragments.tsv.gz",
    "A022-E-002-D" = "./scATAC/fragments_files/HTAN/A022-E-002-D_20200817_fragments.tsv.gz"
  )

  inputFiles5K <- c(
    "A001-C-002-D" = "./scATAC/fragments_files/HTAN/A001-C-002-D_20200702_fragments.tsv.gz",
    "A001-C-124-D-R1" = "./scATAC/fragments_files/HTAN/A001-C-124-D_20200214_fragments.tsv.gz",
    "A001-C-124-D-R2" = "./scATAC/fragments_files/HTAN/A001-C-124-D_20200702_fragments.tsv.gz",
    "A002-C-106-D" = "./scATAC/fragments_files/HTAN/A002-C-106-D_20200702_fragments.tsv.gz",
    "A002-C-204-D" = "./scATAC/fragments_files/HTAN/A002-C-204-D_20200702_fragments.tsv.gz",
    "A008-E-015-D" = "./scATAC/fragments_files/HTAN/A008-E-015-D_20200817_fragments.tsv.gz"
  )

  inputFiles10K <- c(
    "A001-C-202-D" = "./scATAC/fragments_files/HTAN/A001-C-202-D_20200804_fragments.tsv.gz"
  )

  #Create Arrow Files
  ArrowFiles1p5K <- createArrowFiles(
    inputFiles = inputFiles1p5K,
    sampleNames = names(inputFiles1p5K), filterFrags = 1500,
  )
  doubScores <- addDoubletScores(ArrowFiles1p5K, k = 10, knnMethod = "UMAP", LSIMethod = 1)

  ArrowFiles2K <- createArrowFiles(
    inputFiles = inputFiles2K,
    sampleNames = names(inputFiles2K), filterFrags = 2000,
  )

  ArrowFiles3K <- createArrowFiles(
    inputFiles = inputFiles3K,
    sampleNames = names(inputFiles3K), filterFrags = 3000,
  )

  ArrowFiles4K <- createArrowFiles(
    inputFiles = inputFiles4K,
    sampleNames = names(inputFiles4K), filterFrags = 4000,
  )

  ArrowFiles5K <- createArrowFiles(
    inputFiles = inputFiles5K,
    sampleNames = names(inputFiles5K), filterFrags = 5000,
  )

  ArrowFiles10K <- createArrowFiles(
    inputFiles = inputFiles10K,
    sampleNames = names(inputFiles10K), filterFrags = 10000,
  )

  doubScores <- addDoubletScores(ArrowFiles2K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
  doubScores <- addDoubletScores(ArrowFiles3K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
  doubScores <- addDoubletScores(ArrowFiles4K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
  doubScores <- addDoubletScores(ArrowFiles5K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
  doubScores <- addDoubletScores(ArrowFiles10K, k = 10, knnMethod = "UMAP", LSIMethod = 1)


  # HuBMAP Samples
  inputFiles1p5K <- c(
    "B001-A-302-D" = "./scATAC/fragments_files/HuBMAP/B001-A-302-D_fragments.tsv.gz",
    "B001-A-401-D" = "./scATAC/fragments_files/HuBMAP/B001-A-401-D_fragments.tsv.gz",
    "B001-A-406-D" = "./scATAC/fragments_files/HuBMAP/B001-A-406-D_fragments.tsv.gz",
    "B001-A-501-D" = "./scATAC/fragments_files/HuBMAP/B001-A-501-D_fragments.tsv.gz")

  inputFiles2K <- c(
    "B004-A-004-D" = "./scATAC/fragments_files/HuBMAP/B004-A-004-D_20200715_fragments.tsv.gz"
  )

  inputFiles3K <- c(
    "B001-A-301-D" = "./scATAC/fragments_files/HuBMAP/B001-A-301-D_20200804_fragments.tsv.gz",
    "B005-A-301-D" = "./scATAC/fragments_files/HuBMAP/B005-A-301-D_20200917_fragments.tsv.gz",
    "B005-A-501-D" = "./scATAC/fragments_files/HuBMAP/B005-A-501-D_20200917_fragments.tsv.gz"
  )

  inputFiles4K <- c(
    "B004-A-204-D" = "./scATAC/fragments_files/HuBMAP/B004-A-204-D_20200702_fragments.tsv.gz",
    "B004-A-004-D-R2" = "./scATAC/fragments_files/HuBMAP/B004-A-004-D_20200817_fragments.tsv.gz"
  )

  inputFiles5K <- c(
    "B004-A-008-D" = "./scATAC/fragments_files/HuBMAP/B004-A-008-D_20200817_fragments.tsv.gz"
  )

  ArrowFiles1p5K <- createArrowFiles(
    inputFiles = inputFiles1p5K,
    sampleNames = names(inputFiles1p5K), filterFrags = 1500,
  )

  ArrowFiles2K <- createArrowFiles(
    inputFiles = inputFiles2K,
    sampleNames = names(inputFiles2K), filterFrags = 2000,
  )

  ArrowFiles3K <- createArrowFiles(
    inputFiles = inputFiles3K,
    sampleNames = names(inputFiles3K), filterFrags = 3000,
  )

  ArrowFiles4K <- createArrowFiles(
    inputFiles = inputFiles4K,
    sampleNames = names(inputFiles4K), filterFrags = 4000,
  )

  ArrowFiles5K <- createArrowFiles(
    inputFiles = inputFiles5K,
    sampleNames = names(inputFiles5K), filterFrags = 5000,
  )

  doubScores <- addDoubletScores(ArrowFiles1p5K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
  doubScores <- addDoubletScores(ArrowFiles2K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
  doubScores <- addDoubletScores(ArrowFiles3K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
  doubScores <- addDoubletScores(ArrowFiles4K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
  doubScores <- addDoubletScores(ArrowFiles5K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
}


################################################################################################################
#..............................................................................................................#
################################################################################################################
# # 2) Create Project
if (2 in execute_steps){
  # Selecting only the samples we want included
  ArrowPath <- "./projects/ArrowFiles/"
  ArrowFiles <- c("A001-C-002-D.arrow",
    "A001-C-124-D-R1.arrow",
    "A002-C-010-S2.arrow",
    "A002-C-116-S2.arrow",
    "A002-C-212-D.arrow", 
    "A014-C-101-D.arrow",
    "A015-C-010-D-R2.arrow",
    "B004-A-004-D.arrow",
    "A001-C-014-D.arrow",
    "A001-C-124-D-R2.arrow",
    "A002-C-016-D.arrow",
    "A002-C-121-D.arrow",       
    "A014-C-108-D.arrow",
    "A015-C-102-D.arrow",
    "A001-C-023-D-R1.arrow",
    "A001-C-202-D.arrow",
    "A002-C-021-D.arrow",
    "A002-C-121-S2.arrow",      
    "A014-C-111-D.arrow",
    "A015-C-104-D.arrow",
    "B004-A-204-D.arrow",
    "A001-C-023-D-R2.arrow",
    "A001-C-203-D.arrow",
    "A002-C-024-D.arrow",
    "A002-C-201-D.arrow",       
    "A014-C-001-D.arrow",
    "A014-C-114-D.arrow",
    "A015-C-106-D.arrow",
    "A001-C-104-D-R1.arrow",
    "A001-C-207-D.arrow",
    "A002-C-025-D.arrow",
    "A002-C-202-D-OCT.arrow",   
    "A014-C-002-D.arrow",
    "A014-C-201-D.arrow",
    "A015-C-109-D.arrow",
    "B001-A-301-D.arrow",
    "A001-C-104-D-R2.arrow",
    "A001-C-218-D.arrow",
    "A002-C-025-S2.arrow",
    "A014-C-013-D.arrow",
    "A015-C-002-D.arrow",
    "A015-C-203-D.arrow",
    "B001-A-302-D.arrow",
    "A001-C-108-D.arrow",
    "A001-C-223-D.arrow",
    "A002-C-106-D.arrow",
    "A002-C-203-D.arrow",       
    "A014-C-040-D.arrow",
    "A015-C-005-D.arrow",
    "A015-C-204-D.arrow",
    "B001-A-401-D.arrow",
    "A001-C-119-D.arrow",
    "A002-C-010-D.arrow",
    "A002-C-114-D.arrow",
    "A002-C-204-D.arrow",       
    "A014-C-043-D.arrow",
    "A015-C-006-D.arrow",
    "A015-C-206-D.arrow",
    "B001-A-406-D.arrow",
    "F007-D.arrow",
    "A001-C-123-D.arrow",
    "A002-C-010-D-R1.arrow",
    "A002-C-116-D.arrow",
    "A002-C-205-D.arrow",       
    "A014-C-054-D.arrow",
    "A015-C-008-D.arrow",
    "A015-C-208-D.arrow",
    "B001-A-501-D.arrow",
    "A001-C-007-D.arrow",
    "A014-C-052-D.arrow",
    "A014-C-008-D.arrow",
    "A015-C-001-D.arrow",
    "F034-D.arrow",
    "F091-D.arrow",
    "F072B-D.arrow",
    "A008-E-008-D.arrow",
    "A010-E-018-D.arrow",
    "A018-E-013-D.arrow",
    "A015-C-202-D.arrow",
    "A010-E-023-D.arrow",
    "A018-E-020-D.arrow",
    "A022-E-002-D.arrow",
    "B004-A-004-D-R2.arrow",
    "A008-E-015-D.arrow",
    "B004-A-008-D.arrow",
    "CRC-1-8810-D.arrow",
    "CRC-2-15564-D.arrow",
    "CRC-3-11773-D.arrow",
    "CRC-4-8456-D.arrow")


  # #Create ArchRProject
  proj <- ArchRProject(
    ArrowFiles = paste0(ArrowPath, ArrowFiles), 
    outputDirectory = "all_samples_HuBMAP_HTAN"
  )

}

################################################################################################################
#..............................................................................................................#
################################################################################################################
# # 3) Add metadata to project
if (3 in execute_steps){
  
  for (j in 2:dim(metadata)[2]){
    # initialize list
    cellsNamesToAdd <- c()
    annotationToAdd <- c()
    for (i in 1:dim(metadata)[1]){
      idxSample <- BiocGenerics::which(getCellColData(proj, "Sample") %in% metadata[i,"Sample"])
      cellsSample <- proj$cellNames[idxSample[["Sample"]]]
      cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
      annotationToAdd <- append(annotationToAdd, rep(metadata[i,j], length(cellsSample)))
    }

    proj <- addCellColData(ArchRProj = proj, data = paste0(annotationToAdd), cells = paste0(cellsNamesToAdd), name = colnames(metadata)[j], force = TRUE)
  }
  saveArchRProject(ArchRProj = proj, outputDirectory = "all_samples_HuBMAP_HTAN", load = FALSE, overwrite = FALSE)
}


################################################################################################################
#..............................................................................................................#
################################################################################################################
# # 4) QC Plots and filter doublets
if (4 in execute_steps){
  # Make violin plots of QC stats
  p2 <- plotGroups(
      ArchRProj = proj, 
      groupBy = "Sample", 
      colorBy = "cellColData", 
      name = "TSSEnrichment",
      plotAs = "violin",
      alpha = 0.4,
      addBoxPlot = TRUE
  )
  p3 <- plotGroups(
      ArchRProj = proj, 
      groupBy = "Sample", 
      colorBy = "cellColData", 
      name = "nFrags",
      plotAs = "violin",
      alpha = 0.4,
      addBoxPlot = TRUE
  )
  plotPDF(p2,p3, name = "QC-Sample-Statistics-TSS-Sample-Violin.pdf", ArchRProj = proj, addDOC = FALSE, width = 8, height = 6)

  # For reproducibility, can use previously defined set of cells after doublet removal
  if (pre_defined_doublets) {
    proj <- proj[cells_to_include$V1,]
  } else {
    proj <- filterDoublets(proj, filterRatio = 1.2)
  }

  p2 <- plotGroups(
      ArchRProj = proj, 
      groupBy = "Sample", 
      colorBy = "cellColData", 
      name = "TSSEnrichment",
      plotAs = "violin",
      alpha = 0.4,
      addBoxPlot = TRUE
  )
  p3 <- plotGroups(
      ArchRProj = proj, 
      groupBy = "Sample", 
      colorBy = "cellColData", 
      name = "nFrags",
      plotAs = "violin",
      alpha = 0.4,
      addBoxPlot = TRUE
  )
  plotPDF(p2,p3, name = "QC-Sample-Statistics-TSS-Sample-Violin-Doublets-Filtered.pdf", ArchRProj = proj, addDOC = FALSE, width = 18, height = 14)

  pal <- paletteDiscrete(unique(getCellColData(proj)$Sample))
  for (i in names(pal)){
    pal[i] <- "#D51F26"
  }
  p2 <- plotGroups(
      ArchRProj = proj, 
      groupBy = "Sample", 
      name = "TSSEnrichment",
      colorBy = "cellColData", 
      plotAs = "violin",
      alpha = 1,
      addBoxPlot = FALSE, pal = pal)
  p2 <- p2+geom_boxplot(outlier.shape = NA, alpha = 1)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))

  plotPDF(p2, name = "Temp-QC-Sample-Statistics-TSS-Sample-Violin-Doublets-Filtered.pdf", ArchRProj = proj, addDOC = FALSE, width = 30, height = 20)
}



################################################################################################################
#..............................................................................................................#
################################################################################################################
# 5) Initial Dimensionality reductions and clustering
if (5 in execute_steps){
  proj <- addIterativeLSI(
      ArchRProj = proj,
      useMatrix = "TileMatrix", 
      name = "IterativeLSI", 
      iterations = 2, 
      clusterParams = list(
          resolution = c(0.2), 
          sampleCells = 10000, #nCells(proj), 
          n.start = 10
      ), 
      varFeatures = 25000, 
      dimsToUse = 1:30,
      sampleCellsPre = 50000,
      sampleCellsFinal = 50000, force = TRUE
  )

  proj <- addUMAP(
      ArchRProj = proj, 
      reducedDims = "IterativeLSI", 
      name = "UMAP", 
      nNeighbors = 30, 
      minDist = 0.5, 
      metric = "cosine", force = TRUE
  )

  # Identify Clusters from Iterative LSI
  proj <- addClusters(
      input = proj,
      reducedDims = "IterativeLSI",
      method = "Seurat",
      name = "Clusters",
      resolution = 1.7, sampleCells = 50000
  )

  cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Sample))
  library(pheatmap)
  cM <- cM / Matrix::rowSums(cM)
  p <- pheatmap::pheatmap(
      mat = as.matrix(cM), 
      color = paletteContinuous("whiteBlue"), 
      border_color = "black"
  )
  plotPDF(p, name = "Sample-Cluster-Confusion-Matrix", width = 6, height = 6,  ArchRProj = proj, addDOC = FALSE)


  # Plot the UMAP Embedding with Metadata Overlayed such as Experimental Sample and Clusters.
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
  p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
  p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DiseaseState", embedding = "UMAP")
  p8 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
  p9 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")
  p10 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Location", embedding = "UMAP")
  p11 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Donor", embedding = "UMAP")
  plotPDF(p1,p2,p3,p8,p9,p10,p11, name = "UMAP-Samples-Clusters", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)

  p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DiseaseState", pal = c("#89288F", "#D7CEC7", "#D51F26", "#D7CEC7"), embedding = "UMAP")
  plotPDF(p5, name = "UMAP-Disease_State_Only2", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)

  p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DiseaseState", pal = c("#D51F26", "#D7CEC7", "#89288F", "#D7CEC7"), embedding = "UMAP")
  plotPDF(p5, name = "UMAP-Disease_State_Only3", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)

  saveArchRProject(ArchRProj = proj, outputDirectory = "all_samples_HuBMAP_HTAN", load = FALSE, overwrite = FALSE)
}


################################################################################################################
#..............................................................................................................#
################################################################################################################
# 6) Marker Genes plotting and analysis
if (6 in execute_steps){

  markerGenesImmune  <- c(

      "PAX5", "MS4A1", "CD19", "IGLL5", "VPREB3", #B-Cell Trajectory
      "TPSAB1", "HDC", "CTSG", "CMA1", "KRT1", "IL1RAPL1", "GATA2", #Mast
      "SERPINA9", "HRK", "HTR3A", "TCL6", "CD180", "FCRLA", #GC
      "CMA1", "IL1RAPL1", "CD69", #CD69+ Mast
      "KRT1", #CD69- Mast
      "CD207", #DC2
      "KLRF1", "SH2D1B", "SH2D1B", #NKs
      "SSR4", "IGLL5", "IGLL1", "AMPD1",#Plasma
      "CD14", "CLEC9A", "FCGR1A", "LILRB2", "CD209", "CD1E", #Monocytes
      "S100A8", "S100A9", # Inflammatory Monocytes
      "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "TBX21", "IL7R", "CD4", "CD2", #TCells
      "BATF","TNFRSF4", "FOXP3","CTLA4","LAIR2", # Tregs
      "FOLR2","FABP3","PLA2G2D" #Macrophages
    )

  markerGenesStromal  <- c("CD44",
      "TGFB1", "BMP7", "MAP3K2",
      "COL6A1", "CD36", 
      "CD34",
      "FAP", "CBLN2", "SPOCK1", "ACSS3", # Fibroblast
      "SYT10", "SOSTDC1", "DES", "TAGLN", #Myofibroblasts
      "SELP", "ZNF385D", "FAM155A", "GALNT15", "MADCAM1", "CORT", #Post capillary venules
      "COX4I2", "KCNJ8", "HIGD1B", "RGS5", "NOTCH3", "HEYL", "FAM162B", #Pericytes
      "FAM110D", "INHBB", "NPR1", "NOVA2", "GPIHBP1", "SOX17", #endothelial
      "S100A1", # nerves
      "RSPO3", "CCL11", "WNT5B", "BMP4", "CHI3L1", "ACTA2", "WNT2B", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", 
      "MADCAM1", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A"
    )

  markerGenesEpithelial  <- c(
      "DCLK1", #Tuft
      "GNG13","SH2D7","SH2D6","TRPM5","AZGP1","KRT18","BMX","PSTPIP2","LRMP","PTGS1","IL17RB","HCK","PLCG2","ANXA13", # Tuft
      "KLK1","ITLN1","WFDC2","CLCA1","LRRC26","RETNLB","SPINK4","AGR2", # Immature Goblet
      "MUC2", "TFF1", "FCGBP","TBX10","FFAR4","MUC1","SMIM6","CAPN8","SYTL2","LGALS9B","BCAS1", # Goblet
      "FAP", # Fibroblast
      "CA1", # E.Immature_Enterocytes
      "RAB6B", #Enterocytes
      "CRYBA2","SCGN","FEV","CHGA","GCG","SCG5","PCSK1N","PYY","NEUROD1","MS4A8","DDC", #Enteroendocrine
      "CA2", "SI", # absorptive
      "SOX9", "CD34", #progenitor
      "MUC1", "KRT1", # General epithelial
      "LYZ", "DEFA5", # Paneth
      "GP2", "KRT7", # M cells
      "NTRK2","CCL23","CCL20","POLD1","GJB3","KCNE2","RNF207","TNFRSF11A","AKR1C2","SPINK5","NOXO1", # M cells
      "BEST4", "CA7","OTOP2","OTOP3", "MYOM1","MT1G","MT1H", # Best4+ enterocytes
      "DES", #smooth muscle
      "APC", # colon tumor surpressor
      "SMOC2", "RGMB", "LGR5", "ASCL2", #stem
      "FOXA1","EPCAM", "LGR5", "KRT8", "KRT18",
      "COL1A1", "COL1A2", "COL6A1", "COL6A2", "VWF", "PLVAP", "CDH5", "S100B"
    )

  proj <- addImputeWeights(proj)

  p <- plotEmbedding(
      ArchRProj = proj, 
      colorBy = "GeneScoreMatrix", 
      name = markerGenesImmune, 
      embedding = "UMAP",
      imputeWeights = getImputeWeights(proj)
  )
  plotPDF(plotList = p, 
      name = "Plot-UMAP-Marker-Genes-Immune-W-Imputation.pdf", 
      ArchRProj = proj, 
      addDOC = FALSE, width = 5, height = 5)

  p <- plotEmbedding(
      ArchRProj = proj, 
      colorBy = "GeneScoreMatrix", 
      name = markerGenesStromal, 
      embedding = "UMAP",
      imputeWeights = getImputeWeights(proj)
  )
  plotPDF(plotList = p, 
      name = "Plot-UMAP-Marker-Genes-Stromal-W-Imputation.pdf", 
      ArchRProj = proj, 
      addDOC = FALSE, width = 5, height = 5)

  p <- plotEmbedding(
      ArchRProj = proj, 
      colorBy = "GeneScoreMatrix", 
      name = markerGenesEpithelial, 
      embedding = "UMAP",
      imputeWeights = getImputeWeights(proj)
  )
  plotPDF(plotList = p, 
      name = "Plot-UMAP-Marker-Genes-Epithelial-W-Imputation.pdf", 
      ArchRProj = proj, 
      addDOC = FALSE, width = 5, height = 5)
}




