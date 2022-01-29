# Script to create a normal colon epithelial reference for HTAN project
# ml R/3.6.1
# ml python
# ml py-macs2/2.1.1_py27

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Set things up
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
`%notin%` <- Negate(`%in%`)

#Set/Create Working Directory to Folder
setwd("./projects/")

#Load Genome Annotations
addArchRGenome("hg38")

#Set Threads to be used
addArchRThreads()

# Define subscript to help keep track of dimensionality reducitons
subscript = "hubmap_colon_epithelial"

# Read in cells in hubmap normal reference
hubmap_epithelial_cells_doublets_removed <- read.table("./final_hubmap_epithelial_cells_in_HTAN_paper.txt", stringsAsFactors = FALSE)$V1

# Read in metadata
metadata <- read.table("./hubmap_htan_metadata_atac_and_rna_final.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)

# Path to regev RNA data, see https://github.com/cssmillie/ulcerative_colitis for download of seurat object
Regev_RNA_Path <- "/oak/stanford/groups/wjg/wbecker/other/scATAC/other_datasets/Regev_lab_Cell_2019/epithelialRNAse.rds"

# Path to hubmap normal colon reference
hubmap_RNA_Path <- "/oak/stanford/groups/wjg/wbecker/other/scRNA/analysis/epithelial_normal_final/set1_11epithelial_iterativeLSI/diet_clustered_normal_colon_proj_seurat.rds"

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 1) Create project
ArrowPath <- "./projects/ArrowFiles/"
ArrowFiles <- c(
    "B004-A-004-D.arrow",
    "B004-A-204-D.arrow",
    "B001-A-301-D.arrow",
    "B001-A-302-D.arrow",
    "B001-A-401-D.arrow",
    "B001-A-406-D.arrow",
    "B001-A-501-D.arrow",
    "B004-A-004-D-R2.arrow",
    "B004-A-008-D.arrow")

proj <- ArchRProject(
    ArrowFiles = paste0(ArrowPath, ArrowFiles), 
    outputDirectory = "all_samples_HuBMAP_normal_colon_epithelial_cells_doublets_removed"
)

# Subset project with previously defined non doublets
proj <- proj[hubmap_epithelial_cells_doublets_removed,]

saveArchRProject(ArchRProj = proj, outputDirectory = "all_samples_HuBMAP_normal_colon_epithelial_cells_doublets_removed", load = FALSE, overwrite = FALSE)

# proj <- loadArchRProject("all_samples_HuBMAP_normal_colon_epithelial_cells_doublets_removed")

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 2) Add metadata
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
saveArchRProject(ArchRProj = proj, outputDirectory = "all_samples_HuBMAP_normal_colon_epithelial_cells_doublets_removed", load = FALSE, overwrite = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 3) Define marker genes
markerGenes  <- c(
    "DCLK1", #Tuft
    "HTR3C", "HTR3E", "B4GALNT4", "OGDHL", "POU2F3", "TAS1R3", "TNNI1",
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
    "S100A1", # nerves
    "SMOC2", "RGMB", "LGR5", "ASCL2", #stem
    "FOXA1","EPCAM", "LGR5", "KRT8", "KRT18",
    "COL1A1", "COL1A2", "COL6A1", "COL6A2", "VWF", "PLVAP", "CDH5", "S100B"
  )

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 4) Dimensionality reduciton and clustering
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = paste("IterativeLSI", subscript, sep = ""), 
    iterations = 4, 
    clusterParams = list(
        resolution = c(0.1, 0.1, 0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 15000, sampleCellsPre = NULL,
    dimsToUse = 1:25, force = TRUE
)

proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = paste("IterativeLSI", subscript, sep = ""), 
    name = paste("UMAP", subscript, sep = ""), 
    nNeighbors = 30, 
    minDist = 0.4, 
    metric = "cosine", force=TRUE
)

proj <- addClusters(
    input = proj,
    reducedDims = paste("IterativeLSI", subscript, sep = ""),
    method = "Seurat",
    name = paste("Clusters_DR", subscript, sep = ""),
    resolution = 2, force=TRUE, nOutlier = 20, seed = 1
)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = paste("UMAP", subscript, sep = ""))
plotPDF(p1, name = paste(paste("qTESTPlot-UMAP-Sample-Clusters1", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj, outputDirectory = "all_samples_HuBMAP_normal_colon_epithelial_cells_current_doublets_removed", overwrite = FALSE, load = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 5) Plot Gene Scores
proj <- addImputeWeights(proj, reducedDims = paste0("IterativeLSI", subscript), k=9, sampleCells = floor(nCells(proj)/4))

p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = paste("UMAP", subscript, sep = ""),
    imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p, 
    name = paste(paste("Plot-UMAP-IterativeLSI-Marker-Genes-W-Imputation", subscript, sep = "-"), ".pdf", sep = ""), 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = paste("ClustersHarmonyhubmap_colon_epithelial"),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
plotPDF(heatmapGS, name = paste(paste("GeneScores-Marker-CellType-Heatmap", subscript, sep = "-"), ".pdf", sep = ""), width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 6) RNA integration

# Regev lab data
seRNA <- readRDS(Regev_RNA_Path)
proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = paste0("IterativeLSI", subscript),
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "Clusters",
    nameCell = "predictedCell_Regev", #Name of column where cell from scRNA is matched to each cell
    nameGroup = "predictedGroup_Regev", #Name of column where group from scRNA is matched to each cell
    nameScore = "predictedScore_Regev", #Name of column where prediction score from scRNA
    force = TRUE
)
pal <- paletteDiscrete(values = seRNA$Clusters)
p1 <- plotEmbedding(
    proj, 
    colorBy = "cellColData", 
    name = "predictedGroup_Regev", 
    pal = pal, embedding = paste0("UMAP", subscript)
)
plotPDF(p1, name = "Plot-UMAP-RNA-Integration-Regev-data.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Our Data
seRNA <- readRDS(hubmapRNA)
proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = paste0("IterativeLSI", subscript),
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "CellType",
    nameCell = "predictedCell_UnNewRNA", #Name of column where cell from scRNA is matched to each cell
    nameGroup = "predictedGroup_UnNewRNA", #Name of column where group from scRNA is matched to each cell
    nameScore = "predictedScore_UnNewRNA", #Name of column where prediction score from scRNA
    force = TRUE
)
pal <- paletteDiscrete(values = seRNA@meta.data$CellType)
p1 <- plotEmbedding(
    proj, 
    colorBy = "cellColData", 
    name = "predictedGroup_UnNewRNA", 
    pal = pal, embedding = paste0("UMAP", subscript)
)
plotPDF(p1, name = "Plot-UMAP-RNA-Integration-New-data.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 7) Harmony, clustering, and cell type labeling
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = paste("IterativeLSI", subscript, sep = ""),
    name = paste("Harmony", subscript, sep = ""),
    groupBy = "Sample", force = TRUE
)
proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = paste("Harmony", subscript, sep = ""), 
    name = paste("UMAPHarmony", subscript, sep = ""), 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", force = TRUE
)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = paste("UMAPHarmony", subscript, sep = ""))
p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "CellType", embedding = paste("UMAPHarmony", subscript, sep = ""))
p6 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup_UnNewRNA", embedding = paste("UMAPHarmony", subscript, sep = ""))
proj <- addClusters(
    input = proj,
    reducedDims = paste0("Harmony", subscript),
    method = "Seurat",
    name = paste0("ClustersHarmony", subscript),
    resolution = 2.5, nOutlier = 40, force = TRUE
)
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = paste("ClustersHarmony", subscript, sep = ""), embedding = paste("UMAPHarmony", subscript, sep = ""))
plotPDF(p3,p4,p5,p6, name = paste(paste("Plot-UMAP2Harmony-Sample-Clusters-doublet-removed", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


cellsNamesToAdd <- c()
clusterNamesToAdd <- c()

# Add annotations to main project
clusterNames <- sort(unique(proj$ClustersHarmonyhubmap_colon_epithelial))

clusterCellTypes <- c("Goblet",
    "Goblet", #10
    "Goblet",
    "TA2",
    "Stem",
    "TA2",
    "TA2",
    "TA1",
    "Enterocyte Progenitors",
    "Best4+ Enterocytes",
    "Best4+ Enterocytes",
    "Immature Enterocytes", 
    "Enterocytes", 
    "Immature Enterocytes", 
    "Immature Enterocytes", 
    "Enteroendocrine",
    "Immature Goblet",
    "Immature Goblet",
    "Secretory TA")
for (i in 1:length(clusterNames)){
    idxSample <- BiocGenerics::which(getCellColData(proj, paste("ClustersHarmony", subscript, sep = "")) %in% c(clusterNames[i]))
    cellsSample <- proj$cellNames[idxSample[[paste("ClustersHarmony", subscript, sep = "")]]]
    cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
    clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}
write.table(DataFrame(cellsNamesToAdd, clusterNamesToAdd), paste0("HuBMAPCellAnnotations_", subscript, '.tsv'), sep = '\t')
proj <- addCellColData(ArchRProj = proj, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "CellType", force = TRUE)

pal <- c("#D51F26")
names(pal) <- "Best4+ Enterocytes"
pal["Enterocyte Progenitors"] <- "#208A42"
pal["Enterocytes"] <- "#89288F"
pal["Enteroendocrine"] <- "#F47D2B"
pal["Goblet"] <- "#FEE500"
pal["Immature Enterocytes"] <- "#8A9FD1"
pal["Immature Goblet"] <- "#C06CAB"
pal["Secretory TA"] <- "#D24B27"
pal["Stem"] <- "#D8A767"
pal["TA1"] <- "#90D5E4"
pal["TA2"] <- "#89C75F"

p6 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "CellType", embedding = paste("UMAP", subscript, sep = ""), pal = pal)
p7 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "CellType", embedding = paste("UMAPHarmony", subscript, sep = ""), pal = pal)
plotPDF(p6,p7, name = paste(paste("Plot-UMAP-CellTypesNew", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = proj, outputDirectory = "all_samples_HuBMAP_normal_colon_epithelial_cells_current_doublets_removed", overwrite = FALSE, load = FALSE)


# plot confusion matrix heatmap to compare manual and integration labels
cM <- confusionMatrix(proj$predictedGroup_UnNewRNA, proj$CellType)
cM <- data.frame(t(cM) / Matrix::rowSums(t(cM)))
atac_order <- c("Stem","TA2","TA1","Enterocyte Progenitors","Immature Enterocytes","Enterocytes","Best4+ Enterocytes", "Secretory TA","Immature Goblet","Goblet", "Enteroendocrine")  
rna_order <- c("CyclingTA","Stem","TA2","TA1","Enterocyte.Progenitors","Immature.Enterocytes","Enterocytes","Best4..Enterocytes", "Immature.Goblet","Goblet", "Enteroendocrine", "Tuft")
cM <- cM[,rna_order]
cM <- cM[atac_order,]

paletteLength <- 256
myBreaks <- c(seq(0, 0.5, length.out=ceiling(paletteLength/2) + 1), 
              seq(0.51, 1, length.out=floor(paletteLength/2)))
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black", breaks = myBreaks, cluster_cols = FALSE, cluster_rows = FALSE
)
plotPDF(p, name = paste0("manual_vs_rna_cca_labels_percent_manual"), width = 7, height = 7,  ArchRProj = proj, addDOC = FALSE)


