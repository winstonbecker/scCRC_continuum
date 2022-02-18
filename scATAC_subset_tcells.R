# Script to analysze scATAC t-cells from HTAN colon project with ArchR

# The following modules should be loaded:
# ml R/3.6.1
# ml python
# ml py-macs2/2.1.1_py27

# Load libraries
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

# Define locations and project names
project_directory <- "./projects/"
clustered_parent_project_path <- "./all_samples_HuBMAP_HTAN_immune_cells/"
subscript = "tcells"
`%notin%` <- Negate(`%in%`)

# Load RNA data to label t-cells
seRNA3 <- readRDS("./other_datasets/t_cell_seurat_object.rds") # this is the bcc t cell dataset from Yost et al.

#Set/Create Working Directory to Folder
setwd(project_directory)

#Load Genome Annotations
addArchRGenome("hg38")

#Set Threads to be used
addArchRThreads()

# Load previously defined archr project - should be the immune object, and will now subset only the T-cells
proj_immune <- loadArchRProject(path = clustered_parent_project_path)

# Define marker genes and motifs
markerGenes  <- c("PDCD1", "TOX", "GZMB", "CTLA4")

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Subset seurat object
idxSample <- BiocGenerics::which(getCellColData(proj_immune, "Clustersimmune") %in% paste0("C", c(8:15)))
cellsSample <- proj_immune$cellNames[idxSample[["Clustersimmune"]]]

proj_tcells <- subsetArchRProject(
  ArchRProj = proj_immune,
  cells = cellsSample,
  outputDirectory = "all_samples_HuBMAP_HTAN_t_cells"
)

proj_tcells <- loadArchRProject(path = "./all_samples_HuBMAP_HTAN_t_cells/")

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# LSI Projection, UMAP Reduction, and Clustering
proj_tcells <- addIterativeLSI(
    ArchRProj = proj_tcells,
    useMatrix = "TileMatrix", 
    name = paste("IterativeLSI", subscript, sep = ""), 
    iterations = 2, 
    clusterParams = list(
        resolution = c(0.2), 
        sampleCells = NULL, 
        n.start = 10
    ), 
    varFeatures = 25000,
    sampleCellsPre = NULL,
    dimsToUse = 1:30, force = TRUE
)
proj_tcells <- addClusters(
    input = proj_tcells,
    reducedDims = paste("IterativeLSI", subscript, sep = ""),
    method = "Seurat",
    name = paste("Clusters", subscript, sep = ""),
    resolution = 2.5, force=TRUE, nOutlier = 50, seed = 1
)
proj_tcells <- addUMAP(
    ArchRProj = proj_tcells, 
    reducedDims = paste("IterativeLSI", subscript, sep = ""), 
    name = paste("UMAP", subscript, sep = ""), 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", force=TRUE
)
p1 <- plotEmbedding(ArchRProj = proj_tcells, colorBy = "cellColData", name = "Sample", embedding = paste("UMAP", subscript, sep = ""))
p2 <- plotEmbedding(ArchRProj = proj_tcells, colorBy = "cellColData", name = paste("Clusters", subscript, sep = ""), embedding = paste("UMAP", subscript, sep = ""))
plotPDF(p1,p2, name = paste(paste("Plot-UMAP-Sample-Clusters", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_tcells, addDOC = FALSE, width = 5, height = 5)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Plot gene scores
proj_tcells <- addImputeWeights(proj_tcells, reducedDims = paste("IterativeLSI", subscript, sep = ""), k=9, sampleCells = floor(nCells(proj_tcells)))

p <- plotEmbedding(
    ArchRProj = proj_tcells, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = paste("UMAP", subscript, sep = ""),
    imputeWeights = getImputeWeights(proj_tcells)
)
plotPDF(plotList = p, 
    name = paste(paste("Plot-UMAP-Marker-Genes-W-Imputation", subscript, sep = "-"), ".pdf", sep = ""), 
    ArchRProj = proj_tcells, 
    addDOC = FALSE, width = 5, height = 5
)
markersGS <- getMarkerFeatures(
    ArchRProj = proj_tcells, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = paste0("Clusters", subscript),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
plotPDF(heatmapGS, name = paste(paste("GeneScores-Marker-Heatmap-Clusters", subscript, sep = "-"), ".pdf", sep = ""), width = 8, height = 6, ArchRProj = proj_tcells, addDOC = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# RNA integration
seRNA3 <- subset(seRNA3, subset = cluster %in% c("CD8_mem","Th17","Tregs","CD8_act","Naive","Tfh","CD8_eff","NK_cells", "CD8_ex","CD8_ex_act"))

proj_tcells <- addGeneIntegrationMatrix(
    ArchRProj = proj_tcells, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = paste0("IterativeLSI", subscript),
    seRNA = seRNA3,
    addToArrow = FALSE,
    groupRNA = "cluster",
    nameCell = "predictedCell_Yost", #Name of column where cell from scRNA is matched to each cell
    nameGroup = "predictedGroup_Yost", #Name of column where group from scRNA is matched to each cell
    nameScore = "predictedScore_Yost", #Name of column where prediction score from scRNA
    force = TRUE
)

pal <- paletteDiscrete(values = seRNA3@meta.data$cluster)
p1 <- plotEmbedding(
    proj_tcells, 
    colorBy = "cellColData", 
    name = "predictedGroup_Yost", 
    pal = pal, embedding = paste0("UMAP", subscript)
)
plotPDF(p1, name = "Plot-UMAP-RNA-IntegrationLSI_Yost.pdf", ArchRProj = proj_tcells, addDOC = FALSE, width = 5, height = 5)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Add cell types to project
# initialize list
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()

clusterNames <- sort(unique(proj_tcells$Clusterstcells))
clusterCellTypes <- c("ILCs" #1
    ,"CD4+" #10
    ,"CD4+" #11
    ,"CD4+" #12
    ,"CD4+ Tfh PD1+" #13
    ,"Tregs" #14
    ,"Tregs" #15
    ,"CD4+ Activated" #16
    ,"CD8+" #17
    ,"CD4+ Activated" #18
    ,"CD4+ Activated" #19
    ,"Exhausted T cells" #2
    ,"CD4+ Activated" #20
    ,"CD4+ Activated" #21
    ,"CD8+" #22
    ,"Exhausted T cells" #3
    ,"NK" #4
    ,"CD8+" #5
    ,"NK" #6
    ,"CD8+" #7
    ,"CD8+" #8
    ,"CD8+") #9

for (i in 1:length(clusterNames)){
    idxSample <- BiocGenerics::which(getCellColData(proj_tcells, paste("Clusters", subscript, sep = "")) %in% c(clusterNames[i]))
    cellsSample <- proj_tcells$cellNames[idxSample[[paste("Clusters", subscript, sep = "")]]]
    cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
    clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

write.table(DataFrame(cellsNamesToAdd, clusterNamesToAdd), paste0("tCellAnnotations", subscript, '.tsv'), sep = '\t')
proj_tcells <- addCellColData(ArchRProj = proj_tcells, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "CellType", force = TRUE)

p6 <- plotEmbedding(ArchRProj = proj_tcells, colorBy = "cellColData", name = "CellType", embedding = paste("UMAP", subscript, sep = ""))
plotPDF(p6,name = paste(paste("Plot-UMAP-CellTypes", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_tcells, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj_tcells, load = FALSE, overwrite = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
#Create Group Coverage Files that can be used for downstream analysis
proj_tcells <- addGroupCoverages(ArchRProj = proj_tcells, groupBy = "CellType", force = TRUE)

#Call Reproducible Peaks w/ Macs2
pathToMacs2 <- findMacs2()
proj_tcells <- addReproduciblePeakSet(
    ArchRProj = proj_tcells, groupBy = "CellType", force = TRUE, 
    pathToMacs2 = pathToMacs2
)

#Add Peak Matrix
proj_tcells <- addPeakMatrix(ArchRProj = proj_tcells, force = TRUE)

#Save updated project
saveArchRProject(ArchRProj = proj_tcells, load = FALSE, overwrite = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
#Add Deviation Matrix
proj_tcells <- addBgdPeaks(proj_tcells)
proj_tcells <- addDeviationsMatrix(
  ArchRProj = proj_tcells, 
  peakAnnotation = "Motif",
  force = TRUE
)

#Plot some motif deviations
motifs <- c("NR4A2")
variable_motifs <- c("SPI1_322", "TCF4_97", "ID3_38", "MESP2_94", "ESP1_69", "BCL11A_194")
t_exhausted_motifs <- c("NFKB2", "NFKB1", "MAF", "BATF", "IRF4", "NFATC1")

markerMotifs <- getFeatures(proj_tcells, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markervariable_motifs <- getFeatures(proj_tcells, select = paste(variable_motifs, collapse="|"), useMatrix = "MotifMatrix")
markert_exhausted_motifs <- getFeatures(proj_tcells, select = paste(t_exhausted_motifs, collapse="|"), useMatrix = "MotifMatrix")

p <- plotEmbedding(ArchRProj = proj_tcells, colorBy = "MotifMatrix", name = sort(markerMotifs), 
    embedding = paste("UMAP", subscript, sep = ""),imputeWeights = getImputeWeights(proj_tcells))
plotPDF(p, name = paste0("UMAP-Groups-Marker-MotifDeviations-w-Imputation-", subscript, ".pdf"), width = 5, height = 5,  ArchRProj = proj_tcells, addDOC = FALSE)
p <- plotEmbedding(ArchRProj = proj_tcells, colorBy = "MotifMatrix", name = sort(markervariable_motifs), 
    embedding = paste("UMAP", subscript, sep = ""),imputeWeights = getImputeWeights(proj_tcells))
plotPDF(p, name = paste0("UMAP-Groups-marker-variable_motifs-MotifDeviations-w-Imputation-", subscript, ".pdf"), width = 5, height = 5,  ArchRProj = proj_tcells, addDOC = FALSE)
p <- plotEmbedding(ArchRProj = proj_tcells, colorBy = "MotifMatrix", name = sort(markert_exhausted_motifs), 
    embedding = paste("UMAP", subscript, sep = ""),imputeWeights = getImputeWeights(proj_tcells))
plotPDF(p, name = paste0("UMAP-Groups-marker-t_exhausted_motifs-MotifDeviations-w-Imputation-", subscript, ".pdf"), width = 5, height = 5,  ArchRProj = proj_tcells, addDOC = FALSE)

#Save updated project
saveArchRProject(ArchRProj = proj_tcells, load = FALSE, overwrite = FALSE)



