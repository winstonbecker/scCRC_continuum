library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

project_directory <- "/oak/stanford/groups/wjg/wbecker/other/scATAC/reproducibility_check/projects/"
clustered_parent_project_path <- "./all_samples_HuBMAP_HTAN/"
subscript = "immune"
`%notin%` <- Negate(`%in%`)
# RNA data
seRNA3 <- readRDS("/oak/stanford/groups/wjg/wbecker/other/scRNA/analysis/immune_final/immune_normalize_and_scale/diet_clustered_full_colon_proj_seurat.rds")


#Set/Create Working Directory to Folder
setwd(project_directory)

#Load Genome Annotations
addArchRGenome("hg38")

#Set Threads to be used
addArchRThreads()

# Load previously defined archr project
proj <- loadArchRProject(path = clustered_parent_project_path)

##############################################################
# Define subset to explore
idxSample <- BiocGenerics::which(getCellColData(proj, "Clusters") %in% paste0("C", 31:35))
cellsSample <- proj$cellNames[idxSample[["Clusters"]]]

proj_immune <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsSample,
  outputDirectory = "all_samples_HuBMAP_HTAN_immune_cells"
)

saveArchRProject(ArchRProj = proj_immune, outputDirectory = "all_samples_HuBMAP_HTAN_immune_cells", load = FALSE, overwrite = FALSE)


##############################################################
# Define marker genes and motifs

markerGenes  <- c(
    "PAX5", "CD3D", "CTLA4","TOX"
  )

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# # LSI Projection and Clustering
proj_immune <- addIterativeLSI(
    ArchRProj = proj_immune,
    useMatrix = "TileMatrix", 
    name = paste("IterativeLSI", subscript, sep = ""), 
    iterations = 2, 
    clusterParams = list(
        resolution = c(0.2), 
        sampleCells = 20000, 
        n.start = 10
    ), 
    varFeatures = 25000,
    sampleCellsPre = NULL,
    dimsToUse = 1:30, force = TRUE
)

proj_immune <- addClusters(
    input = proj_immune,
    reducedDims = paste("IterativeLSI", subscript, sep = ""),
    method = "Seurat",
    name = paste("Clusters", subscript, sep = ""),
    resolution = 1.7, force=TRUE, nOutlier = 50, seed = 1
)

proj_immune <- addUMAP(
    ArchRProj = proj_immune, 
    reducedDims = paste("IterativeLSI", subscript, sep = ""), 
    name = paste("UMAP", subscript, sep = ""), 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", force=TRUE
)

p1 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = "Sample", embedding = paste("UMAP", subscript, sep = ""))
p2 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = paste("Clusters", subscript, sep = ""), embedding = paste("UMAP", subscript, sep = ""))
plotPDF(p1,p2, name = paste(paste("Plot-UMAP-Sample-Clusters", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_immune, addDOC = FALSE, width = 5, height = 5)


# Remove 3, 5, and 8--downstream analysis shows that these are likely doublets/not clearly defined clusters
idxSample <- BiocGenerics::which(getCellColData(proj_immune, "Clustersimmune") %notin% paste0("C", c(3,5,8)))
cellsSample <- proj_immune$cellNames[idxSample[["Clustersimmune"]]]
proj_immune <- proj_immune[cellsSample, ]


proj_immune <- addIterativeLSI(
    ArchRProj = proj_immune,
    useMatrix = "TileMatrix", 
    name = paste("IterativeLSI", subscript, sep = ""), 
    iterations = 2, 
    clusterParams = list(
        resolution = c(0.2), 
        sampleCells = 20000, 
        n.start = 10
    ), 
    varFeatures = 25000,
    sampleCellsPre = NULL,
    dimsToUse = 1:30, force = TRUE
)

proj_immune <- addClusters(
    input = proj_immune,
    reducedDims = paste("IterativeLSI", subscript, sep = ""),
    method = "Seurat",
    name = paste("Clusters", subscript, sep = ""),
    resolution = 1.7, force=TRUE, nOutlier = 50, seed = 1
)

proj_immune <- addUMAP(
    ArchRProj = proj_immune, 
    reducedDims = paste("IterativeLSI", subscript, sep = ""), 
    name = paste("UMAP", subscript, sep = ""), 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", force=TRUE
)


p1 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = "Sample", embedding = paste("UMAP", subscript, sep = ""))
p2 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = paste("Clusters", subscript, sep = ""), embedding = paste("UMAP", subscript, sep = ""))
plotPDF(p1,p2, name = paste(paste("Plot-UMAP-Sample-Clusters", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_immune, addDOC = FALSE, width = 5, height = 5)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# # Gene scores
proj_immune <- addImputeWeights(proj_immune, reducedDims = paste("IterativeLSI", subscript, sep = ""), k=9, sampleCells = floor(nCells(proj_immune)/4))

p <- plotEmbedding(
    ArchRProj = proj_immune, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = paste("UMAP", subscript, sep = ""),
    imputeWeights = getImputeWeights(proj_immune)
)

plotPDF(plotList = p, 
    name = paste(paste("Plot-UMAP-Marker-Genes-W-Imputation", subscript, sep = "-"), ".pdf", sep = ""), 
    ArchRProj = proj_immune, 
    addDOC = FALSE, width = 5, height = 5)

p <- plotEmbedding(
    ArchRProj = proj_immune, 
    colorBy = "GeneScoreMatrix", 
    name = markersImmunoglobulins, 
    embedding = paste("UMAP", subscript, sep = ""),
    imputeWeights = getImputeWeights(proj_immune)
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
plotPDF(heatmapGS, name = paste(paste("GeneScores-Marker-Heatmap-Clusters", subscript, sep = "-"), ".pdf", sep = ""), width = 8, height = 6, ArchRProj = proj_immune, addDOC = FALSE)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# RNA integration
proj_immune <- addGeneIntegrationMatrix(
    ArchRProj = proj_immune, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = paste0("IterativeLSI", subscript),
    seRNA = seRNA3,
    addToArrow = FALSE,
    groupRNA = "CellType",
    nameCell = "predictedCell_Un", #Name of column where cell from scRNA is matched to each cell
    nameGroup = "predictedGroup_Un", #Name of column where group from scRNA is matched to each cell
    nameScore = "predictedScore_Un", #Name of column where prediction score from scRNA
    force = TRUE
)

pal <- paletteDiscrete(values = seRNA3@meta.data$CellType)

pal <- c("#D51F26")
names(pal) <- "CD4+"
#pal["CD4+ Activated"] <- "#272E6A"
#pal["CD4+ Tfh PD1+"] <- "#208A42"
pal["CD8+"] <- "#89288F"
pal["DC"] <- "#F47D2B"
#pal["Exhausted T cells"] <- "#FEE500"
pal["GC"] <- "#8A9FD1"
#pal["Inflammatory Fibroblasts"] <- "#C06CAB"
pal["ILCs"] <- "#C06CAB"
pal["Macrophages"] <- "#D8A767"
pal["Mast"] <- "#89C75F"
#pal["Mast2"] <- "#89C75F"
pal["Memory B"] <- "#F37B7D"
pal["Niave B"] <- "#9983BD"
pal["NK"] <- "#D24B27"
pal["Tregs"] <- "#3BBCA8"
pal["Plasma"] <- "#6E4B9E"

p1 <- plotEmbedding(
    proj_immune, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    pal = pal, embedding = paste0("UMAP", subscript)
)
plotPDF(p1, name = "Plot-UMAP-RNA-IntegrationLSI_OurData.pdf", ArchRProj = proj_immune, addDOC = FALSE, width = 5, height = 5)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Add cell types to project
# initialize list
cellsNamesToAdd <- c()
clusterNamesToAdd <- c()

clusterNames <- paste0("C", c(1:7,16:26))
clusterCellTypes <- c("Macrophages","Macrophages","DC","Macrophages","Macrophages","Mast2","Mast",
    "Plasma", "Plasma", "Plasma", "Plasma", "Plasma", "Plasma", "Plasma", #16-22
    "GC", "Memory B", "Memory B", "Niave B") #23-26


for (i in 1:length(clusterNames)){
    idxSample <- BiocGenerics::which(getCellColData(proj_immune, paste("Clusters", subscript, sep = "")) %in% c(clusterNames[i]))
    cellsSample <- proj_immune$cellNames[idxSample[[paste("Clusters", subscript, sep = "")]]]
    cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
    clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
}

write.table(DataFrame(cellsNamesToAdd, clusterNamesToAdd), paste0("Bcell_Myeloid_CellAnnotations", subscript, '.tsv'), sep = '\t')

# Add in t cell annotations if done seperately
tcells <- read.table("/oak/stanford/groups/wjg/wbecker/other/scATAC/HuBMAP_HTAN_ENCODE_Only/archr_all_samples_final/cell_type_annotations/tCellAnnotationstcells.tsv", sep = '\t', stringsAsFactors = FALSE)
cellsNamesToAdd <- append(cellsNamesToAdd, tcells$cellsNamesToAdd)
clusterNamesToAdd <- append(clusterNamesToAdd, tcells$clusterNamesToAdd)
write.table(DataFrame(cellsNamesToAdd, clusterNamesToAdd), paste0("ImmuneCellAnnotations", subscript, '.tsv'), sep = '\t')


proj_immune <- addCellColData(ArchRProj = proj_immune, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "CellType", force = TRUE)

p6 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = "CellType", embedding = paste("UMAP", subscript, sep = ""))
plotPDF(p6, name = paste(paste("Plot-UMAP-CellTypes", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_immune, addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = proj_immune, load = FALSE, overwrite = FALSE)



