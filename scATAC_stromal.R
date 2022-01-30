# script to analyze stromal scATAC cells

# The following should be loaded:
# ml R/3.6.1
# ml python
# ml py-macs2/2.1.1_py27

# Steps in script
# 0) Load and subset previously defined archr project
# 4) LSI Computation and Clustering
# 5) Gene score plots
# 6) RNA Integration - LSI, Harmony, Constrained LSI, Constrained Harmony (need to set clusters for constrained versions)
# 7) Add cell types and make cell fraction plots
# 8) Gene score heatmap with cell types
# 9) Browser tracks
# 10) Call Peaks
# 11) Add deviations matrix
# 12) Peak to gene linkages
# 13) Positive TF regulators
# 14) trajectory
# 14) Significance around RUNX1
# 15) Additional evidence for preCAFs
# 16) Add harmony and show it has very little effect on clustering

# Define the steps above you want to execute
execute_steps <- c(6,7,8,9,10,11)

library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
`%notin%` <- Negate(`%in%`)

# Set/Create Working Directory to Folder
setwd("/oak/stanford/groups/wjg/wbecker/other/scATAC/PostReview/projects/")

# Load Genome Annotations
addArchRGenome("hg38")

# Set Threads to be used
addArchRThreads()

# Things to set for subseting projects
subscript = "stromal"
new_project_save_name <- "all_samples_HuBMAP_HTAN_stromal_cells_current"
use_Regev_RNA <- FALSE

# Read set of cells
final_cell_list <- read.table("./final_set_stromal_cells.txt", stringsAsFactors = FALSE)$V1

# Load Seurat object containing RNA for RNA integration
seRNA_Regev <- readRDS("./fibRNAse.rds")
seRNA_htan <- readRDS("./stromal_cells_rna.rds")

# Define some example marker genes
markerGenes  <- c("CD44",
    "TGFB1", "BMP7", "MAP3K2",
    "COL6A1", 
    "CBLN2", "SPOCK1", "ACSS3", # Fibroblast
    "SYT10", "SOSTDC1", "DES", "TAGLN", #Myofibroblasts
    "CHI3L1","MMP3","PLAU","MMP1","TRAFD1","GBP1", # Inflammatory fibroblasts
    "SELP", "ZNF385D", "FAM155A", "GALNT15", "MADCAM1", "CORT", #Post capillary venules
    "MCAM", "COX4I2", "KCNJ8", "HIGD1B", "RGS5", "NOTCH3", "HEYL", "FAM162B", #Pericytes (MCAM-MCAM)
    "FAM110D", "INHBB", "NPR1", "NOVA2", "GPIHBP1", "SOX17", #endothelial 
    "PLVAP","CD36","DYSF","NRP1","SH3BP5","EXOC3L2","FABP5","VWA1","BAALC","PRSS23","RAPGEF4","APLN","HTRA1", #microvascular
    "S100A1", # nerves
    "CCL11", "WNT5B", "BMP4", "CHI3L1", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", 
    "MADCAM1", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A", "SOX6", "CCL8", "FABP5",
    "BMP2", "WNT5A", "BCL2", "OGN",
    "CD74", "CCL19",
    "CD34", "GREM1", "RSPO3", "WNT2B", # stem cell
    "SEMA3B", "SEMA3E", "SEMA3C", #epithelial growth
    "BMP5", "BMP4", "BMP2", #BMP Signaling
    "FBLN1", "PCOLCE2", "MFAP5", #extracellular matrix
    "COL4A5", "COL4A6", #basement membrane
    "CXCL12", "CXCL14", #immune attraction
    "CCL11", "CCL13", #eos recriuitment
    "CCL2", "CCL8", # myeloid recruitment
    "MYH11", "TAGLN", "ACTA2", "TPM4", #myofibroblasts
    "TWIST1", "WNT2", "FAP", #inflammatory/CAF
    "CXCL1", "CXCL2", "CYR61", "IL1B", "IL6",
    "HGF" #"VEGF", "GAF6"
)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 0) Load and subset previously defined archr project

if (0 %in% execute_steps){
    # Load previously defined archr project
    proj <- loadArchRProject(path = "./all_samples_HuBMAP_HTAN/")    
    
    # Define subset to explore
    idxSample <- BiocGenerics::which(getCellColData(proj, "Clusters") %in% c("C28", "C29", "C30"))
    cellsSample <- proj$cellNames[idxSample[["Clusters"]]]

    proj_stromal <- subsetArchRProject(
      ArchRProj = proj,
      cells = cellsSample,
      outputDirectory = new_project_save_name, dropCells = FALSE
    )

    # remove doublet/low-quality cells
    proj_stromal <- proj_stromal[final_cell_list,]

    saveArchRProject(ArchRProj = proj_stromal, load = FALSE)#, overwrite = FALSE)
} else {
    proj_stromal <- loadArchRProject(path = new_project_save_name)
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 4) Dim reduction and UMAP
if (4 %in% execute_steps){
    proj_stromal <- addIterativeLSI(
        ArchRProj = proj_stromal,
        useMatrix = "TileMatrix", 
        name = paste("IterativeLSI", subscript, sep = ""), 
        iterations = 3, 
        clusterParams = list(
            resolution = c(0.1, 0.2), 
            sampleCells = NULL, 
            n.start = 10
        ), 
        varFeatures = 20000, sampleCellsPre = NULL,
        dimsToUse = 1:30, force = TRUE
    )

    proj_stromal <- addClusters(
        input = proj_stromal,
        reducedDims = paste("IterativeLSI", subscript, sep = ""),
        method = "Seurat",
        name = paste("Clusters", subscript, sep = ""),
        resolution = 1.1, force=TRUE, nOutlier = 50, seed = 1
    )

    proj_stromal <- addUMAP(
        ArchRProj = proj_stromal, 
        reducedDims = paste("IterativeLSI", subscript, sep = ""), 
        name = paste("UMAP", subscript, sep = ""), 
        nNeighbors = 35, #40
        minDist = 0.5, #0.6
        metric = "cosine", force=TRUE
    )

    p1 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Sample", embedding = paste("UMAP", subscript, sep = ""))
    p2 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = paste("Clusters", subscript, sep = ""), embedding = paste("UMAP", subscript, sep = ""))
    p3 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Clusters", embedding = paste("UMAP", subscript, sep = ""))
    p4 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "DiseaseState", embedding = paste("UMAP", subscript, sep = ""))
    plotPDF(p1,p2,p3,p4, name = paste(paste("Plot-UMAP-Sample-Clusters-doublet-cluster-removed", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5) 
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 5) Plot Gene scores on embedding
if (5 %in% execute_steps){
    proj_stromal <- addImputeWeights(proj_stromal, reducedDims = paste("IterativeLSI", subscript, sep = ""), sampleCells = floor(nCells(proj_stromal)))

    p <- plotEmbedding(
        ArchRProj = proj_stromal, 
        colorBy = "GeneScoreMatrix", 
        name = markerGenes, 
        embedding = paste("UMAP", subscript, sep = ""),
        imputeWeights = getImputeWeights(proj_stromal)
    )

    plotPDF(plotList = p, 
        name = paste(paste("Plot-UMAP-Marker-Genes-W-Imputation", subscript, sep = "-"), ".pdf", sep = ""), 
        ArchRProj = proj_stromal, 
        addDOC = FALSE, width = 5, height = 5)
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 6) RNA Integration
if (6 %in% execute_steps){
  if (use_Regev_RNA){
    seRNA <- seRNA_Regev
    # Unconstrained integration
    proj_stromal <- addGeneIntegrationMatrix(
        ArchRProj = proj_stromal, 
        useMatrix = "GeneScoreMatrix",
        matrixName = "GeneIntegrationMatrix",
        reducedDims = paste0("IterativeLSI", subscript),
        seRNA = seRNA,
        addToArrow = FALSE,
        groupRNA = "Clusters",
        nameCell = "predictedCell_Un_Reg", #Name of column where cell from scRNA is matched to each cell
        nameGroup = "predictedGroup_Un_Reg", #Name of column where group from scRNA is matched to each cell
        nameScore = "predictedScore_Un_Reg", #Name of column where prediction score from scRNA
        force = TRUE
    )
    pal <- paletteDiscrete(values = colData(seRNA)$Clusters)
    p1 <- plotEmbedding(
        proj_stromal, 
        colorBy = "cellColData", 
        name = "predictedGroup_Un_Reg", 
        pal = pal, embedding = paste0("UMAP", subscript)
    )
    plotPDF(p1, name = "Plot-UMAP-RNA-Integration_Reg.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)

    clustGoblet <- c("C1", "C2", "C3", "C23")
    clustNonSecretory <- paste0("C", c(4:22,24:27))
    cGoblets = "Endothelial|Post\\-capillary\\ Venules|Microvascular|Pericytes"
    cNonSec = "WNT5B\\+\\ 2|WNT5B\\+\\ 1|WNT2B\\+\\ Fos\\-hi|WNT2B\\+\\ Fos\\-lo\\ 1|WNT2B\\+\\ Fos\\-lo\\ 2|RSPO3\\+|Glia|Inflammatory\\ Fibroblasts|Myofibroblasts"
    rnaGoblets <- colnames(seRNA)[grep(cGoblets, colData(seRNA)$Clusters)]
    rnaNonSecretory <- colnames(seRNA)[grep(cNonSec, colData(seRNA)$Clusters)]
    groupList <- SimpleList(
        Secretory = SimpleList(
            ATAC = proj_stromal$cellNames[proj_stromal$Clustersstromal %in% clustGoblet],
            RNA = rnaGoblets
        ),
        NonSecretory = SimpleList(
            ATAC = proj_stromal$cellNames[proj_stromal$Clustersstromal %in% clustNonSecretory],
            RNA = rnaNonSecretory
        )
    )
    proj_stromal <- addGeneIntegrationMatrix(
        ArchRProj = proj_stromal, 
        useMatrix = "GeneScoreMatrix",
        matrixName = "GeneIntegrationMatrix",
        reducedDims = paste0("IterativeLSI", subscript),
        seRNA = seRNA,
        addToArrow = FALSE, 
        groupList = groupList, #Constrain List
        groupRNA = "Clusters",
        nameCell = "predictedCell",
        nameGroup = "predictedGroup",
        nameScore = "predictedScore",
        force = TRUE
    )
    pal <- paletteDiscrete(values = colData(seRNA)$Clusters)
    p2 <- plotEmbedding(
        proj_stromal, 
        colorBy = "cellColData", 
        name = "predictedGroup", 
        pal = pal, embedding = paste0("UMAP", subscript)
    )
    plotPDF(p2, name = "Plot-UMAP-RNA-Integration-Constrained-Final.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)
  }
}


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 7) Add cell types and make cell fraction plots
if (7 %in% execute_steps){
    # Add cell annotations
    clusterCellTypes <- c("Lymphatic Endothelial Cells","Myofibroblasts 2",
      "Myofibroblasts 2","Inflammatory Fibroblasts",
      "Villus Fibroblasts WNT5B+","Pericytes",
      "Crypt Fibroblasts 3","Crypt Fibroblasts 3",
        "Crypt Fibroblasts 2","Cancer Associated Fibroblasts",
        "Crypt Fibroblasts RSPO3+","Endothelial",
        "Crypt Fibroblasts 1","Crypt Fibroblasts 1",
        "Crypt Fibroblasts 1","Endothelial",
        "Glia","Crypt Fibroblasts RSPO3+",
        "Unknown","Myofibroblasts 1",
        "Myofibroblasts 1","Myofibroblasts GREM1+")
    clusterNames <- sort(unique(proj_stromal$Clustersstromal))
    cellsNamesToAdd <- c()
    clusterNamesToAdd <- c()
    for (i in 1:length(clusterNames)){
        idxSample <- BiocGenerics::which(getCellColData(proj_stromal, paste("Clusters", subscript, sep = "")) %in% c(clusterNames[i]))
        cellsSample <- proj_stromal$cellNames[idxSample[[paste("Clusters", subscript, sep = "")]]]
        cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
        clusterNamesToAdd <- append(clusterNamesToAdd, rep(clusterCellTypes[i], length(cellsSample)))
    }

    write.table(DataFrame(cellsNamesToAdd, clusterNamesToAdd), paste0("cellAnnotationsRepeat", subscript, '.tsv'), sep = '\t')
    proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "CellType", force = TRUE)
    
    # Plot cell annotations
    pal <- c("#272E6A")
    names(pal) <- "Cancer Associated Fibroblasts"
    pal["Crypt Fibroblasts 1"] <- "#208A42"
    pal["Crypt Fibroblasts 2"] <- "#89288F"
    pal["Crypt Fibroblasts 3"] <- "#F47D2B"
    pal["Crypt Fibroblasts RSPO3+"] <- "#FEE500"
    pal["Endothelial"] <- "#8A9FD1"
    pal["Glia"] <- "#C06CAB"
    pal["Inflammatory Fibroblasts"] <- "#0C727C"
    pal["Lymphatic endothelial cells"] <- "#D8A767"
    pal["Myofibroblasts 1"] <- "#90D5E4"
    pal["Myofibroblasts 2"] <- "#89C75F"
    pal["Myofibroblasts GREM1+"] <- "#F37B7D"
    pal["Pericytes"] <- "#D24B27"
    pal["Unknown"] <- "#3BBCA8"
    pal["Villus Fibroblasts WNT5B+"] <- "#6E4B9E"
    p1 <- plotEmbedding(proj_stromal, colorBy = "cellColData", name = "CellType", embedding = paste0("UMAP", subscript), pal = pal)
    plotPDF(p1, name = "Plot-UMAP-CellType.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)


    seRNA <- seRNA_htan
    # Now do constrained integration to make sure the pericytes map to each other
    clustGoblet <- c("Pericytes")
    cGoblets = "Pericytes"
    rnaGoblets <- colnames(seRNA)[grep(cGoblets, seRNA@meta.data$CellType)]
    rnaNonSecretory <- colnames(seRNA)[!grepl(cGoblets, seRNA@meta.data$CellType)]
    groupList <- SimpleList(
        Secretory = SimpleList(
            ATAC = proj_stromal$cellNames[proj_stromal$CellType %in% clustGoblet],
            RNA = rnaGoblets
        ),
        NonSecretory = SimpleList(
            ATAC = proj_stromal$cellNames[proj_stromal$CellType %ni% clustGoblet],
            RNA = rnaNonSecretory
        )
    )
    proj_stromal <- addGeneIntegrationMatrix(
        ArchRProj = proj_stromal, 
        useMatrix = "GeneScoreMatrix",
        matrixName = "GeneIntegrationMatrix",
        reducedDims = paste0("IterativeLSI", subscript),
        seRNA = seRNA,
        addToArrow = TRUE,
        groupList = groupList, #Constrain List
        groupRNA = "CellType",
        nameCell = "predictedCell_C_OurData", #Name of column where cell from scRNA is matched to each cell
        nameGroup = "predictedGroup_C_OurData", #Name of column where group from scRNA is matched to each cell
        nameScore = "predictedScore_C_OurData",
        force = TRUE
    )

    #pal <- c("#D51F26")
    #names(pal) <- "Adipocytes"
    pal <- c("#272E6A")
    names(pal) <- "Cancer Associated Fibroblasts"
    #pal["Cancer Associated Fibroblasts"] <- "#272E6A"
    pal["Crypt Fibroblasts 1"] <- "#208A42"
    pal["Crypt Fibroblasts 2"] <- "#89288F"
    pal["Crypt Fibroblasts 3"] <- "#F47D2B"
    pal["Crypt Fibroblasts RSPO3+"] <- "#FEE500"
    pal["Endothelial"] <- "#8A9FD1"
    pal["Glia"] <- "#C06CAB"
    pal["Lymphatic endothelial cells"] <- "#D8A767"
    pal["Myofibroblasts 1"] <- "#90D5E4"
    pal["Myofibroblasts 2"] <- "#89C75F"
    pal["Myofibroblasts GREM1+"] <- "#F37B7D"
    pal["Neurons"] <- "#9983BD"
    pal["Pericytes"] <- "#D24B27"
    pal["Unknown"] <- "#3BBCA8"
    pal["Villus Fibroblasts WNT5B+"] <- "#6E4B9E"
    p2 <- plotEmbedding(
        proj_stromal, 
        colorBy = "cellColData", 
        name = "predictedGroup_C_OurData_Reproduce", 
        pal = pal, embedding = "UMAPstromal"
    )
    plotPDF(p2, name = "Reproduce-Plot-UMAP-RNA-Harmony-Constrained-Integration.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)

    #############################################################################################
    # Figure 1 plot cell type divided by disease state
    cM <- confusionMatrix(paste0(proj_stromal$CellType), paste0(proj_stromal$Sample))
    cM <- data.frame(cM / Matrix::rowSums(cM))
    samplenames <- c()
    cellnames <- c()
    values <- c()
    diseaseStates <- c()
    for (sample in colnames(cM)){
        diseaseState <- unique(getCellColData(proj_stromal)[getCellColData(proj_stromal)$Sample == gsub("\\.", "-", sample),]$DiseaseState)
        for (celltype in rownames(cM)){
            samplenames <- c(samplenames, sample)
            cellnames <- c(cellnames, celltype)
            diseaseStates <- c(diseaseStates, diseaseState)
            values <- c(values, cM[celltype,sample])
        }
    }
    data <- data.frame(samplenames,cellnames,values, diseaseStates)
    data <- data[order(sapply(data$diseaseStates, function(x) which(x == c("Normal", "Unaffected", "Polyp", "Adenocarcinoma")))), ]
    colOrder <-  c("Crypt Fibroblasts 1",
       "Crypt Fibroblasts 2",
       "Crypt Fibroblasts 3",
       "Crypt Fibroblasts RSPO3+",
       "Villus Fibroblasts WNT5B+",
       "Cancer Associated Fibroblasts",
       "Inflammatory Fibroblasts",
       "Myofibroblasts 1",
       "Myofibroblasts 2",
       "Myofibroblasts GREM1+",
       "Pericytes",
       "Lymphatic Endothelial Cells",
       "Endothelial",
       "Glia",
       "Unknown"
       )                

    data <- data[order(sapply(data$cellnames, function(x) which(x == colOrder))), ]
    data$cellnames <- factor(data$cellnames, levels = colOrder)

    counts <- table(getCellColData(proj_stromal)$CellType)[colOrder]
    axis_labels <- paste0(levels(data$cellnames), " (N = ", counts, ")")

    samplenames <- c()
    diseaseStates <- c()
    sampleCounts <- data.frame(table(paste0(proj_stromal$Sample)))
    sampleCounts$Freq <- sampleCounts$Freq/sum(sampleCounts$Freq)
    colnames(sampleCounts) <- c("samplenames", "values")
    sampleCounts$cellnames = "allcells"
    for (sample in sampleCounts$samplenames){
        diseaseState <- unique(data[data$samplenames == gsub("-", "\\.", sample),]$diseaseStates)
        samplenames <- c(samplenames, gsub("-", "\\.", sample))
        diseaseStates <- c(diseaseStates, as.character(diseaseState))
    }
    sampleCounts$samplenames = samplenames
    sampleCounts$diseaseStates = diseaseStates

    data <- rbind(data,sampleCounts)
    axis_labels <- c(axis_labels, "All Cells")

    data$samplenames <- factor(data$samplenames, levels = unique(data$samplenames))
    p <- ggplot(data, aes(fill=samplenames, y=values, x=cellnames)) + 
        geom_bar(position="stack", stat="identity") + theme_ArchR() + 
      ylab("Fraction Cells") + scale_x_discrete(labels= axis_labels) +
      xlab("Sample") + theme(axis.text.x = element_text(angle = 90)) + theme(axis.text=element_text(size=20,hjust=0.95,vjust=0.2),
            axis.title=element_text(size=20)) +
      scale_fill_manual("legend", values = c(
                            colorRampPalette(c("#007849","#1E392A"))(9), 
                            colorRampPalette(c("#008F95", "#0375B4","#062F4F"))(22),
                            colorRampPalette(c("#94618E", "#6E3667","#813772"))(52),
                            colorRampPalette(c("#E24E42", "#fe3401","#B82601"))(6)))
    plotPDF(p, name = paste(paste("Sample-CellType-Bar-Fraction_Samplev3", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_stromal, addDOC = FALSE, width = 15, height = 15)

    #############################################################################################
    # plot sample on the x axis divided by cell type
    cM <- confusionMatrix(paste0(getCellColData(proj_stromal)$CellType), paste0(getCellColData(proj_stromal)$Sample))
    cM <- data.frame(t(cM) / Matrix::rowSums(t(cM)))
    samplenames <- c()
    cellnames <- c()
    values <- c()
    diseaseStates <- c()
    for (sample in rownames(cM)){
        diseaseState <- unique(getCellColData(proj_stromal)[getCellColData(proj_stromal)$Sample == sample,]$DiseaseState)
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

    # Set appropriate order
    new_levels <-  c("Crypt.Fibroblasts.1",
       "Crypt.Fibroblasts.2",
       "Crypt.Fibroblasts.3",
       "Crypt.Fibroblasts.RSPO3.",
       "Villus.Fibroblasts.WNT5B.",
       "Cancer.Associated.Fibroblasts",
       "Inflammatory.Fibroblasts",
       "Myofibroblasts.1",
       "Myofibroblasts.2",
       "Myofibroblasts.GREM1.",
       "Pericytes",
       "Lymphatic.Endothelial.Cells",
       "Endothelial",
       "Glia",
       "Unknown"
       )

    data$cellnames <- factor(data$cellnames, levels = new_levels)

    test <- data.frame(table(getCellColData(proj_stromal)$Sample))
    rownames(test) <- test$Var1

    p <- ggplot(data, aes(fill=cellnames, y=values, x=samplenames)) + 
        geom_bar(position="stack", stat="identity") + theme_ArchR() + 
      ylab("Fraction Cells") +
      xlab("Sample") + theme(axis.text.x = element_text(angle = 90)) +
      scale_x_discrete(labels=paste0(levels(data$samplenames), " (n = ",  test[levels(data$samplenames),]$Freq, ")")) +
      scale_fill_manual("legend", values = c("#208A42", "#89288F" ,"#F47D2B" ,"#FEE500", "#6E4B9E" ,"#272E6A" ,"#0C727C", "#90D5E4", "#89C75F", "#F37B7D", "#D24B27", "#D8A767", "#8A9FD1", "#C06CAB","#3BBCA8" ))
    plotPDF(p, name = paste0("Sample-CellType-Bar-Fraction_Sample_Stromal", ".pdf"), ArchRProj = proj_stromal, addDOC = FALSE, width = 25, height = 20)

    #############################################################################################
    # Make celltype fraction plot with celltypes on x axis and stacked bars colored by sample/donor
    cM <- confusionMatrix(paste0(getCellColData(proj_stromal)$CellType), paste0(getCellColData(proj_stromal)$Sample))
    cM <- data.frame(cM / Matrix::rowSums(cM))
    samplenames <- c()
    cellnames <- c()
    values <- c()
    donors <- c()
    colnames(cM) <- unique(getCellColData(proj_stromal)$Sample)
    for (sample in unique(getCellColData(proj_stromal)$Sample)){
        Donor <- unique(getCellColData(proj_stromal)[getCellColData(proj_stromal)$Sample == sample,]$Donor)
        for (celltype in rownames(cM)){
            samplenames <- c(samplenames, sample)
            cellnames <- c(cellnames, celltype)
            donors <- c(donors, Donor)
            values <- c(values, cM[celltype,sample])
        }
    }
    data <- data.frame(samplenames,cellnames,values, donors)
    data <- data[order(sapply(data$donors, function(x) which(x == c("B001", "B004", "EP", "A001", "A002", "A008", "A010", "A014", "A015", "A018", "A022", "CRC1", "CRC2", "CRC3", "CRC4")))), ]

    colOrder <-  c("Crypt Fibroblasts 1",
       "Crypt Fibroblasts 2",
       "Crypt Fibroblasts 3",
       "Crypt Fibroblasts RSPO3+",
       "Villus Fibroblasts WNT5B+",
       "Cancer Associated Fibroblasts",
       "Inflammatory Fibroblasts",
       "Myofibroblasts 1",
       "Myofibroblasts 2",
       "Myofibroblasts GREM1+",
       "Pericytes",
       "Lymphatic Endothelial Cells",
       "Endothelial",
       "Glia",
       "Unknown"
       )        

    data$cellnames <- factor(data$cellnames, levels = colOrder)

    data <- data[order(sapply(data$cellnames, function(x) which(x == colOrder))), ]
    data$cellnames <- factor(data$cellnames, levels = colOrder)

    counts <- table(getCellColData(proj_stromal)$CellType)[colOrder]
    axis_labels <- paste0(levels(data$cellnames), " (N = ", counts, ")")

    data$samplenames <- factor(data$samplenames, levels = unique(data$samplenames))
    p <- ggplot(data, aes(fill=samplenames, y=values, x=cellnames)) + 
        geom_bar(position="stack", stat="identity") + theme_ArchR() + 
      ylab("Fraction Cells") + scale_x_discrete(labels= axis_labels)+
      xlab("Sample") + theme(axis.text.x = element_text(angle = 90)) + theme(axis.text=element_text(size=20,hjust=0.95,vjust=0.2),
            axis.title=element_text(size=20)) +
      scale_fill_manual("legend", values = c(
                            colorRampPalette(c("#007849","#1E392A"))(5), 
                            colorRampPalette(c("#008F95", "#0375B4","#062F4F"))(4),
                            colorRampPalette(c("#94618E", "#6E3667","#813772"))(4),
                            colorRampPalette(c("#E24E42", "#fe3401","#B82601"))(17),
                            colorRampPalette(c("#fbd2b6", "#F47D2B","#913f08"))(20),
                            colorRampPalette(c("#fff7b3", "#FEE500"))(2),
                            colorRampPalette(c("#dae1f1", "#8A9FD1"))(2),
                            colorRampPalette(c("#e0b8d6", "#C06CAB","#6b2e5c"))(13),
                            colorRampPalette(c("#b1e7df", "#3BBCA8","#1f6157"))(15),
                            colorRampPalette(c("#c1e8f0", "#90D5E4"))(2),
                            colorRampPalette(c("#3BBCA8"))(1),
                            colorRampPalette(c("#F37B7D"))(1),
                            colorRampPalette(c("#9983BD"))(1),
                            colorRampPalette(c("#D24B27"))(1),
                            colorRampPalette(c("#D8A767"))(1)))

    plotPDF(p, name = paste0("Sample-CellType-Bar-Fraction_Sample-ByDonor", ".pdf"), ArchRProj = proj_stromal, addDOC = FALSE, width = 8, height = 12)

    ##############################################################################################################################
    ##############################################################################################################################
    # Boxplots and Wilcoxin tests for differential abundance
    cM <- t(as.matrix(confusionMatrix(paste0(getCellColData(proj_stromal)$CellType), paste0(getCellColData(proj_stromal)$Sample))))
    cM["A001-C-023-D-R1",] <- cM["A001-C-023-D-R1",]+cM["A001-C-023-D-R2",]
    cM["A001-C-104-D-R1",] <- cM["A001-C-104-D-R1",]+cM["A001-C-104-D-R2",]
    cM["A001-C-124-D-R1",] <- cM["A001-C-124-D-R1",]+cM["A001-C-124-D-R2",]
    cM["A002-C-010-D",] <- cM["A002-C-010-D",]+cM["A002-C-010-D-R1",]+cM["A002-C-010-S2",]
    cM["A002-C-025-D",] <- cM["A002-C-025-D",]+cM["A002-C-025-S2",]
    cM["A002-C-116-D",] <- cM["A002-C-116-D",]+cM["A002-C-116-S2",]
    cM["A002-C-121-D",] <- cM["A002-C-121-D",]+cM["A002-C-121-S2",]
    cM["B004-A-004-D",] <- cM["B004-A-004-D",]+cM["B004-A-004-D-R2",]
    cM <- cM[rownames(cM) %ni% c("A001-C-023-D-R2","A001-C-104-D-R2","A001-C-124-D-R2","A002-C-010-S2","A002-C-010-D-R1","A002-C-025-S2","A002-C-116-S2","A002-C-121-S2", "B004-A-004-D-R2"),] #combine and remove replicate samples
    cM <- cM[rowSums(cM)>50,] #must have >50 cells
    cM <- cM/rowSums(cM)
    samples <- unique(getCellColData(proj_stromal,c("DiseaseState", "Sample")))

    rownames(samples) <- samples$Sample
    samples <- samples[samples$Sample%in% rownames(cM),]

    cellTypes <- unique(getCellColData(proj_stromal, "CellType")$CellType)
    sample1 <- "Normal"
    allPvalues <- data.frame(matrix(ncol = length(cellTypes), nrow = 0))
    colnames(allPvalues) <- cellTypes
    for (sample2 in c("Unaffected", "Polyp", "Adenocarcinoma")){
      pvals <- c()
      for (cellTypeCompared in cellTypes){
        pvals <- c(pvals, wilcox.test(cM[rownames(samples[samples$DiseaseState == sample1,, drop = FALSE]),cellTypeCompared], 
            cM[rownames(samples[samples$DiseaseState == sample2,, drop = FALSE]),cellTypeCompared], alternative = "two.sided")$p.value)
      }
      allPvalues[sample2,] <- pvals
    }
     
    data <- as.data.frame(cbind(cM, samples[rownames(cM),]))
    data$DiseaseState <- factor(data$DiseaseState, levels = c("Normal", "Unaffected", "Polyp", "Adenocarcinoma"))
    colnames(data) <- c("Crypt.Fibroblasts.1","Inflammatory.Fibroblasts","Myofibroblasts.2","Endothelial",
        "Crypt.Fibroblasts.RSPO3.","Pericytes", "Myofibroblasts.1","Villus.Fibroblasts.WNT5B.",
        "Myofibroblasts.GREM1.","Crypt.Fibroblasts.3", "Lymphatic.Endothelial.Cells","Cancer.Associated.Fibroblasts",
        "Glia","Crypt.Fibroblasts.2", "Unknown","DiseaseState","Sample")
    library(ggpubr)

    for (cellTypeCompared in c("Crypt.Fibroblasts.1","Inflammatory.Fibroblasts","Myofibroblasts.2","Endothelial",
        "Crypt.Fibroblasts.RSPO3.","Pericytes", "Myofibroblasts.1","Villus.Fibroblasts.WNT5B.",
        "Myofibroblasts.GREM1.","Crypt.Fibroblasts.3", "Lymphatic.Endothelial.Cells","Cancer.Associated.Fibroblasts",
        "Glia","Crypt.Fibroblasts.2", "Unknown")){
      pdf(paste0(paste("cell-type-boxplot", cellTypeCompared, sep = "-"), ".pdf"), width = 4, height = 4, onefile=F)
      p <- ggplot(data, aes_string(x = "DiseaseState", y = cellTypeCompared)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) + geom_jitter(color="black", size=2, alpha=0.9)+ theme_ArchR() +
        stat_compare_means(comparisons = list(c("Normal", "Unaffected"), c("Normal", "Polyp"), c("Normal", "Adenocarcinoma")), label.y = max(data[,cellTypeCompared])*c(1.1, 1.2, 1.3))
        #geom_signif(comparisons = list(c("Normal", "Unaffected"), c("Normal", "Polyp"), c("Normal", "Adenocarcinoma")), map_signif_level=TRUE)
      print(p)
      dev.off()
    }

    for (cellTypeCompared in c("Crypt.Fibroblasts.1","Inflammatory.Fibroblasts","Myofibroblasts.2","Endothelial",
        "Crypt.Fibroblasts.RSPO3.","Pericytes", "Myofibroblasts.1","Villus.Fibroblasts.WNT5B.",
        "Myofibroblasts.GREM1.","Crypt.Fibroblasts.3", "Lymphatic.Endothelial.Cells","Cancer.Associated.Fibroblasts",
        "Glia","Crypt.Fibroblasts.2", "Unknown")){
      dataTemp <- data[,c(cellTypeCompared, "DiseaseState")]
      colnames(dataTemp) <- c("cellTypeCompared", "DiseaseState")
      wilcox <- compare_means(cellTypeCompared ~ DiseaseState, ref.group = "Normal", p.adjust.method = "bonferroni", method='wilcox.test', data = dataTemp)
      pdf(paste0(paste("cell-type-boxplot-padj-", cellTypeCompared, sep = "-"), ".pdf"), width = 4, height = 4, onefile=F)
      p <- ggplot(data, aes_string(x = "DiseaseState", y = cellTypeCompared)) +
          geom_boxplot(position = position_dodge(), outlier.shape = NA) + geom_jitter(color="black", size=2, alpha=0.9) + stat_pvalue_manual(wilcox, label = "p.adj", y.position = max(data[,cellTypeCompared])*c(1.1, 1.2, 1.3))+ theme_ArchR()
      print(p)
      dev.off()
    }

    for (cellTypeCompared in c("Crypt.Fibroblasts.1","Inflammatory.Fibroblasts","Myofibroblasts.2","Endothelial",
        "Crypt.Fibroblasts.RSPO3.","Pericytes", "Myofibroblasts.1","Villus.Fibroblasts.WNT5B.",
        "Myofibroblasts.GREM1.","Crypt.Fibroblasts.3", "Lymphatic.Endothelial.Cells","Cancer.Associated.Fibroblasts",
        "Glia","Crypt.Fibroblasts.2", "Unknown")){
      dataTemp <- data[,c(cellTypeCompared, "DiseaseState")]
      dataTemp$DiseaseState[dataTemp$DiseaseState == "Normal"] <- "Unaffected"
      dataFull <- data
      dataFull$DiseaseState[data$DiseaseState == "Normal"] <- "Unaffected"
      colnames(dataTemp) <- c("cellTypeCompared", "DiseaseState")
      wilcox <- compare_means(cellTypeCompared ~ DiseaseState, ref.group = "Unaffected", p.adjust.method = "bonferroni", method='wilcox.test', data = dataTemp)
      pdf(paste0(paste("cell-type-boxplot-padj-normal-and-unaffected-", cellTypeCompared, sep = "-"), ".pdf"), width = 4, height = 4, onefile=F)
      p <- ggplot(dataFull, aes_string(x = "DiseaseState", y = cellTypeCompared)) +
          geom_boxplot(position = position_dodge(), outlier.shape = NA) + geom_jitter(color="black", size=2, alpha=0.9) + stat_pvalue_manual(wilcox, label = "p.adj", y.position = max(data[,cellTypeCompared])*c(1.1, 1.2))+ theme_ArchR()
      print(p)
      dev.off()
    }

    for (cellTypeCompared in c("Crypt.Fibroblasts.1","Inflammatory.Fibroblasts","Myofibroblasts.2","Endothelial",
        "Crypt.Fibroblasts.RSPO3.","Pericytes", "Myofibroblasts.1","Villus.Fibroblasts.WNT5B.",
        "Myofibroblasts.GREM1.","Crypt.Fibroblasts.3", "Lymphatic.Endothelial.Cells","Cancer.Associated.Fibroblasts",
        "Glia","Crypt.Fibroblasts.2", "Unknown")){
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

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 10) Call Peaks
if (10 %in% execute_steps){
    proj_stromal <- addGroupCoverages(ArchRProj = proj_stromal, groupBy = "CellType", force = TRUE)
    pathToMacs2 <- findMacs2()

    #Call Reproducible Peaks w/ Macs2
    proj_stromal <- addReproduciblePeakSet(
        ArchRProj = proj_stromal, groupBy = "CellType", force = TRUE, 
        pathToMacs2 = pathToMacs2
    )

    #Add Peak Matrix
    proj_stromal <- addPeakMatrix(ArchRProj = proj_stromal, force = TRUE)
    saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)

    #Get marker features
    markersPeaks <- getMarkerFeatures(
        ArchRProj = proj_stromal, 
        useMatrix = "PeakMatrix", 
        groupBy = "CellType",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon"
    )
    markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

    order <- c("Cancer Associated Fibroblasts", "Inflammatory Fibroblasts",  "Villus Fibroblasts WNT5B+", 
    "Crypt Fibroblasts 1","Crypt Fibroblasts 2","Crypt Fibroblasts 3","Crypt Fibroblasts RSPO3+",
    "Myofibroblasts 1","Myofibroblasts 2","Myofibroblasts GREM1+","Pericytes","Endothelial", 
    "Lymphatic Endothelial Cells", "Glia","Unknown")

    #Visualize Markers as a heatmap
    heatmapPeaks <- plotMarkerHeatmap(
      seMarker = markersPeaks[,order], 
      cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
      clusterCols = FALSE, transpose = TRUE
    )
    plotPDF(heatmapPeaks, name = paste("Peak-Marker-Heatmap", subscript, sep = "-"), width = 12, height = 5, ArchRProj = proj_stromal, addDOC = FALSE)

    #Motif Search in Peak Set and add to Peak Annotations
    if("Motif" %ni% names(proj_stromal@peakAnnotation)){
      proj_stromal <- addMotifAnnotations(ArchRProj = proj_stromal, motifSet = "cisbp", name = "Motif")
    }
    saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)

    #Identify Motif Enrichments
    enrichMotifs <- peakAnnoEnrichment(
        seMarker = markersPeaks,
        ArchRProj = proj_stromal,
        peakAnnotation = "Motif",
        cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
      )

    heatmapEM <- enrichHeatmap(enrichMotifs[,order], n = 7, transpose = TRUE, clusterCols = FALSE)
    plotPDF(heatmapEM, name = paste("Motifs-Enrich-Heatmap", subscript, sep = "-"), width = 12, height = 8, ArchRProj = proj_stromal, addDOC = FALSE)
}


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 11) Add deviations matrix
if (11 %in% execute_steps){
    proj_stromal <- addBgdPeaks(proj_stromal)
    proj_stromal <- addDeviationsMatrix(
      ArchRProj = proj_stromal, 
      peakAnnotation = "Motif",
      force = TRUE
    )
    saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)
    plotVarDev <- getVarDeviations(proj_stromal, plot = TRUE)
    plotPDF(plotVarDev, name = paste0("Variable-Motif-Deviation-Scores-", subscript), width = 5, height = 5, ArchRProj = proj_stromal, addDOC = FALSE)
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 12) Peak to gene linkages
if (12 %in% execute_steps){
    proj_stromal <- addPeak2GeneLinks(
        ArchRProj = proj_stromal,
        reducedDims = "IterativeLSIstromal"
    )

    pal <- c("#272E6A")
    names(pal) <- "Cancer Associated Fibroblasts"
    pal["Crypt Fibroblasts 1"] <- "#208A42"
    pal["Crypt Fibroblasts 2"] <- "#89288F"
    pal["Crypt Fibroblasts 3"] <- "#F47D2B"
    pal["Crypt Fibroblasts RSPO3+"] <- "#FEE500"
    pal["Endothelial"] <- "#8A9FD1"
    pal["Glia"] <- "#C06CAB"
    pal["Inflammatory Fibroblasts"] <- "#0C727C"
    pal["Lymphatic endothelial cells"] <- "#D8A767"
    pal["Myofibroblasts 1"] <- "#90D5E4"
    pal["Myofibroblasts 2"] <- "#89C75F"
    pal["Myofibroblasts GREM1+"] <- "#F37B7D"
    pal["Pericytes"] <- "#D24B27"
    pal["Unknown"] <- "#3BBCA8"
    pal["Villus Fibroblasts WNT5B+"] <- "#6E4B9E"

    p2g_matricies <- plotPeak2GeneHeatmap(ArchRProj = proj_stromal, groupBy = "CellType", k=15, palGroup = pal, returnMatrices = FALSE, nPlot = 100000)
    plotPDF(p2g_matricies, name = paste(paste("Peak2-Marker-Heatmap", subscript, sep = "-"), ".pdf", sep = ""), width = 8, height = 8, ArchRProj = proj_stromal, addDOC = FALSE)
    
    markerGenes2 <- c("WNT2", "VCAN", "SERPINE1", "RUNX1", "KIF26B", "PLXDC2", "THBS2")
    p <- plotBrowserTrack(
        ArchRProj = proj_stromal, 
        groupBy = "CellType", 
        geneSymbol = markerGenes2, 
        pal = pal,
        upstream = 50000,
        downstream = 50000,
        loops = getPeak2GeneLinks(proj_stromal)
    )
    plotPDF(plotList = p, 
        name = paste0("Plot-Tracks-Marker-Genes-with-Peak2Gene", subscript, ".pdf"), 
        ArchRProj = proj_stromal, 
        addDOC = FALSE, width = 5, height = 5)

    saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)


    p2g_matricies <- plotPeak2GeneHeatmap(ArchRProj = proj_stromal, groupBy = "CellType", k=15, palGroup = pal, returnMatrices = TRUE, nPlot = 100000)
    # Generate a dummy markerTest object
    markerTest <- getMarkerFeatures(
      ArchRProj = proj_stromal, 
      useMatrix = "PeakMatrix",
      groupBy = "CellType",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = "Pericytes",
      bgdGroups = "Glia"
    )

    for (i in 1:15){
      # create dummy summarized experiment where nothing is significant
      clusterMarkerTest <- markerTest
      assays(clusterMarkerTest)$Log2FC[,] <- 0
      assays(clusterMarkerTest)$FDR[,] <- 0

      #peaks from marker test
      new <- rowData(clusterMarkerTest)
      allpeaks <- paste0(new$seqnames, ":", new$start, "-", new$end)

      # get peaks in current cluster and then find which of the entire peak list are in the current cluster
      peaks <- p2g_matricies$Peak2GeneLinks[p2g_matricies$ATAC$kmeansId == i,]$peak
      clusters <- which(allpeaks %in% peaks)
      
      # Set Log2FC as significant for peaks in current cluster and compute enrichment of motifs in those peaks
      assays(clusterMarkerTest)$Log2FC[clusters,] <- 10000
      motifsUp <- peakAnnoEnrichment(
        seMarker = clusterMarkerTest,
        ArchRProj = proj_stromal,
        peakAnnotation = "Motif",
        cutOff = "Log2FC > 9999"
      )
      # plot enriched tfs
      df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
      df <- df[order(df$mlog10Padj, decreasing = TRUE),]
      df$rank <- seq_len(nrow(df))
      ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
        geom_point(size = 1) +
        ggrepel::geom_label_repel(data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),size = 1.5,nudge_x = 2,color = "black") + 
        theme_ArchR() + 
        ylab("-log10(P-adj) Motif Enrichment") + 
        xlab("Rank Sorted TFs Enriched") +
        scale_color_gradientn(colors = paletteContinuous(set = "comet"))
      plotPDF(ggUp, name = paste0("peak2gene_motif_enrichment_cluster", i), width = 5, height = 5, ArchRProj = proj_stromal, addDOC = FALSE)

      # save adjusted p values
      new_additionUp <- assays(motifsUp)$mlog10Padj
      colnames(new_additionUp) = paste0("C", i)
      if (i==1){
        save_motifs_up <- new_additionUp
      }
      if (i>1){
        save_motifs_up <- cbind(new_additionUp, save_motifs_up)
      }
    }

    df <- DataFrame(save_motifs_up)
    df <- df[apply(df[,-1], 1, function(x) !all(x==0)),]
    df <- df[apply(df[,-1], 1, function(x) sum(x>30)>0),]
    paletteLength <- 256
    myBreaks <- c(seq(0, 50, length.out=ceiling(paletteLength/2) + 1), 
                  seq(50.1, 100, length.out=floor(paletteLength/2)))
    p <- pheatmap::pheatmap(
        mat = as.matrix(df), 
        color = paletteContinuous("whiteBlue"), 
        border_color = "black", breaks = myBreaks, cluster_cols = TRUE
    )
    plotPDF(p, name = paste0("Motif-Enrichment_peak2gene-", "-kmeans-clusters"), width = 20, height = 20,  ArchRProj = proj_stromal, addDOC = FALSE)

}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 13) Positive TF regulators
if (13 %in% execute_steps){
    # compute TF regulators following archr tutorial https://www.archrproject.com/bookdown/identification-of-positive-tf-regulators.html
    seGroupMotif <- getGroupSE(ArchRProj = proj_stromal, useMatrix = "MotifMatrix", groupBy = "CellType")
    seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
    rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
        rowMaxs(assay(seZ) - assay(seZ)[,x])
    }) %>% Reduce("cbind", .) %>% rowMaxs

    corGIM_MM <- correlateMatrices(
        ArchRProj = proj_stromal,
        useMatrix1 = "GeneIntegrationMatrix",
        useMatrix2 = "MotifMatrix",
        reducedDims = "IterativeLSIstromal"
    )
    corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
    corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
    corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
    corGIM_MM$TFRegulator <- "NO"
    corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
    sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
    p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator, label = GeneIntegrationMatrix_name)) +
        geom_point() +
        geom_text(data = data.frame(subset(corGIM_MM, corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))), hjust=-0.2, vjust=0.2)+
        theme_ArchR() +
        geom_vline(xintercept = 0, lty = "dashed") + 
        scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
        xlab("Correlation To Gene Expression") +
        ylab("Max TF Motif Delta")+
        xlim(-1.3, 1.3) +
        scale_y_continuous(
          expand = c(0,0), 
          limits = c(0, max(corGIM_MM$maxDelta)*1.1)
    )

    plotPDF(p, name = paste(paste("positive_TF_regulators_gene_integration_matrix", subscript, sep = "-"), ".pdf", sep = ""), width = 8, height = 8, ArchRProj = proj_stromal, addDOC = FALSE)

    # motifs
    motifs <- c(
      "CEBPB","ERG","ETS1","ETS2","FLI1","FOXF1", "FOXF2","FOXN3-AS2","HIVEP1","HIVEP3","HOXD13","MAFF",   
      "MEF2C","MSC-AS1","NFIB","NFKB1","NFKB2","NR2F1-AS1", "RELB","RUNX1","RUNX2","SPIB","SPIC","TCF21")
    motifs <- c("CEBPB","RUNX1","RUNX2")
    motifPositions <- getPositions(proj_stromal)
    motifNames <- getFeatures(proj_stromal, "MotifMatrix")
    markerMotifs <- unlist(lapply(motifs, function(x) grep(x, motifNames, value = TRUE)))
    proj_stromal <- addImputeWeights(proj_stromal, reducedDims = paste("IterativeLSI", subscript, sep = ""), sampleCells = floor(nCells(proj_stromal)))

    p <- plotEmbedding(
        ArchRProj = proj_stromal, 
        colorBy = "MotifMatrix", 
        name = sort(markerMotifs), 
        embedding = paste("UMAP", subscript, sep = ""),
        imputeWeights = getImputeWeights(proj_stromal)
    )
    plotPDF(p, name = paste("Groups-Marker-MotifDeviations-w-Imputation-TF-regulators", subscript, sep = "-"), width = 5, height = 5,  ArchRProj = proj_stromal, addDOC = FALSE)

    p <- plotEmbedding(
        ArchRProj = proj_stromal, 
        colorBy = "GeneIntegrationMatrix", 
        name = motifs,
        pal = paletteContinuous("blueYellow"),
        embedding = paste("UMAP", subscript, sep = ""),
        imputeWeights = getImputeWeights(proj_stromal)
    )
    plotPDF(p, name = paste0("Gene-Integration-TF-regulators-", subscript, ".pdf"), width = 5, height = 5,  ArchRProj = proj_stromal, addDOC = FALSE)

    # Compute marker genes and marker motifs to provide info for 
    markerTest <- getMarkerFeatures(
        ArchRProj = proj_stromal, 
        useMatrix = "GeneIntegrationMatrix",
        groupBy = "CellType",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)")
    )

    markerMotifTest <- getMarkerFeatures(
        ArchRProj = proj_stromal, 
        useMatrix = "MotifMatrix",
        groupBy = "CellType",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)")
    )
    # Plot some violins
    pal <- c("#272E6A")
    names(pal) <- "Cancer Associated Fibroblasts"
    pal["Crypt Fibroblasts 1"] <- "#208A42"
    pal["Crypt Fibroblasts 2"] <- "#89288F"
    pal["Crypt Fibroblasts 3"] <- "#F47D2B"
    pal["Crypt Fibroblasts RSPO3+"] <- "#FEE500"
    pal["Endothelial"] <- "#8A9FD1"
    pal["Glia"] <- "#C06CAB"
    pal["Inflammatory Fibroblasts"] <- "#0C727C"
    pal["Lymphatic Endothelial Cells"] <- "#D8A767"
    pal["Myofibroblasts 1"] <- "#90D5E4"
    pal["Myofibroblasts 2"] <- "#89C75F"
    pal["Myofibroblasts GREM1+"] <- "#F37B7D"
    pal["Pericytes"] <- "#D24B27"
    pal["Unknown"] <- "#3BBCA8"
    pal["Villus Fibroblasts WNT5B+"] <- "#6E4B9E"

    motifs <- c("CEBPB")
    #motifs <- c("RUNX1")
    #motifs <- c("RUNX2")
    markerMotifs <- unlist(lapply(motifs, function(x) grep(x, motifNames, value = TRUE)))
    #motifOnly <- substr(markerMotifs[grepl("dev", markerMotifs)], 12,nchar(markerMotifs[grepl("dev", markerMotifs)]))
    #print(colnames(markerMotifTest)[assays(markerMotifTest[rowData(markerMotifTest)$name == motifOnly,])$FDR<0.05 & assays(markerMotifTest[rowData(markerMotifTest)$name == motifOnly,])$Log2FC>1])
    p <- plotGroups(ArchRProj = proj_stromal, groupBy = "CellType", colorBy = "MotifMatrix", pal = pal,
        name = sort(markerMotifs),plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE, imputeWeights = NULL)
    plotPDF(p, name = paste0("Violins-Groups-Marker-MotifDeviations-wo-Imputation-", motifs, ".pdf"), ArchRProj = proj_stromal, addDOC = FALSE, width = 6, height = 6)

    
    print(colnames(markerTest)[assays(markerTest[rowData(markerTest)$name == motifs,])$FDR<0.01 & assays(markerTest[rowData(markerTest)$name == motifs,])$Log2FC>1])
    p <- plotGroups(ArchRProj = proj_stromal, groupBy = "CellType", colorBy = "GeneIntegrationMatrix", name = motifs, pal = pal,
        plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE, imputeWeights = NULL)
    plotPDF(p, name = paste0("Violins-Groups-Gene-Integration-TF-regulators-wo-Imputation-", motifs, ".pdf"), ArchRProj = proj_stromal, addDOC = FALSE, width = 6, height = 6)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 14) trajectory
if (14 %in% execute_steps){
    trajectory <- c("Villus Fibroblasts WNT5B+", "Inflammatory Fibroblasts", "Cancer Associated Fibroblasts")
    proj_stromal <- addTrajectory(
        ArchRProj = proj_stromal, 
        name = "caf_trajectory", 
        groupBy = "CellType",
        trajectory = trajectory, 
        embedding = paste0("UMAP", subscript), 
        force = TRUE
    )

    p <- plotTrajectory(proj_stromal, trajectory = "caf_trajectory", embedding = paste0("UMAP", subscript), colorBy = "cellColData", name = "caf_trajectory")
    plotPDF(p, name = "Plot-caf_trajectory-Traj-UMAP.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)

    trajPM  <- getTrajectory(ArchRProj = proj_stromal, name = "caf_trajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)
    p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"), labelTop = 100)
    trajMM  <- getTrajectory(ArchRProj = proj_stromal, name = "caf_trajectory", useMatrix = "MotifMatrix", log2Norm = TRUE)
    p6 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), labelTop = 100)
    plotPDF(p4, p6, name = "Plot-caf_trajectory-Heatmaps.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 6, height = 10)
    
    trajGIM <- getTrajectory(ArchRProj = proj_stromal, name = "caf_trajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
    p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"), labelTop = 20, varCutOff = 0.6)
    plotPDF(p3, name = "Plot-caf_trajectory-GIM-Heatmap.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 6, height = 10)

    # get transitions
    test <- getCellColData(proj_stromal)[!is.na(getCellColData(proj_stromal)$caf_trajectory),c("CellType", "caf_trajectory")]
    as.matrix(test[order(test$caf_trajectory),"CellType"])
    new_df <- data.frame(cbind(test[order(test$caf_trajectory),"CellType"], c(1:3043), rep(1,3043)))
    new_df$X2 <- as.numeric(as.character(new_df$X2))
    p <- ggplot(new_df, aes(X2, X3)) + geom_tile(aes(fill = X1)) + 
        scale_fill_manual(values= c("#272E6A", "#0C727C", "#6E4B9E"))+theme_ArchR()+theme(axis.text.x = element_text(angle = 90))
    plotPDF(p, name = paste0("caf_trajectory_hm.pdf"), width = 5, height = 5, ArchRProj = proj_stromal, addDOC = FALSE)


    trajectory <- c("Cancer Associated Fibroblasts", "Villus Fibroblasts WNT5B+", "Inflammatory Fibroblasts")
    proj_stromal <- addTrajectory(
        ArchRProj = proj_stromal, 
        name = "control_caf_trajectory", 
        groupBy = "CellType",
        trajectory = trajectory, 
        reducedDims = paste0("IterativeLSI", subscript),
        #embedding = paste0("UMAP", subscript), 
        force = TRUE
    )

    p <- plotTrajectory(proj_stromal, trajectory = "control_caf_trajectory", embedding = paste0("UMAP", subscript), colorBy = "cellColData", name = "control_caf_trajectory")
    plotPDF(p, name = "Plot-control_caf_trajectory-Traj-UMAP.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)

    trajPM  <- getTrajectory(ArchRProj = proj_stromal, name = "control_caf_trajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)
    p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"), labelTop = 100)
    trajMM  <- getTrajectory(ArchRProj = proj_stromal, name = "control_caf_trajectory", useMatrix = "MotifMatrix", log2Norm = TRUE)
    p6 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), labelTop = 100)
    plotPDF(p4, p6, name = "Plot-control_caf_trajectory-Heatmaps.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 6, height = 10)
    
    trajGIM <- getTrajectory(ArchRProj = proj_stromal, name = "control_caf_trajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)
    p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"), labelTop = 20, varCutOff = 0.6)
    plotPDF(p3, name = "Plot-control_caf_trajectory-GIM-Heatmap.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 6, height = 10)

    # get transitions
    test <- getCellColData(proj_stromal)[!is.na(getCellColData(proj_stromal)$control_caf_trajectory),c("CellType", "control_caf_trajectory")]
    as.matrix(test[order(test$control_caf_trajectory),"CellType"])
    new_df <- data.frame(cbind(test[order(test$control_caf_trajectory),"CellType"], c(1:3043), rep(1,3043)))
    new_df$X2 <- as.numeric(as.character(new_df$X2))
    p <- ggplot(new_df, aes(X2, X3)) + geom_tile(aes(fill = X1)) + 
        scale_fill_manual(values= c("#272E6A", "#0C727C", "#6E4B9E"))+theme_ArchR()+theme(axis.text.x = element_text(angle = 90))
    plotPDF(p, name = paste0("caf_control_trajectory_hm.pdf"), width = 5, height = 5, ArchRProj = proj_stromal, addDOC = FALSE)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 14) Significance around RUNX1
if (14 %in% execute_steps){
    markersGS <- getMarkerFeatures(
        ArchRProj = proj_stromal, 
        useMatrix = "GeneScoreMatrix", 
        groupBy = "CellType",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
    )

    RUNX1GS <- markersGS[rowData(markersGS)$name == "RUNX1",]

    markersPeaks <- getMarkerFeatures(
        ArchRProj = proj_stromal, 
        useMatrix = "PeakMatrix", 
        groupBy = "CellType",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon"
    )

    markersPeaksRUNX1 <- markersPeaks[rowData(markersPeaks)$seqnames == "chr21",]
    markersPeaksRUNX1 <- markersPeaksRUNX1[rowData(markersPeaksRUNX1)$start >= (35049344-50000) & rowData(markersPeaksRUNX1)$start <= (35049344+50000) ,]
    colSums(assays(markersPeaksRUNX1)$FDR<0.1 & assays(markersPeaksRUNX1)$Log2FC>1)

    # Plot some violins
    pal <- c("#272E6A")
    names(pal) <- "Cancer Associated Fibroblasts"
    pal["Crypt Fibroblasts 1"] <- "#208A42"
    pal["Crypt Fibroblasts 2"] <- "#89288F"
    pal["Crypt Fibroblasts 3"] <- "#F47D2B"
    pal["Crypt Fibroblasts RSPO3+"] <- "#FEE500"
    pal["Endothelial"] <- "#8A9FD1"
    pal["Glia"] <- "#C06CAB"
    pal["Inflammatory Fibroblasts"] <- "#0C727C"
    pal["Lymphatic Endothelial Cells"] <- "#D8A767"
    pal["Myofibroblasts 1"] <- "#90D5E4"
    pal["Myofibroblasts 2"] <- "#89C75F"
    pal["Myofibroblasts GREM1+"] <- "#F37B7D"
    pal["Pericytes"] <- "#D24B27"
    pal["Unknown"] <- "#3BBCA8"
    pal["Villus Fibroblasts WNT5B+"] <- "#6E4B9E"
    p <- plotGroups(ArchRProj = proj_stromal, groupBy = "CellType", colorBy = "GeneScoreMatrix", name = c("WNT2", "RUNX1"), pal = pal,
        plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE, imputeWeights = NULL)
    plotPDF(p, name = paste0("plotGroups-GeneScoreMatrix-Wnt2-RUNX1", ".pdf"), ArchRProj = proj_stromal, addDOC = FALSE, width = 8, height = 10)


    p <- plotBrowserTrack(
        ArchRProj = proj_stromal, 
        groupBy = "CellType", 
        geneSymbol = markerGenes2, 
        features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)[c("Cancer Associated Fibroblasts", 
            "Inflammatory Fibroblasts", "Crypt Fibroblasts 1", "Crypt Fibroblasts 2", "Crypt Fibroblasts 3", "Crypt Fibroblasts RSPO3+",
            "Myofibroblasts 1", "Myofibroblasts 2", "Myofibroblasts GREM1+", "Villus Fibroblasts WNT5B+")],
        pal = pa
        upstream = 50000,
        downstream = 50000,
        loops = getPeak2GeneLinks(proj_stromal)
    )
    plotPDF(plotList = p, name = paste0("All-Fib-markerPeaks-On-Plot-Tracks-Marker-Genes-with-Peak2Gene", subscript, ".pdf"), ArchRProj = proj_stromal, addDOC = FALSE, width = 8, height = 10)


}


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 15) Additional evidence for preCAFs
if (15 %in% execute_steps){
    ############################################################################################################
    # First global measure will be correlation in average accessibility between different cell types
    peakSE <- getGroupSE(
      ArchRProj = proj_stromal,
      useMatrix = "PeakMatrix",
      groupBy = "CellType",
      divideN = TRUE,
      scaleTo = NULL,
      threads = getArchRThreads(),
      verbose = TRUE,
      logFile = createLogFile("getGroupSE")
    )
    #cor(assays(peakSE)$PeakMatrix)
    p <- pheatmap::pheatmap(
      mat = cor(assays(peakSE)$PeakMatrix), 
      color = paletteContinuous("solarExtra"), show_rownames = T, 
      border_color = NA, cluster_rows = FALSE, cluster_cols = FALSE #breaks = myBreaks
    )
    plotPDF(p, name = "cell_type_pearson_correlation.pdf", ArchRProj = proj_stromal, addDOC = FALSE, height=7,width=5)


    ############################################################################################################
    # Next we will identify a CAF Signature and examine how that varies across cell types
    markerTestPeak <- getMarkerFeatures(
        ArchRProj = proj_stromal, 
        useMatrix = "PeakMatrix",
        groupBy = "CellType",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        useGroups = "Cancer Associated Fibroblasts" #,bgdGroups = paste0(sampleNameB, "-", typeB)
    )
    markerTestSub <- markerTestPeak[(assays(markerTestPeak)$Log2FC>1) & (assays(markerTestPeak)$FDR <= 0.1),]

    # Geneate a peak matirix of just the peaks of interest
    PeakMatrix <- getMatrixFromProject(
        ArchRProj = proj_stromal,
        useMatrix = "PeakMatrix",
        useSeqnames = c("chr1"),
        verbose = TRUE,
        binarize = FALSE,
        threads = getArchRThreads(),
        logFile = createLogFile("getMatrixFromProject")
    )
    PeakMatrixKeep <- PeakMatrix[rowData(markerTestSub)[rowData(markerTestSub)$seqnames == "chr1",]$idx,]
    for (chrom in c(paste0("chr", 2:22), "chrX")){
        PeakMatrix <- getMatrixFromProject(
          ArchRProj = proj_stromal,
          useMatrix = "PeakMatrix",
          useSeqnames = c(chrom),
          verbose = TRUE,
          binarize = FALSE,
          threads = getArchRThreads(),
          logFile = createLogFile("getMatrixFromProject")
        )
        PeakMatrixKeep <- rbind(PeakMatrixKeep, PeakMatrix[rowData(markerTestSub)[rowData(markerTestSub)$seqnames == chrom,]$idx,])
    }
    # Add the scores to the project
    proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = colSums(assays(PeakMatrixKeep)$PeakMatrix), cells = paste0(names(colSums(assays(PeakMatrixKeep)$PeakMatrix))), name = "CAF_score", force = TRUE)
    proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = getCellColData(proj_stromal)$CAF_score/getCellColData(proj_stromal)$nFrags, cells = paste0(rownames(getCellColData(proj_stromal))), name = "CAF_score_depth_normalized", force = TRUE)

    pal <- c("#272E6A")
    names(pal) <- "Cancer Associated Fibroblasts"
    pal["Crypt Fibroblasts 1"] <- "#208A42"
    pal["Crypt Fibroblasts 2"] <- "#89288F"
    pal["Crypt Fibroblasts 3"] <- "#F47D2B"
    pal["Crypt Fibroblasts RSPO3+"] <- "#FEE500"
    pal["Endothelial"] <- "#8A9FD1"
    pal["Glia"] <- "#C06CAB"
    pal["Inflammatory Fibroblasts"] <- "#0C727C"
    pal["Lymphatic Endothelial Cells"] <- "#D8A767"
    pal["Myofibroblasts 1"] <- "#90D5E4"
    pal["Myofibroblasts 2"] <- "#89C75F"
    pal["Myofibroblasts GREM1+"] <- "#F37B7D"
    pal["Pericytes"] <- "#D24B27"
    pal["Unknown"] <- "#3BBCA8"
    pal["Villus Fibroblasts WNT5B+"] <- "#6E4B9E"

    order <- c("Cancer Associated Fibroblasts", "Inflammatory Fibroblasts",  "Villus Fibroblasts WNT5B+", 
    "Crypt Fibroblasts 1","Crypt Fibroblasts 2","Crypt Fibroblasts 3","Crypt Fibroblasts RSPO3+",
    "Myofibroblasts 1","Myofibroblasts 2","Myofibroblasts GREM1+","Pericytes","Endothelial", 
    "Lymphatic Endothelial Cells", "Glia","Unknown")

    # Violin plots of the CAF scores
    p1 <- plotGroups(
      ArchRProj = proj_stromal, 
      groupBy = "CellType", 
      colorBy = "cellColData", 
      name = "CAF_score_depth_normalized",
      plotAs = "violin", pal = pal
     )
    plotPDF(p1, name = "CAF-Scores-depth-normalized.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 10, height = 7)

    p1 <- plotGroups(
      ArchRProj = proj_stromal, 
      groupBy = "CellType",  
      colorBy = "cellColData", 
      name = "CAF_score",
      plotAs = "violin"
     )
    plotPDF(p1, name = "CAF-Scores-unnormalized.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 15, height = 10)

    proj_temp <- proj_stromal[getCellColData(proj_stromal)$DiseaseState == "Adenocarcinoma",]
    p1 <- plotGroups(
      ArchRProj = proj_temp, 
      groupBy = "CellType", 
      colorBy = "cellColData", 
      name = "CAF_score_depth_normalized",
      plotAs = "violin"
     )
    plotPDF(p1, name = "CAF-Scores-depth-normalized-excluding-CRC-samples.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 15, height = 10)

    ############################################################################################################
    # Additional marker evidence for preCAFs
    # MoreCAF markers
      markerTestGIM <- getMarkerFeatures(
        ArchRProj = proj_stromal, 
        useMatrix = "GeneIntegrationMatrix",
        groupBy = "CellType",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)"),
    )
    markerTestGIM_CAFs <- markerTestGIM[,"Cancer Associated Fibroblasts"]
    markerTestGIM_pCAFs <- markerTestGIM[,"Inflammatory Fibroblasts"]
    markerTestShared <- markerTestGIM[(assays(markerTestGIM_CAFs)$Log2FC>1.5) & (assays(markerTestGIM_CAFs)$FDR <= 0.1) & (assays(markerTestGIM_pCAFs)$Log2FC>1.5) & (assays(markerTestGIM_pCAFs)$FDR <= 0.1) & (rowSums(assays(markerTestGIM)$Log2FC>1)==2),]

    proj_stromal <- addImputeWeights(proj_stromal, reducedDims = paste("IterativeLSI", subscript, sep = ""), sampleCells = floor(nCells(proj_stromal)))

    p <- plotEmbedding(
        ArchRProj = proj_stromal, 
        colorBy = "GeneScoreMatrix", 
        name = c("MMP3", "INHBA", "KIF26B", "TRPV4", "LRRC15", "LOX"), 
        embedding = paste("UMAP", subscript, sep = ""),
        imputeWeights = getImputeWeights(proj_stromal)
    )

    plotPDF(plotList = p, 
        name = paste(paste("Plot-UMAP-CAF-Marker-Gene-Scores-W-Imputation", subscript, sep = "-"), ".pdf", sep = ""), 
        ArchRProj = proj_stromal, 
        addDOC = FALSE, width = 5, height = 5)

    pal <- c("#272E6A")
    names(pal) <- "Cancer Associated Fibroblasts"
    pal["Crypt Fibroblasts 1"] <- "#208A42"
    pal["Crypt Fibroblasts 2"] <- "#89288F"
    pal["Crypt Fibroblasts 3"] <- "#F47D2B"
    pal["Crypt Fibroblasts RSPO3+"] <- "#FEE500"
    pal["Endothelial"] <- "#8A9FD1"
    pal["Glia"] <- "#C06CAB"
    pal["Inflammatory Fibroblasts"] <- "#0C727C"
    pal["Lymphatic Endothelial Cells"] <- "#D8A767"
    pal["Myofibroblasts 1"] <- "#90D5E4"
    pal["Myofibroblasts 2"] <- "#89C75F"
    pal["Myofibroblasts GREM1+"] <- "#F37B7D"
    pal["Pericytes"] <- "#D24B27"
    pal["Unknown"] <- "#3BBCA8"
    pal["Villus Fibroblasts WNT5B+"] <- "#6E4B9E"

    motifs <- c("MMP3", "INHBA", "KIF26B", "TRPV4", "LRRC15", "LOX")
    p <- plotGroups(ArchRProj = proj_stromal, groupBy = "CellType", colorBy = "GeneIntegrationMatrix", name = motifs, pal = pal,
        plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE, imputeWeights = NULL)
    plotPDF(p, name = paste0("Violins-Groups-Gene-Integration-wo-Imputation-", "CAF-and-pCAF-Marker-Gene", ".pdf"), ArchRProj = proj_stromal, addDOC = FALSE, width = 6, height = 6)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 16) Add harmony and show it has very little effect on clustering
if (16 %in% execute_steps){
    proj_stromal <- addHarmony(
        ArchRProj = proj_stromal,
        reducedDims = paste("IterativeLSI", subscript, sep = ""),
        name = paste("Harmony", subscript, sep = ""),
        groupBy = "Sample", force = TRUE
    )
    proj_stromal <- addUMAP(
        ArchRProj = proj_stromal, 
        reducedDims = paste("Harmony", subscript, sep = ""), 
        name = paste("UMAPHarmony", subscript, sep = ""), 
        nNeighbors = 30, 
        minDist = 0.5, 
        metric = "cosine", force = TRUE
    )

        pal <- c("#272E6A")
        names(pal) <- "Cancer Associated Fibroblasts"
        pal["Crypt Fibroblasts 1"] <- "#208A42"
        pal["Crypt Fibroblasts 2"] <- "#89288F"
        pal["Crypt Fibroblasts 3"] <- "#F47D2B"
        pal["Crypt Fibroblasts RSPO3+"] <- "#FEE500"
        pal["Endothelial"] <- "#8A9FD1"
        pal["Glia"] <- "#C06CAB"
        pal["Inflammatory Fibroblasts"] <- "#0C727C"
        pal["Lymphatic Endothelial Cells"] <- "#D8A767"
        pal["Myofibroblasts 1"] <- "#90D5E4"
        pal["Myofibroblasts 2"] <- "#89C75F"
        pal["Myofibroblasts GREM1+"] <- "#F37B7D"
        pal["Pericytes"] <- "#D24B27"
        pal["Unknown"] <- "#3BBCA8"
        pal["Villus Fibroblasts WNT5B+"] <- "#6E4B9E"
        
    p5 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "CellType", pal = pal, embedding = paste("UMAPHarmony", subscript, sep = ""))
    plotPDF(p5, name = paste0("UMAP-Harmony-CellType", ".pdf"), ArchRProj = proj_stromal, addDOC = FALSE, width = 6, height = 6)
}


