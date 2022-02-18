# Script for downstream analysis of scATAC epithelial cells
# WRB 2020-2021

# The following modules should be loaded:
# ml R/3.6.1
# ml python
# ml py-macs2/2.1.1_py27

# The following steps are executed in this script
# 0) Load and subset previously defined archr project-
# 1) Add metadata-
# 2) Add annotations-
# 3) Call Peaks-
# 4) pairwise comparisons of polyps to normal-
# 5) pairwise comparisons of polyps to all unaffecteds
# 7) compute malignancy pseudotimes and plot results of differential analyses
# 8) Plot tracks

execute_steps <- c(0,1,2,3,4,5,6)

library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ComplexHeatmap)
library(ggpubr)
`%notin%` <- Negate(`%in%`)

#Set/Create Working Directory to Folder
setwd("./projects/")

#Load Genome Annotations
addArchRGenome("hg38")

#Set Threads to be used
addArchRThreads()

# Things to set for subseting projects
subscript <- "all_epithelial"
previous_project_name <- "./all_samples_HuBMAP_HTAN/"
previous_clusters <- c(1:27)
new_project_save_name <- "all_epithelial_HTAN_HuBMAP_final"
integrate_regev_data <- FALSE
integrate_our_RNA <- TRUE
metadata <- read.table("./hubmap_htan_metadata_atac_and_rna_final.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)

# Load annotation files: necessary if running step 2
prefix <- "./cell_type_annotations/"
polypCellAnotations <- read.table(paste0(prefix, "all_samples_cell_types_projections.csv"), sep = ',', row.names = 1, header = TRUE)
hubmapOnlyAnotations <- read.table(paste0(prefix, "HuBMAPCellAnnotations_hubmap_colon_epithelial.tsv"), sep = '\t', row.names = 1, header = TRUE)


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 0) Load and subset previously defined archr project
if (0 %in% execute_steps){
    proj <- loadArchRProject(path = previous_project_name)

    # Define subset to explore
    idxSample <- BiocGenerics::which(getCellColData(proj, "Clusters") %in% paste0("C", previous_clusters))
    cellsSample <- proj$cellNames[idxSample[["Clusters"]]]

    proj_stromal <- subsetArchRProject(
      ArchRProj = proj,
      cells = cellsSample,
      outputDirectory = new_project_save_name, dropCells = FALSE
    )

    saveArchRProject(ArchRProj = proj_stromal, outputDirectory = new_project_save_name, load = FALSE, overwrite = FALSE)
    
} else {
    proj_stromal <- loadArchRProject(path = new_project_save_name)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 1) Add metadata - most likely already done
if (1 %in% execute_steps){
  for (j in 2:dim(metadata)[2]){
    # initialize list
    cellsNamesToAdd <- c()
    annotationToAdd <- c()
    for (i in 1:dim(metadata)[1]){
      idxSample <- BiocGenerics::which(getCellColData(proj_stromal, "Sample") %in% metadata[i,"Sample"])
      cellsSample <- proj_stromal$cellNames[idxSample[["Sample"]]]
      cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
      annotationToAdd <- append(annotationToAdd, rep(metadata[i,j], length(cellsSample)))
    }

    proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = paste0(annotationToAdd), cells = paste0(cellsNamesToAdd), name = colnames(metadata)[j], force = TRUE)
  }
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 2) Add annotations
if (2 %in% execute_steps){
    # Get cell names and cluster names for projected samples
    cellsNamesToAdd <- polypCellAnotations$Cell
    clusterNamesToAdd <- polypCellAnotations$CellType

    # Since hubmap samples were used as the reference, we will remove those from the projected cell list
    non_normal_names <- rownames(getCellColData(proj_stromal)[getCellColData(proj_stromal)$DiseaseState != "Normal",])
    clusterNamesToAdd <- paste0(clusterNamesToAdd)[paste0(cellsNamesToAdd) %in% non_normal_names]
    cellsNamesToAdd <- paste0(cellsNamesToAdd)[paste0(cellsNamesToAdd) %in% non_normal_names]

    # Now add hubmap cell annotations
    cellsNamesToAdd <- c(cellsNamesToAdd, paste0(hubmapOnlyAnotations$cellsNamesToAdd))
    clusterNamesToAdd <- c(clusterNamesToAdd, paste0(hubmapOnlyAnotations$clusterNamesToAdd))

    # Make sure that all of these cells are actually in the project, which should be true
    clusterNamesToAdd <- paste0(clusterNamesToAdd)[paste0(cellsNamesToAdd) %in% rownames(getCellColData(proj_stromal))]
    cellsNamesToAdd <- paste0(cellsNamesToAdd)[paste0(cellsNamesToAdd) %in% rownames(getCellColData(proj_stromal))]

    # Add cell types to the cell col data
    proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = paste0(clusterNamesToAdd), cells = paste0(cellsNamesToAdd), name = "CellType", force = TRUE)

    # Add cell type disease state to metadata as well
    proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = paste(getCellColData(proj_stromal)$CellType, getCellColData(proj_stromal)$DiseaseState, "."), cells = paste0(rownames(getCellColData(proj_stromal, "CellType"))), name = "CellTypeDiseaseState", force = TRUE)

    # Remove any cells from the project that do not have an associated cell type (should be hubmap doublet cluster)
    proj_stromal <- proj_stromal[!is.na(getCellColData(proj_stromal, "CellType")$CellType),]
}
    
##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 3) Make cell type fraction plots
# Make celltype fraction plot with celltypes on x axis and stacked bars colored by sample/disease state (Figure 1)
# Make celltype fraction plot with celltypes on x axis and stacked bars colored by sample/donor
# metadata heatmap
# Boxplots and Wilcoxin tests for differential abundance
if (3 %in% execute_steps){
    ##############################################################################################################################
    ##############################################################################################################################
    # Make celltype fraction plot with celltypes on x axis and stacked bars colored by sample/disease state (Figure 1)
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
    colOrder <-  c("Stem","TA2","TA1","Enterocyte Progenitors","Immature Enterocytes","Enterocytes","Secretory TA",
       "Immature Goblet","Goblet","Enteroendocrine","Best4+ Enterocytes")
    data <- data[order(sapply(data$cellnames, function(x) which(x == colOrder))), ]
    data$cellnames <- factor(data$cellnames, levels = colOrder)

    counts <- table(getCellColData(proj_stromal)$CellType)[colOrder]
    axis_labels <- paste0(colOrder, " (N = ", counts, ")")

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
    sampleCounts <- sampleCounts[order(sapply(sampleCounts$diseaseStates, function(x) which(x == c("Normal", "Unaffected", "Polyp", "Adenocarcinoma")))), ]

    data <- rbind(data,sampleCounts)
    axis_labels <- c(axis_labels, "All Cells")

    data$samplenames <- factor(data$samplenames, levels = unique(data$samplenames))
    p <- ggplot(data, aes(fill=samplenames, y=values, x=cellnames)) + 
        geom_bar(position="stack", stat="identity") + theme_ArchR() + 
      ylab("Fraction Cells") + scale_x_discrete(labels= axis_labels)+
      xlab("Sample") + theme(axis.text.x = element_text(angle = 90)) + theme(axis.text=element_text(size=20,hjust=0.95,vjust=0.2),
            axis.title=element_text(size=20)) +
      scale_fill_manual("legend", values = c(
                            colorRampPalette(c("#007849","#1E392A"))(9), 
                            colorRampPalette(c("#008F95", "#0375B4","#062F4F"))(22),
                            colorRampPalette(c("#94618E", "#6E3667","#813772"))(52),
                            colorRampPalette(c("#E24E42", "#fe3401","#B82601"))(6)))
    plotPDF(p, name = paste(paste("Sample-CellType-Bar-Fraction_Samplev2", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_stromal, addDOC = FALSE, width = 20, height = 20)

    ##############################################################################################################################
    ##############################################################################################################################
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
    data <- data[order(sapply(data$donors, function(x) which(x == c("B001", "B004", "F", "A001", "A002", "A008", "A010", "A014", "A015", "A018", "A022", "CRC1", "CRC2", "CRC3", "CRC4")))), ]

    colOrder <-  c("Stem",
       "TA2",
       "TA1",
       "Enterocyte Progenitors",
       "Immature Enterocytes",
       "Enterocytes",
       "Secretory TA",
       "Immature Goblet",
       "Goblet",
       "Enteroendocrine",
       "Best4+ Enterocytes")   

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
    #metadata heatmap
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

    metadata <- read.table("./hubmap_htan_metadata_atac_and_rna_final.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)
    metadata <- metadata[metadata$Sample %in% data$samplenames,]
    rownames(metadata) <- metadata$Sample
    metadata <- metadata[levels(data$samplenames),c("Sample", "DiseaseState", "FAP", "Size", "Location", "Dysplasia", "HGD", "Pathologist0_HGD")]

    metadata[metadata == "y"] <- "Y"
    metadata[metadata == "n"] <- "N"
    metadata[metadata == ""] <- "na"

    library(ggplot2); library(reshape2)
    metadata3 <- melt(data.frame(metadata), id.var = "Sample")
    metadata3$Sample <- factor(metadata3$Sample, levels = unique(data$samplenames))

    p <- ggplot(metadata3, aes(Sample, variable)) + geom_tile(aes(fill = value),
       colour = "white") + scale_fill_manual(values= c("Normal"="#1E392A", "Unaffected" = "#0375B4", "Polyp" = "#6E3667", 
        "Adenocarcinoma" = "#fe3401","N"="#FEE500", "Y"="#8A9FD1", "na"="#D5D5D5", "Large"="#fbd2b6","Medium"="#F47D2B","Small"="#913f08",
        "Ascending"="#E0B8D6", "Hepatic flexure"="#D399C4", "Transverse"="#C67BB3", "Descending"="#AE5F9B", "Sigmoid"="#8D467B", "Rectum"="#6B2E5C"))+ theme_ArchR()+ theme(axis.text.x = element_text(angle = 90))
    plotPDF(p, name = paste0("metadata_heatmap.pdf"), width = 15, height = 15, ArchRProj = proj_stromal, addDOC = FALSE)

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
    colnames(data) <- c("Stem","SecretoryTA","Goblet","TA2","ImmatureGoblet",
      "TA1","EnterocyteProgenitors", "ImmatureEnterocytes","Enteroendocrine", 
      "Enterocytes","Best4posEnterocytes","DiseaseState","Sample")

    for (cellTypeCompared in c("Stem","SecretoryTA","Goblet","TA2","ImmatureGoblet","TA1","EnterocyteProgenitors", "ImmatureEnterocytes","Enteroendocrine", "Enterocytes","Best4posEnterocytes")){
      pdf(paste0(paste("cell-type-boxplot", cellTypeCompared, sep = "-"), ".pdf"), width = 4, height = 4, onefile=F)
      p <- ggplot(data, aes_string(x = "DiseaseState", y = cellTypeCompared)) +
        geom_boxplot(position = position_dodge(), outlier.shape = NA) + geom_jitter(color="black", size=2, alpha=0.9)+ theme_ArchR() +
        stat_compare_means(comparisons = list(c("Normal", "Unaffected"), c("Normal", "Polyp"), c("Normal", "Adenocarcinoma")), label.y = max(data[,cellTypeCompared])*c(1.1, 1.2, 1.3))
        #geom_signif(comparisons = list(c("Normal", "Unaffected"), c("Normal", "Polyp"), c("Normal", "Adenocarcinoma")), map_signif_level=TRUE)
      print(p)
      dev.off()
    }

    for (cellTypeCompared in c("Stem","SecretoryTA","Goblet","TA2","ImmatureGoblet",
      "TA1","EnterocyteProgenitors", "ImmatureEnterocytes","Enteroendocrine", 
      "Enterocytes","Best4posEnterocytes")){
      dataTemp <- data[,c(cellTypeCompared, "DiseaseState")]
      colnames(dataTemp) <- c("cellTypeCompared", "DiseaseState")
      wilcox <- compare_means(cellTypeCompared ~ DiseaseState, ref.group = "Normal", p.adjust.method = "bonferroni", method='wilcox.test', data = dataTemp)
      pdf(paste0(paste("cell-type-boxplot-padj-", cellTypeCompared, sep = "-"), ".pdf"), width = 4, height = 4, onefile=F)
      p <- ggplot(data, aes_string(x = "DiseaseState", y = cellTypeCompared)) +
          geom_boxplot(position = position_dodge(), outlier.shape = NA) + geom_jitter(color="black", size=2, alpha=0.9) + stat_pvalue_manual(wilcox, label = "p.adj", y.position = max(data[,cellTypeCompared])*c(1.1, 1.2, 1.3))+ theme_ArchR()
      print(p)
      dev.off()
    }

    for (cellTypeCompared in c("Stem","SecretoryTA","Goblet","TA2","ImmatureGoblet",
      "TA1","EnterocyteProgenitors", "ImmatureEnterocytes","Enteroendocrine", 
      "Enterocytes","Best4posEnterocytes")){
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

    for (cellTypeCompared in c("Stem","SecretoryTA","Goblet","TA2","ImmatureGoblet",
      "TA1","EnterocyteProgenitors", "ImmatureEnterocytes","Enteroendocrine", 
      "Enterocytes","Best4posEnterocytes")){
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
# 3) Call Peaks
if (3 %in% execute_steps){
    # Make groupings for peak calling
    # get current sample names, cell types, and disease states
    samples <- getCellColData(proj_stromal, "SimplifiedSampleName") # maybe add a peak call sample name to this...
    sample_names_only <- unique(samples$SimplifiedSampleName)
    disease_state <- getCellColData(proj_stromal, "DiseaseState")
    cell_types <- getCellColData(proj_stromal, "CellType")
    
    # Set the things that we want to call peaks in a sample specific way if we have enough cells
    sample_specific_epithelial_cell_types <- c("Stem",
         "TA2",
         "TA1",
         "Enterocyte Progenitors",
         "Immature Enterocytes",
         "Enterocytes",
         "Secretory TA",
         "Immature Goblet",
         "Goblet")

    # Set the things that we will call in a disease state specific way because they are rare
    disease_specific_epithelial_cell_types <- c(
         "Enteroendocrine",
         "Best4+ Enterocytes",
         "Tuft")

    # For epithelial cells we will add the sample name or disease state as specificed above
    cell_types[cell_types$CellType %in% sample_specific_epithelial_cell_types,] <- paste(samples[cell_types$CellType %in% sample_specific_epithelial_cell_types,], cell_types[cell_types$CellType %in% sample_specific_epithelial_cell_types,], sep = "-")
    cell_types[cell_types$CellType %in% disease_specific_epithelial_cell_types,] <- paste(disease_state[cell_types$CellType %in% disease_specific_epithelial_cell_types,], cell_types[cell_types$CellType %in% disease_specific_epithelial_cell_types,], sep = "-")

    # Some of these groupings dont have enough cells, so we will get most uncommon new groupings (i.e. dont call peaks on less than 200 cells)
    grouping_counts <- data.frame(table(cell_types))
    rownames(grouping_counts) <- grouping_counts$Var1
    not_enough_cells <- rownames(grouping_counts[grouping_counts$Freq < 300,])

    # define a mapping to iteratively combine goupings for a given sample
    test_epithelial_cell_types <- c("Enterocytes",
         "Enterocyte Progenitors",
         "TA1",
         "TA2",
         "Stem",
         "Secretory TA",
         "Goblet")

    combine_with_epithelial_cell_types <- c("Immature Enterocytes",
         "Immature Enterocytes",
         "TA2",
         "Stem",
         "TA2",
         "Immature Goblet",
         "Immature Goblet")

    # iterate through samples and merge cell types as necessary, as indicated above
    for (i in 1:length(test_epithelial_cell_types)){
      for (j in sample_names_only){
        if (paste(j,test_epithelial_cell_types[i],sep = "-") %in% not_enough_cells){
          cell_types[cell_types$CellType == paste(j,test_epithelial_cell_types[i],sep = "-"),] <- paste(j,combine_with_epithelial_cell_types[i],sep = "-")
        }
      }
    }

    # merge the less common enteroendocrine cells and best4+ enterocytes with the more common sample types
    cell_types[cell_types$CellType %in% "Normal-Enteroendocrine",] <- "Unaffected-Enteroendocrine" 
    cell_types[cell_types$CellType %in% "Adenocarcinoma-Enteroendocrine",] <- "Polyp-Enteroendocrine" 
    cell_types[cell_types$CellType %in% "Adenocarcinoma-Best4+ Enterocytes",] <- "Polyp-Best4+ Enterocytes" 

    # redefine the groupings to see where we still need to combine things
    grouping_counts <- data.frame(table(cell_types))
    rownames(grouping_counts) <- grouping_counts$Var1
    not_enough_cells <- rownames(grouping_counts[grouping_counts$Freq < 300,])

    # Now combine remaining rare cell groupings by patient
    cell_types[cell_types$CellType %in% (not_enough_cells[grepl("*Immature Enterocytes", not_enough_cells) & grepl("A002*", not_enough_cells)]),] <- "A002-Immature Enterocytes" 
    cell_types[cell_types$CellType %in% (not_enough_cells[grepl("*Immature Enterocytes", not_enough_cells) & grepl("F*", not_enough_cells)]),] <- "F-Immature Enterocytes" 
    cell_types[cell_types$CellType %in% (not_enough_cells[grepl("*Immature Enterocytes", not_enough_cells) & grepl("A014*", not_enough_cells)]),] <- "A014-Immature Enterocytes"
    cell_types[cell_types$CellType %in% (not_enough_cells[grepl("*Immature Goblet", not_enough_cells) & grepl("A014_Una*", not_enough_cells)]),] <- "A002-Unaffected-Immature Goblet" 

    # And merge the last few that still have less than 200 into a grouping from the same patient
    cell_types[cell_types$CellType == "F091-Immature Goblet",] <- "F034-Immature Goblet" 
    cell_types[cell_types$CellType == "F091-TA2",] <- "F034-TA2" 
    cell_types[cell_types$CellType == "A015-C-001-Immature Goblet",] <- "A015-C-002-Immature Goblet" 
    cell_types[cell_types$CellType == "A015-C-001-TA2",] <- "A015-C-002-TA2"
    cell_types[cell_types$CellType == "CRC-1-8810-Immature Goblet",] <- "A001-C-007-Immature Goblet"
    cell_types[cell_types$CellType == "CRC-2-15564-Immature Goblet",] <- "A001-C-007-Immature Goblet"
    cell_types[cell_types$CellType == "CRC-3-11773-Immature Goblet",] <- "A001-C-007-Immature Goblet"
    cell_types[cell_types$CellType == "A014-C-008-Immature Goblet",] <- "A014-C-108-Immature Goblet"
    cell_types[cell_types$CellType == "A002-C-021-Immature Goblet",] <- "A002-C-016-Immature Goblet"

    # Now add the peak calling groups to the project
    proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = cell_types$CellType, cells = rownames(cell_types), name = "PeakCallingGroup", force = TRUE)

    #Create Group Coverage Files that can be used for downstream analysis
    proj_stromal <- addGroupCoverages(ArchRProj = proj_stromal, groupBy = "PeakCallingGroup", force = TRUE)

    # Find macs2
    pathToMacs2 <- findMacs2()

    #Call Reproducible Peaks w/ Macs2
    proj_stromal <- addReproduciblePeakSet(
        ArchRProj = proj_stromal, groupBy = "PeakCallingGroup", force = TRUE, 
        pathToMacs2 = pathToMacs2
    )

    #Add Peak Matrix
    proj_stromal <- addPeakMatrix(ArchRProj = proj_stromal, force = TRUE)

    saveArchRProject(ArchRProj = proj_stromal, outputDirectory = new_project_save_name, load = FALSE, overwrite = FALSE)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 4) pairwise comparisons of polyps to true normals
if (4 %in% execute_steps){
    # Add motif annotations if needed
    if("Motif" %ni% names(proj_stromal@peakAnnotation)){
      proj_stromal <- addMotifAnnotations(ArchRProj = proj_stromal, motifSet = "cisbp", name = "Motif")
    }
    system("wget https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra-Human-Motifs.rds")
    motifPWMs <- readRDS("Vierstra-Human-Motifs.rds")
    proj_stromal <- addMotifAnnotations(proj_stromal, motifPWMs = motifPWMs, name = "Vierstra")

    cellTypes <- unique(getCellColData(proj_stromal, "CellType")$CellType)

    # Add CellTypeDiseaseState if not present
    proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = paste(getCellColData(ArchRProj = proj_stromal, select = "DiseaseState")$DiseaseState, getCellColData(ArchRProj = proj_stromal, select = "CellType")$CellType, sep ='-'), cells = rownames(getCellColData(ArchRProj = proj_stromal, select = "CellType")), name = "CellTypeDiseaseState", force = TRUE)

    # Read metadata file
    metadata <- read.table("./hubmap_htan_metadata_atac_and_rna_final.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)

    # Set polyp names and matched normals for pairwise comparisons
    polypNames <-  metadata[metadata$DiseaseState != "Normal",]$SimplifiedSampleName

    # Get sample names
    samples <- getCellColData(ArchRProj = proj_stromal, select = "SimplifiedSampleName")

    # Simplify cell types
    NewCellType <- getCellColData(ArchRProj = proj_stromal, select = "CellType")
    proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = paste(samples$SimplifiedSampleName, NewCellType$CellType, sep ='-'), cells = rownames(getCellColData(ArchRProj = proj_stromal, select = "CellType")), name = "CellTypeSample", force = TRUE)

    cellTypesToCompare <- c("Stem", "TA2")
    cellTypesToCompare <- c("Stem")
    for (j in 1:length(cellTypesToCompare)){
        typeA <- cellTypesToCompare[j]
        typeB <- cellTypesToCompare[j]
        for (i in 1:length(polypNames)){
          sampleNameA <- polypNames[i]
          sampleNameB <- "NormalColon"
          if (sum(proj_stromal$CellTypeSample == paste0(sampleNameA, "-", typeA))>200){
            markerTest <- getMarkerFeatures(
              ArchRProj = proj_stromal, 
              useMatrix = "PeakMatrix",
              groupBy = "CellTypeSample",
              testMethod = "wilcoxon",
              bias = c("TSSEnrichment", "log10(nFrags)"),
              useGroups = paste0(sampleNameA, "-", typeA),
              bgdGroups = paste0(sampleNameB, "-", typeB)
            )
            motifsUp <- peakAnnoEnrichment(seMarker = markerTest,ArchRProj = proj_stromal,peakAnnotation = "Motif",cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
            motifsDo <- peakAnnoEnrichment(seMarker = markerTest,ArchRProj = proj_stromal,peakAnnotation = "Motif",cutOff = "FDR <= 0.1 & Log2FC <= -0.5")

            # Add to an object to save all of the differentials and motif enrichments within those differentials
            new_additionUp <- assays(motifsUp)$mlog10Padj
            colnames(new_additionUp) = paste0(colnames(new_additionUp), "Stem")
            new_additionDown <- assays(motifsDo)$mlog10Padj
            colnames(new_additionDown) = paste0(colnames(new_additionDown), "Stem")
            if (i==1){
              save_motifs_up <- new_additionUp
              save_motifs_down <- new_additionDown
              metadata(markerTest)$Params$ArchRProj <- NULL
              fullMarkerTest <- markerTest
            }
            if (i>1){
              save_motifs_up <- cbind(new_additionUp, save_motifs_up)
              save_motifs_down <- cbind(new_additionDown, save_motifs_down)
              metadata(markerTest)$Params$ArchRProj <- NULL
              fullMarkerTest <- cbind(fullMarkerTest, markerTest)
            }
          }
        }

        saveRDS(save_motifs_up, paste0("./individual_sample_motif_enrichments/pairwise_polyp_to_normal_comparisons_up_", typeA, ".rds"))
        saveRDS(save_motifs_down, paste0("./individual_sample_motif_enrichments/pairwise_polyp_to_normal_comparisons_down_", typeA, ".rds"))
        saveRDS(fullMarkerTest, paste0("./individual_sample_motif_enrichments/pairwise_polyp_to_normal_comparisons_markerTestFile_", typeA, ".rds"))
        
        peaks_down <- data.frame(colSums((assays(fullMarkerTest)$Log2FC < -0.5) & (assays(fullMarkerTest)$FDR<0.1)))
        colnames(peaks_down) <- "num_diff"
        peaks_down$sample <- colnames(fullMarkerTest)
        peaks_down$change <- rep("down", length(colnames(fullMarkerTest)))

        peaks_up <- data.frame(colSums((assays(fullMarkerTest)$Log2FC > 0.5) & (assays(fullMarkerTest)$FDR<0.1)))
        colnames(peaks_up) <- "num_diff"
        peaks_up$sample <- colnames(fullMarkerTest)
        peaks_up$change <- rep("up", length(colnames(fullMarkerTest)))

        peaks_diff <- rbind(peaks_down, peaks_up)

        p <- ggplot(peaks_diff, aes(y=num_diff, x=sample, fill = change)) + 
            geom_bar(position="stack", stat="identity") + theme_ArchR() + 
          ylab("Fraction Cells") +
          xlab("Sample") + theme(axis.text.x = element_text(angle = 90)) +
          scale_fill_manual("legend", values = c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC"))
        plotPDF(p, name = paste(paste("Peaks-change-stem-relative-to-normal-stem", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_stromal, addDOC = FALSE, width = 15, height = 15)
    }
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 5) pairwise comparisons of polyps to all unaffecteds
if (5 %in% execute_steps){
    # Add CellTypeDiseaseState if not present
    proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = paste(getCellColData(ArchRProj = proj_stromal, select = "DiseaseState")$DiseaseState, getCellColData(ArchRProj = proj_stromal, select = "CellType")$CellType, sep ='-'), cells = rownames(getCellColData(ArchRProj = proj_stromal, select = "CellType")), name = "CellTypeDiseaseState", force = TRUE)

    # Read metadata file
    metadata <- read.table("./hubmap_htan_metadata_atac_and_rna_final.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)

    # Set polyp names and matched normals for pairwise comparisons
    polypNames <-  metadata[metadata$DiseaseState != "Normal",]$DifferentialGroup
    polypNames <-  polypNames[!grepl("Unaffected", polypNames)]

    # Get sample names
    cellsNamesToAdd <- c()
    annotationToAdd <- c()
    for (i in 1:dim(metadata)[1]){
      idxSample <- BiocGenerics::which(getCellColData(proj_stromal, "Sample") %in% metadata[i,"Sample"])
      cellsSample <- proj_stromal$cellNames[idxSample[["Sample"]]]
      cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
      annotationToAdd <- append(annotationToAdd, rep(metadata[i,"DifferentialGroup"], length(cellsSample)))
    }
    proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = paste0(annotationToAdd), cells = paste0(cellsNamesToAdd), name = "DifferentialGroup", force = TRUE)

    samples <- getCellColData(ArchRProj = proj_stromal, select = "DifferentialGroup")

    # Simplify cell types
    NewCellType <- getCellColData(ArchRProj = proj_stromal, select = "CellType")
    proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = paste(samples$DifferentialGroup, NewCellType$CellType, sep ='-'), cells = rownames(getCellColData(ArchRProj = proj_stromal, select = "CellType")), name = "CellTypeSample", force = TRUE)

    cellTypesToCompare <- c("Stem", "TA2")
    cellTypesToCompare <- c("Stem")
    for (j in 1:length(cellTypesToCompare)){
        typeA <- cellTypesToCompare[j]
        typeB <- cellTypesToCompare[j]
        for (i in 1:length(polypNames)){
          sampleNameA <- polypNames[i]
          sampleNameB <- "Unaffected"
          if (sum(proj_stromal$CellTypeSample == paste0(sampleNameA, "-", typeA))>200){
            markerTest <- getMarkerFeatures(
              ArchRProj = proj_stromal, 
              useMatrix = "PeakMatrix",
              groupBy = "CellTypeSample",
              testMethod = "wilcoxon",
              bias = c("TSSEnrichment", "log10(nFrags)"),
              useGroups = paste0(sampleNameA, "-", typeA),
              bgdGroups = paste0(sampleNameB, "-", typeB)
            )

            motifsUp <- peakAnnoEnrichment(seMarker = markerTest,ArchRProj = proj_stromal,peakAnnotation = "Motif",cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
            motifsUpVierstra <- peakAnnoEnrichment(seMarker = markerTest,ArchRProj = proj_stromal,peakAnnotation = "Vierstra",cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
            motifsDo <- peakAnnoEnrichment(seMarker = markerTest,ArchRProj = proj_stromal,peakAnnotation = "Motif",cutOff = "FDR <= 0.1 & Log2FC <= -0.5")
            motifsDoVierstra <- peakAnnoEnrichment(seMarker = markerTest,ArchRProj = proj_stromal,peakAnnotation = "Vierstra",cutOff = "FDR <= 0.1 & Log2FC <= -0.5")

            # Add to an object to save all of the differentials and motif enrichments within those differentials
            new_additionUp <- assays(motifsUp)$mlog10Padj
            colnames(new_additionUp) = paste0(colnames(new_additionUp), "Stem")
            new_additionDown <- assays(motifsDo)$mlog10Padj
            colnames(new_additionDown) = paste0(colnames(new_additionDown), "Stem")
            new_additionUpVierstra <- assays(motifsUpVierstra)$mlog10Padj
            colnames(new_additionUpVierstra) = paste0(colnames(new_additionUpVierstra), "Stem")
            new_additionDownVierstra <- assays(motifsDoVierstra)$mlog10Padj
            colnames(new_additionDownVierstra) = paste0(colnames(new_additionDownVierstra), "Stem")
            if (i==1){
              save_motifs_up <- new_additionUp
              save_motifs_down <- new_additionDown
              save_motifs_upVierstra <- new_additionUpVierstra
              save_motifs_downVierstra <- new_additionDownVierstra
              metadata(markerTest)$Params$ArchRProj <- NULL
              fullMarkerTest <- markerTest
            }
            if (i>1){
              save_motifs_up <- cbind(new_additionUp, save_motifs_up)
              save_motifs_down <- cbind(new_additionDown, save_motifs_down)
              save_motifs_upVierstra <- cbind(new_additionUpVierstra, save_motifs_upVierstra)
              save_motifs_downVierstra <- cbind(new_additionDownVierstra, save_motifs_downVierstra)
              metadata(markerTest)$Params$ArchRProj <- NULL
              fullMarkerTest <- cbind(fullMarkerTest, markerTest)
            }
          }
        }

        saveRDS(save_motifs_up, paste0("./individual_sample_motif_enrichments/pairwise_polyp_to_allunaffected_comparisons_up_", typeA, ".rds"))
        saveRDS(save_motifs_down, paste0("./individual_sample_motif_enrichments/pairwise_polyp_to_allunaffected_comparisons_down_", typeA, ".rds"))
        saveRDS(save_motifs_upVierstra, paste0("./individual_sample_motif_enrichments/pairwise_polyp_to_allunaffected_comparisons_up_Vierstra_", typeA, ".rds"))
        saveRDS(save_motifs_downVierstra, paste0("./individual_sample_motif_enrichments/pairwise_polyp_to_allunaffected_comparisons_down_Vierstra_", typeA, ".rds"))
        saveRDS(fullMarkerTest, paste0("./individual_sample_motif_enrichments/pairwise_polyp_to_allunaffected_comparisons_markerTestFile_", typeA, ".rds"))
        
        peaks_down <- data.frame(colSums((assays(fullMarkerTest)$Log2FC < -0.5) & (assays(fullMarkerTest)$FDR<0.1)))
        colnames(peaks_down) <- "num_diff"
        peaks_down$sample <- colnames(fullMarkerTest)
        peaks_down$change <- rep("down", length(colnames(fullMarkerTest)))

        peaks_up <- data.frame(colSums((assays(fullMarkerTest)$Log2FC > 0.5) & (assays(fullMarkerTest)$FDR<0.1)))
        colnames(peaks_up) <- "num_diff"
        peaks_up$sample <- colnames(fullMarkerTest)
        peaks_up$change <- rep("up", length(colnames(fullMarkerTest)))

        peaks_diff <- rbind(peaks_down, peaks_up)

        p <- ggplot(peaks_diff, aes(y=num_diff, x=sample, fill = change)) + 
            geom_bar(position="stack", stat="identity") + theme_ArchR() + 
          ylab("Fraction Cells") +
          xlab("Sample") + theme(axis.text.x = element_text(angle = 90)) +
          scale_fill_manual("legend", values = c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC"))
        plotPDF(p, name = paste(paste("Peaks-change-stem-relative-to-allunaffected-stem", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_stromal, addDOC = FALSE, width = 15, height = 15)
    }
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 7) compute malignancy pseudotimes and plot results of differential analyses
if (7 %in% execute_steps){
  typeA <- "Stem"
  #typeA <- "TA2"
  backgrounds <- c("normal", "allunaffected")
  for (background in backgrounds){
    fullMarkerTest <- readRDS(paste0("./individual_sample_motif_enrichments/pairwise_polyp_to_", background, "_comparisons_markerTestFile_", typeA, ".rds"))
    # remove samples w/ less than 250 cells and duplicated samples from full test object
    meta <- getCellColData(proj_stromal)
    more_than_250 <- names(which(table(meta[meta$CellType == typeA,]$SimplifiedSampleName)>250))
    # if you want to include the unaffected samples with some dysplasia, uncomment the following line
    #more_than_250 <- names(which(table(meta[meta$CellType == typeA,]$DifferentialGroup)>250))
    fullMarkerTest250 <- fullMarkerTest[,colnames(fullMarkerTest) %in% paste0(more_than_250, "-", typeA)]
    fullMarkerTest250 <- fullMarkerTest250[,!duplicated(colnames(fullMarkerTest250))]

    # set NA values to the nonsignficant value of 1
    assays(fullMarkerTest250)$FDR[is.na(assays(fullMarkerTest250)$FDR)] <- 1.0
    # set criteria for peaks to plot
    fullMarkerTest250 <- fullMarkerTest250[rowSums(assays(fullMarkerTest250)$FDR<=0.05 & assays(fullMarkerTest250)$Log2FC<=-1.5)>1 |
                              rowSums(assays(fullMarkerTest250)$FDR<=0.05 & assays(fullMarkerTest250)$Log2FC>=1.5)>1,]

    # now compute the pcs on the log2fc of the significant peaks
    pcs <- prcomp(assays(fullMarkerTest250)$Log2FC)
    summary(pcs) # to get proportions of variance
    pc_df <- data.frame(pcs$rotation)
    pc_df$sample <- rownames(pc_df)

    # add the metadata for plotting
    sample_names <- substr(rownames(pc_df),1,nchar(rownames(pc_df))-(nchar(typeA)+1))
    metadata <- read.table("./hubmap_htan_metadata_atac_and_rna_final.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)
    
    cM <- as.data.frame(confusionMatrix(paste0(proj_stromal$CellType), paste0(proj_stromal$SimplifiedSampleName)))
    cM <- t(cM)/colSums(cM)
    f_stem <- cM[,"Stem"]
    percent_stem <- c()
    for (sample in metadata$SimplifiedSampleName){
      if (sample %in% names(f_stem)){
        percent_stem <- c(percent_stem, f_stem[sample])
      } else {
        percent_stem <- c(percent_stem, NA)
      }
    }
    metadata$percent_stem <- percent_stem
    
    metadata <- metadata[,c("SimplifiedSampleName","NearestUnaffected","Donor","FAP","DiseaseState", "Size","Location","Polyp.Type","Dysplasia","HGD","LGD"
      ,"Percent_HGD","Percent_LGD","Percent_sample_wo_neoplasia_overall","Percent_sample_w_cancer_overall","Percent_sample_w_any_degree_dysplasia_overall"
      ,"Pathologist0_PolypType","Pathologist0_HGD","Pathologist0_LGD","Pathologist0_Percent_of_neoplastic_cell_nuclei_as_a_total_of_all_cell_nuclei_in_tumor_ONLY_area"
      ,"Pathologist0_Percent_nondysplastic_stroma","Pathologist0_Percent_Entire_Tissue_Dysplastic","percent_stem" )]
    metadata <- metadata[metadata$SimplifiedSampleName %in% sample_names,]
    metadata <- metadata[!duplicated(metadata),]
    # remove specific partial duplicates
    if (background == "allunaffected"){
      metadata <- metadata[rownames(metadata) %ni% c(),]
    }
    if (background == "normal"){
      metadata <- metadata[rownames(metadata) %ni% c(12,24,58),] # note this may change based on the metadata file you are using
    }
    rownames(metadata) <- metadata$SimplifiedSampleName
    metadata <- metadata[sample_names,]

    #rownames(pc_df) <- metadata
    pc_df <- cbind(pc_df, metadata)
    scalef <- 1
    pc_df$PC1 <- pc_df$PC1*scalef

    # plot the first two pcs
    p <- ggplot(pc_df, aes(x=PC1, y=PC2, color=DiseaseState)) +
      geom_point(size=2) + scale_color_manual(values=c("#D51F26", "#89288F", "#208A42")) + theme_ArchR()
    plotPDF(p, name = paste0("All_Polyp_", typeA, "_vs_", background, "_", typeA, "_Peak_Change_Heatmap_250_cell_cutoff-l2fc2cutoff-atleast2_affected_only_PCA.pdf"), width = 5, height = 5, ArchRProj = proj_stromal, addDOC = FALSE)

    # fit a spline
    require(splines)
    if (background == "normal"){
      fit<-lm(PC2 ~ bs(PC1,knots = c(-0.11)),data = pc_df )
      age.grid<-seq(from=-2000*scalef, to = -300*scalef)/10000
      if (typeA == "TA2"){
        fit<-lm(PC2 ~ bs(PC1,knots = c(0.15)),data = pc_df )
        age.grid<-seq(from=0*scalef, to = 2500*scalef)/10000
      }
    } else {
      fit<-lm(PC2 ~ bs(PC1,knots = c(0.1)),data = pc_df )
      age.grid<-seq(from=-1000*scalef, to = 4000*scalef)/10000
    }

    splinefit <- data.frame(age.grid,predict(fit,newdata = list(PC1=age.grid)))
    colnames(splinefit) <- c("PC1", "PC2")

    # plot it
    p <- ggplot(pc_df, aes(x=PC1, y=PC2, color=DiseaseState)) +
      geom_point(size=2) + scale_color_manual(values=c("#D51F26", "#89288F", "#208A42")) + 
      xlab(paste0("PC1 (",100*summary(pcs)$importance["Proportion of Variance","PC1"], "%)")) + ylab(paste0("PC2 (",100*summary(pcs)$importance["Proportion of Variance","PC2"], "%)")) + theme_ArchR()
    p <- p + geom_line(data=splinefit, colour="#CC0000")
    plotPDF(p, name = paste0("All_Polyp_", typeA, "_vs_", background, "_", typeA, "_Peak_Change_Heatmap_250_cell_cutoff-l2fc2cutoff-atleast2_affected_only_PCA.pdf"), width = 5, height = 5, ArchRProj = proj_stromal, addDOC = FALSE)

    # get x values of nearest points--this will order from "least to most cancerous"
    x_vals <- c()
    for (i in 1:length(pc_df[,1])){
      x_vals <- c(x_vals,splinefit$PC1[which.min(sqrt((pc_df$PC1[i]-splinefit$PC1)**2+(pc_df$PC2[i]-splinefit$PC2)**2))])
    }

    pc_df$nearest_spline_x_vals <- x_vals
    pc_df <- pc_df[order(pc_df$nearest_spline_x_vals),]

    saveRDS(pc_df, paste0("./differential_PCA_pseudotime/pc_df_", typeA, "_vs_", background,".rds"))
    pc_df <- readRDS(paste0("./differential_PCA_pseudotime/pc_df_Stem_vs_", "normal",".rds"))
    #pc_df <- pc_df[order(pc_df$nearest_spline_x_vals, decreasing = TRUE),]

    p <- ggplot(pc_df, aes(x=nearest_spline_x_vals, y=Percent_sample_w_any_degree_dysplasia_overall, color=DiseaseState)) +
      geom_point(size=2) + scale_color_manual(values=c("#D51F26", "#89288F", "#208A42")) + theme_ArchR()
    plotPDF(p, name = paste0("All_Polyp_", typeA, "_cells_vs_", background, "_", typeA, "cells_continuum_vs_any_dysplasia.pdf"), width = 5, height = 5, ArchRProj = proj_stromal, addDOC = FALSE)

  }

  backgrounds <- c("normal", "allunaffected")
  for (background in backgrounds){
    fullMarkerTest <- readRDS(paste0("./individual_sample_motif_enrichments/pairwise_polyp_to_", background, "_comparisons_markerTestFile_", typeA, ".rds"))

    # remove samples w/ less than 250 cells and duplicated samples from full test object
    meta <- getCellColData(proj_stromal)
    more_than_250 <- names(which(table(meta[meta$CellType == typeA,]$SimplifiedSampleName)>250))
    # if you want to include the unaffected samples with some dysplasia, uncomment the following line
    #more_than_250 <- names(which(table(meta[meta$CellType == typeA,]$DifferentialGroup)>250))
    fullMarkerTest250 <- fullMarkerTest[,colnames(fullMarkerTest) %in% paste0(more_than_250, "-", typeA)]
    fullMarkerTest250 <- fullMarkerTest250[,!duplicated(colnames(fullMarkerTest250))]

    # set NA values to the nonsignficant value of 1
    assays(fullMarkerTest250)$FDR[is.na(assays(fullMarkerTest250)$FDR)] <- 1.0
    # set criteria for peaks to plot
    fullMarkerTest250 <- fullMarkerTest250[rowSums(assays(fullMarkerTest250)$FDR<=0.05 & assays(fullMarkerTest250)$Log2FC<=-1.5)>1 |
                              rowSums(assays(fullMarkerTest250)$FDR<=0.05 & assays(fullMarkerTest250)$Log2FC>=1.5)>1,]
   
    pc_df <- readRDS(paste0("./differential_PCA_pseudotime/pc_df_Stem_vs_", "normal",".rds"))

    pc_df <- pc_df[rownames(pc_df) %in% colnames(fullMarkerTest250),]
    
    # plot a bar graph of number of differential peaks for all samples
    peaks_down <- data.frame(colSums((assays(fullMarkerTest250)$Log2FC <= -1.5) & (assays(fullMarkerTest250)$FDR<=0.05)))
    colnames(peaks_down) <- "num_diff"
    peaks_down$sample <- colnames(fullMarkerTest250)
    peaks_down$change <- rep("down", length(colnames(fullMarkerTest250)))

    peaks_up <- data.frame(colSums((assays(fullMarkerTest250)$Log2FC >= 1.5) & (assays(fullMarkerTest250)$FDR<=0.05)))
    colnames(peaks_up) <- "num_diff"
    peaks_up$sample <- colnames(fullMarkerTest250)
    peaks_up$change <- rep("up", length(colnames(fullMarkerTest250)))

    peaks_diff <- rbind(peaks_down, peaks_up)
    peaks_diff <- peaks_diff[!duplicated(peaks_diff),]

    peaks_diff <- peaks_diff[order(sapply(peaks_diff$sample, function(x) which(x == rownames(pc_df)))), ]
    peaks_diff$sample <- factor(peaks_diff$sample, levels = rownames(pc_df))

    p <- ggplot(peaks_diff, aes(x=sample, y=num_diff, fill = change)) + 
        geom_bar(position="stack", stat="identity") + theme_ArchR() + 
      ylab("Fraction Cells") +
      xlab("Sample") + theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_manual("legend", values = c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC"))
    plotPDF(p, name = paste(paste("Peaks-change-", typeA, "-relative-to-", background, "-", typeA, "-greater-than-250-cells-l2fc1p5-fdr0p05-", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_stromal, addDOC = FALSE, width = 15, height = 15)

  # make a heatmap where they are ordered by how cancerous they are
  set.seed(10)
  test <- kmeans(assays(fullMarkerTest250)$Log2FC, 10, iter.max = 500)#, algorithm = "Lloyd")

  #save the cluster identities of the peaks
  rowData(fullMarkerTest250)$cluster <- test$cluster
  saveRDS(rowData(fullMarkerTest250), paste0("./kmeans_peak_clusters_stem_vs_", background,".rds"))

  # now find motifs enriched in different clusters

  # iterate through clusters to find enriched motifs
  for (i in 1:10){
    clusters <- test$cluster
    clusters <- clusters[clusters == i]
    clusterMarkerTest <- fullMarkerTest[,paste0("A001-C-002-", typeA)]
    #clusterMarkerTest <- fullMarkerTest250[names(clusters),"A001-C-002-Stem"]
    # create dummy test where only the rows in the cluster are significant
    assays(clusterMarkerTest)$Log2FC[,] <- 0
    assays(clusterMarkerTest)$Log2FC[names(clusters),] <- 10000
    motifsUp <- peakAnnoEnrichment(
      seMarker = clusterMarkerTest,
      ArchRProj = proj_stromal,
      peakAnnotation = "Vierstra",
      cutOff = "Log2FC >9999"
    )

    df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
    df <- df[order(df$mlog10Padj, decreasing = TRUE),]
    df$rank <- seq_len(nrow(df))

    # ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    #   geom_point(size = 1) +
    #   ggrepel::geom_label_repel(data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),size = 1.5,nudge_x = 2,color = "black") + 
    #   theme_ArchR() + 
    #   ylab("-log10(P-adj) Motif Enrichment") + 
    #   xlab("Rank Sorted TFs Enriched") +
    #   scale_color_gradientn(colors = paletteContinuous(set = "comet"))

    # plotPDF(ggUp, name = paste0("All_Polyp_", typeA, "_vs_", background, "_", typeA, "_Peak_Change_Heatmap_250_cell_cutoff-l2fc2cutoff-atleast2_cluster", i), width = 5, height = 5, ArchRProj = proj_stromal, addDOC = FALSE)

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
  df <- df[apply(df, 1, function(x) !all(x==0)),]
  df <- df[apply(df, 1, function(x) sum(x>20)>0),]
  paletteLength <- 256
  myBreaks <- c(seq(0, 50, length.out=ceiling(paletteLength/2) + 1), 
                seq(50.1, 200, length.out=floor(paletteLength/2)))
  p <- pheatmap::pheatmap(
      mat = as.matrix(df), 
      color = paletteContinuous("whiteBlue"), 
      border_color = "black", breaks = myBreaks, cluster_cols = TRUE
  )
  plotPDF(p, name = paste0("Motif-Enrichment_Differential_Peak_vs-background-", background, "-in-kmeans-clusters-", typeA), width = 20, height = 30,  ArchRProj = proj_stromal, addDOC = FALSE)

  motif_names <- names(apply(save_motifs_up, 1, max)[order(apply(save_motifs_up, 1, max), decreasing = TRUE)])
  motif_names_df <- t(data.frame(strsplit(motif_names, ":")))
  #motif_names_df <- motif_names_df[order(as.numeric(motif_names_df[,3])),]
  #non-redundant df
  df <- DataFrame(save_motifs_up[motif_names[order(as.numeric(motif_names_df[,3]))],])
  df <- df[apply(df, 1, function(x) !all(x==0)),]
  df <- df[apply(df, 1, function(x) sum(x>20)>0),]
  paletteLength <- 256
  myBreaks <- c(seq(0, 50, length.out=ceiling(paletteLength/2) + 1), 
                seq(50.1, 200, length.out=floor(paletteLength/2)))
  p <- pheatmap::pheatmap(
      mat = as.matrix(df), 
      color = paletteContinuous("whiteBlue"), 
      border_color = "black", breaks = myBreaks, cluster_cols = TRUE, cluster_rows = FALSE
  )
  plotPDF(p, name = paste0("Motif-Enrichment_Differential_Peak_", background, "_order_motif_clusters-", typeA), width = 20, height = 30,  ArchRProj = proj_stromal, addDOC = FALSE)


  # test different n for kmeans
  set.seed(10)
  nclust <- 15
  test <- kmeans(assays(fullMarkerTest250)$Log2FC, nclust, iter.max = 500)#, algorithm = "Lloyd")

  #save the cluster identities of the peaks
  rowData(fullMarkerTest250)$cluster <- test$cluster

  # now find motifs enriched in different clusters
  # iterate through clusters to find enriched motifs
  for (i in 1:nclust){
    clusters <- test$cluster
    clusters <- clusters[clusters == i]
    clusterMarkerTest <- fullMarkerTest[,paste0("A001-C-002-", typeA)]
    #clusterMarkerTest <- fullMarkerTest250[names(clusters),"A001-C-002-Stem"]
    # create dummy test where only the rows in the cluster are significant
    assays(clusterMarkerTest)$Log2FC[,] <- 0
    assays(clusterMarkerTest)$Log2FC[names(clusters),] <- 10000
    motifsUp <- peakAnnoEnrichment(
      seMarker = clusterMarkerTest,
      ArchRProj = proj_stromal,
      peakAnnotation = "Vierstra",
      cutOff = "Log2FC >9999"
    )

    df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
    df <- df[order(df$mlog10Padj, decreasing = TRUE),]
    df$rank <- seq_len(nrow(df))

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
  df <- df[apply(df, 1, function(x) !all(x==0)),]
  df <- df[apply(df, 1, function(x) sum(x>20)>0),]
  paletteLength <- 256
  myBreaks <- c(seq(0, 50, length.out=ceiling(paletteLength/2) + 1), 
                seq(50.1, 200, length.out=floor(paletteLength/2)))
  p <- pheatmap::pheatmap(
      mat = as.matrix(df), 
      color = paletteContinuous("whiteBlue"), 
      border_color = "black", breaks = myBreaks, cluster_cols = TRUE
  )
  plotPDF(p, name = paste0(nclust, "-kmeans-clusters-", typeA, "Motif-Enrichment_Differential_Peak_vs-background-", background), width = 20, height = 30,  ArchRProj = proj_stromal, addDOC = FALSE)

  motif_names <- names(apply(save_motifs_up, 1, max)[order(apply(save_motifs_up, 1, max), decreasing = TRUE)])
  motif_names_df <- t(data.frame(strsplit(motif_names, ":")))
  #non-redundant df
  df <- DataFrame(save_motifs_up[motif_names[order(as.numeric(motif_names_df[,3]))],])
  df <- df[apply(df, 1, function(x) !all(x==0)),]
  df <- df[apply(df, 1, function(x) sum(x>20)>0),]
  paletteLength <- 256
  myBreaks <- c(seq(0, 50, length.out=ceiling(paletteLength/2) + 1), 
                seq(50.1, 200, length.out=floor(paletteLength/2)))
  p <- pheatmap::pheatmap(
      mat = as.matrix(df), 
      color = paletteContinuous("whiteBlue"), 
      border_color = "black", breaks = myBreaks, cluster_cols = TRUE, cluster_rows = FALSE
  )
  plotPDF(p, name = paste0(nclust, "-kmeans-clusters-", "Motif-Enrichment_Differential_Peak_", background, "_order_motif_clusters-", typeA), width = 20, height = 30,  ArchRProj = proj_stromal, addDOC = FALSE)

  # sample vs celltype dist
  pc_df <- readRDS(paste0("./differential_PCA_pseudotime/pc_df_Stem_vs_", "normal",".rds"))

  cM <- confusionMatrix(paste0(getCellColData(proj_stromal)$CellType), paste0(getCellColData(proj_stromal)$SimplifiedSampleName))
  cM <- data.frame(t(cM) / Matrix::rowSums(t(cM)))
  samplenames <- c()
  cellnames <- c()
  values <- c()
  diseaseStates <- c()
  for (sample in substr(rownames(pc_df),1,nchar(rownames(pc_df))-5)){
      diseaseState <- unique(getCellColData(proj_stromal)[getCellColData(proj_stromal)$SimplifiedSampleName == sample,]$DiseaseState)
      for (celltype in colnames(cM)){
          samplenames <- c(samplenames, sample)
          cellnames <- c(cellnames, celltype)
          diseaseStates <- c(diseaseStates, diseaseState)
          values <- c(values, cM[sample, celltype])
      }
  }
  data <- data.frame(samplenames,cellnames,values, diseaseStates)
  #data <- data[order(samplenames), ]
  #data <- data[order(sapply(data$diseaseStates, function(x) which(x == c("Normal", "Unaffected", "Polyp", "Adenocarcinoma")))), ]

  data$samplenames <- factor(data$samplenames, levels = unique(data$samplenames))

  # Set appropriate order
  new_levels <-  c("Stem",
     "TA2",
     "TA1",
     "Enterocyte.Progenitors",
     "Immature.Enterocytes",
     "Enterocytes",
     "Secretory.TA",
     "Immature.Goblet",
     "Goblet",
     "Enteroendocrine",
     "Best4..Enterocytes")  

  data$cellnames <- factor(data$cellnames, levels = new_levels)

  test <- data.frame(table(getCellColData(proj_stromal)$Sample))
  rownames(test) <- test$Var1


  p <- ggplot(data, aes(fill=cellnames, y=values, x=samplenames)) + 
      geom_bar(position="stack", stat="identity") + theme_ArchR() + 
    ylab("Fraction Cells") +
    xlab("Sample") + theme(axis.text.x = element_text(angle = 90)) +
    scale_x_discrete(labels=paste0(levels(data$samplenames), " (n = ",  test[levels(data$samplenames),]$Freq, ")")) +
    scale_fill_manual("legend", values = c("Best4..Enterocytes"="#D51F26", "Enterocyte.Progenitors" = "#208A42", "Enterocytes" = "#89288F", 
      "Enteroendocrine" = "#F47D2B", "Goblet" = "#FEE500", "Immature.Enterocytes" = "#8A9FD1", "Immature.Goblet" = "#C06CAB", "Secretory.TA" = "#D24B27", "Stem" = "#D8A767", "TA1" = "#90D5E4", "TA2" = "#89C75F"))
  plotPDF(p, name = paste0("Sample-CellType-Bar-Fraction_Sample_Epithelial_continuuum-ordered", ".pdf"), ArchRProj = proj_stromal, addDOC = FALSE, width = 25, height = 20)


  pc_df$nearest_spline_x_vals <- pc_df$nearest_spline_x_vals*-1
  cM <- as.data.frame(confusionMatrix(paste0(proj_stromal$CellType), paste0(proj_stromal$SimplifiedSampleName)))
  cM <- t(cM)/colSums(cM)
  for (celltype in colnames(cM)){
    f_celltype <- cM[,celltype]
    percent_celltype <- c()
    for (sample in pc_df$SimplifiedSampleName){
      if (sample %in% names(f_celltype)){
        percent_celltype <- c(percent_celltype, f_celltype[sample])
      } else {
        percent_celltype <- c(percent_celltype, NA)
      }
    }
    pc_df$percent_celltype <- percent_celltype

    # require(splines)
    # scalef <- 1
    # fit<-lm(percent_celltype ~ bs(nearest_spline_x_vals,knots = c(-0.11)),data = pc_df )
    # age.grid<-seq(from=-2000*scalef, to = -300*scalef)/10000
    # splinefit <- data.frame(age.grid,predict(fit,newdata = list(nearest_spline_x_vals=age.grid)))
    # colnames(splinefit) <- c("nearest_spline_x_vals","percent_celltype")

    p <- ggplot(pc_df, aes(x=nearest_spline_x_vals, y=percent_celltype, color=DiseaseState)) +
    geom_point(colour="white", shape=21, size = 5.5, aes(fill = factor(DiseaseState)), stroke = 1) + # geom_point(size=4) + 
    scale_fill_manual(values=c("#D51F26", "#89288F", "#208A42")) + theme_ArchR()
    #p <- p+ geom_line(data=splinefit, colour="#CC0000") 
    plotPDF(p, name = paste0("percent_", celltype, "_vs_pseudotime_no_spline_2.pdf"), width = 5, height = 5, ArchRProj = proj_stromal, addDOC = FALSE)   
  }
  library(dplyr)

  for (celltype in colnames(cM)){
    f_celltype <- cM[,celltype]
    percent_celltype <- c()
    for (sample in pc_df$SimplifiedSampleName){
      if (sample %in% names(f_celltype)){
        percent_celltype <- c(percent_celltype, f_celltype[sample])
      } else {
        percent_celltype <- c(percent_celltype, NA)
      }
    }
    pc_df$percent_celltype <- percent_celltype
    pc_df_summary <- pc_df %>% # the names of the new data frame and the data frame to be summarised
      group_by(DiseaseState) %>%   # the grouping variable
      summarise(mean_PL = mean(percent_celltype),  # calculates the mean of each group
                sd_PL = sd(percent_celltype), # calculates the standard deviation of each group
                n_PL = n(),  # calculates the sample size per group
                SE_PL = sd(percent_celltype)/sqrt(n())) # calculates the standard error of each group

    p <- ggplot(pc_df_summary, aes(DiseaseState, mean_PL)) + 
                   geom_col() +  
                   geom_errorbar(aes(ymin = mean_PL - SE_PL, ymax = mean_PL + SE_PL), width=0.2) +
                   theme_ArchR()
    plotPDF(p, name = paste0("percent_", celltype, "_barchart.pdf"), width = 5, height = 5, ArchRProj = proj_stromal, addDOC = FALSE)   
  }

}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 8) Plot tracks
if (8 %in% execute_steps){
  proj_stromal <- addCellColData(proj_stromal, data = paste(getCellColData(proj_stromal, "SimplifiedSampleName")$SimplifiedSampleName, getCellColData(proj_stromal, "CellType")$CellType, sep = "-"), name = "CellTypeSample", cells = rownames(getCellColData(proj_stromal)), force = TRUE)
  pc_df <- readRDS(paste0("./differential_PCA_pseudotime/pc_df_Stem_vs_", "normal",".rds"))

  samples_included <- c("NormalColon-Stem", rev(rownames(pc_df)[rownames(pc_df) %in% paste0(unique(getCellColData(proj_stromal, "SimplifiedSampleName")$SimplifiedSampleName), "-Stem")]))#c("NormalColon-Stem", "A001_Unaffected_Descending-Stem", "A002_Unaffected_Ascending-Stem", "A015_Unaffected_Descending-Stem", "A002_Unaffected_Descending-Stem", rownames(pc_df))
  CellTypeSample <- getCellColData(proj_stromal, "CellTypeSample")$CellTypeSample
  samples_included_new <- c()
  # add a number to the beginning of the name so they plot in the right order
  for (i in 1:length(samples_included)){
    CellTypeSample[CellTypeSample == samples_included[i]] <- paste(i,samples_included[i], sep = "-")
    samples_included_new <- c(samples_included_new, paste(i,samples_included[i], sep = "-"))
  }
  proj_stromal <- addCellColData(proj_stromal, data = CellTypeSample, name = "CellTypeSampleNew", cells = rownames(getCellColData(proj_stromal)), force = TRUE)

  idxSample <- BiocGenerics::which(getCellColData(proj_stromal, "CellTypeSampleNew") %in% samples_included_new)
  cellsSample <- proj_stromal$cellNames[idxSample[["CellTypeSampleNew"]]]
  projTmp <- proj_stromal[cellsSample,]
  markerGenesHyper <- c("ITGA4","NR5A2", "BMP3", "CIDEB", "GRASP")

  plotTracks <- plotBrowserTrack(ArchRProj = projTmp, geneSymbol = markerGenesHyper, groupBy = "CellTypeSampleNew")
  plotPDF(plotList = plotTracks, name = "Plot-Tracks-Hypermethylated", width = 6, height = 8, ArchRProj = projTmp, addDOC = FALSE)

  markerGenesHyper <- c("GRASP")

  plotTracks <- plotBrowserTrack(ArchRProj = projTmp, geneSymbol = markerGenesHyper, groupBy = "CellTypeSampleNew", upstream = 15000, downstream = 15000)
  plotPDF(plotList = plotTracks, name = "Plot-Tracks-GRASP", width = 6, height = 8, ArchRProj = projTmp, addDOC = FALSE)
}


