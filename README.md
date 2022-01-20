# scCRC_continuum

Repository to host code produced for "Single-cell analyses reveal a continuum of cell state and composition changes in the malignant transformation of polyps to colorectal cancer" by Becker*, Nevins*, et al 2021. The peer-revieved article is available at (ADD LINK WHEN PUBLISHED) and the preprint is available at https://www.biorxiv.org/content/10.1101/2021.03.24.436532v1.  

![image](https://user-images.githubusercontent.com/15204322/147711771-0a5e3292-095c-443a-a0ef-9d8bf7afb181.png)

**Files/Data:**  
A metadata file is included here that contains primarily the sample specific information included in the supplemental tables as well as some grouping used for calling differential genes/peaks.

The RNA scripts are set up to take the seurat objects for each compartment as input. 

Seurat objects can be found here: ADD LINK TO SEURAT OBJECTS. This includes the following objects, which contain annotated cells from the three compartments with likely doublet and low quality clusters removed.  
epithelial_filtered.rds  
immune_filtered.rds  
stromal_filtered.rds  

The raw data for unaffacted, polyp, and CRC samples will be hosted on the the HTAN data portal (https://htan-portal-nextjs.now.sh/) under the PRE-CANCER ATLAS: FAMILIAL ADENOMATOUS POLYPOSIS project. The raw data for normal colon samples will be hosted on the HuBMAP data portal under the Stanford TMC (https://portal.hubmapconsortium.org/search?group_name[0]=Stanford%20TMC&entity_type[0]=Dataset). 

**Scripts:**

*Preprocessing with cell ranger*

cell_ranger_sbatch_scripts contains run_cell_ranger_ATAC_example.sbatch and run_cell_ranger_RNA_example.sbatch, which are examples of how cell ranger was used to generate expression matricies and fragments files from fastqs. 

*snRNA Analysis Scripts*

scRNA_initial_clustering.R  
-make seurat object will all cells  
-run doublet finder on individual samples  
-divide into three groups (immune, epithelial, and stromal) for downstream analysis  

snRNA_epithelial_analysis.R  
-LSI dimensionality reduction and clustering of normal epithelial cells  
-projection of unaffected, polyp, and CRC epithelial cells into normal manifold  
-computation of differentials and malignancy continuum  

snRNA_stromal_analysis.R  
-standard Seurat 3 (Stuart et al. 2019) workflow for dimensionality reduction and clustering of stomal cells from all samples

snRNA_immune_analysis.R  
-standard Seurat 3 (Stuart et al. 2019) workflow for dimensionality reduction and clustering of stomal cells from all samples

*scATAC Analysis Scripts*

scATAC_initial_clustering.R  
-initial clustering of all scATAC cells with ArchR (Granja et al. 2021)
-note that some aspects of the analysis are non-determenistic (https://www.archrproject.com/bookdown/iterative-latent-semantic-indexing-lsi.html). We have provided lists of cells included in each category to help make steps of the analysis reproducible.

scATAC_subset_immune.R
-clustering and annotation of scATAC immune cells

scATAC_subset_tcells.R
-clustering and annotation of scATAC t-cells

*Methylation Analysis Scripts*


**Questions/Comments:**

For comments/questions, please either add a GitHub issue or email wbecker@stanford.edu.
