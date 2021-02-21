# scCRC_continuum

Repository to host code produced for "Single-cell analyses reveal a continuum of cell state and composition changes in the malignant transformation of polyps to colorectal cancer" by Becker*, Nevins*, et al 2021.

**Files:**  
A metadata file is included here that contains primarily the sample specific information included in the supplemental tables as well as some grouping used for calling differential genes/peaks.

The RNA scripts are set up to take the seurat objects for each compartment as input. Seurat objects can be found here: ADD LINK TO SEURAT OBJECTS. 

The ArchR projects are much larger (~100 GB) so are harder to share, but files containing the cells in each compartment are included, which allows you to quickly subset the fragments files to create the ArchR projects.

The raw data for unaffacted, polyp, and CRC samples will be hosted on the the HTAN data portal (https://htan-portal-nextjs.now.sh/) under the PRE-CANCER ATLAS: FAMILIAL ADENOMATOUS POLYPOSIS project. The raw data for normal colon samples will be hosted on the HuBMAP data portal under the Stanford TMC (https://portal.hubmapconsortium.org/search?group_name[0]=Stanford%20TMC&entity_type[0]=Dataset). 

**Scripts:**  
snRNA_epithelial_analysis.R  
-LSI dimensionality reduction and clustering of normal epithelial cells  
-projection of unaffected, polyp, and CRC epithelial cells into normal manifold  
-computation of differentials and malignancy continuum  

snRNA_stromal_analysis.R  
-standard Seurat 3 (Stuart et al. 2019) workflow for dimensionality reduction and clustering of stomal cells from all samples
