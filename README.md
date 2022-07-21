# Single-cell-multi-omics-integration

Single-cell multi-omics data analysis pipeline for defining the subclone architecture and clonal phenotype from the primary gastric and metastatic colon cancers. We conducted scDNA-seq and scRNA-seq assays on thousands of single cells derived from the same specimen. We constructed the subclone architecture by scDNA-seq data and characterized the distinct large CNVs among subclones. Then we analyzed the paired scRNA-seq gene expression profile to determine cell types and transcriptome features. Utilizing both results, we assigned the scRNA-seq gene expression to subclones derived in scDNA-seq based on gene dosage effect. Finally, we identified the differential phenotype pathways among subclones by GSVA. There are three major parts included in the whole pipeline:
## scDNA-seq analysis:
 1. Cell quality control (QC)
 2. Identified cellular components 
 3. Constructed subclones 
 
 ```
 Rscript Preprocessing_data_scDNA.R
 Rscript run_pipeline_scDNA.R
 Rscript Plot_umap_scDNA.R
 Rscript Plot_evolution_scDNA.R
 Rscript Plot_subclones_heatmap_scDNA.R
 ```
## scRNA-seq analysis: 
 1. Assigned cell cycle 
 2. Annotated cell types 
 3. Extracted the epithelial cell
## Integration analysis 
 1. Assigned epithelial single cells (G0G1 phase) from scRNA-seq into subclones of scDNA-seq
 2. Evaluated the subclone assignment by inferCNV package 
 3. Discovered phenotype biology of subclones  

## Data
The scDNA-seq and scRNA-seq datasets generated for this study are available in NCBI's dbGAP repositories, accession numbers phs001711 and phs001818. 

## Reference
[Single cell multi-omic mapping of subclonal architecture and pathway phenotype in primary gastric and metastatic colon cancers, bioRxiv](https://www.biorxiv.org/content/10.1101/2022.07.03.498616v1)

## Contact 
### xiangqi@stanford.edu 
