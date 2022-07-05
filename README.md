# Single-cell-multi-omics-integration

Single-cell multi-omics data analysis pipeline for integrating subclone structures from paired genome and transcriptome data in primary gastric and metastatic colon cancers. There are three major parts included in the whole pipeline:
* (1) scDNA-seq analysis workflow:
* Cell quality control (QC)
* Identified cellular components 
* Constructed subclones 
* (2) scRNA-seq analysis: 
*  Assigned cell cycle 
*  Annotated cell types 
*  Extracted the epithelial cell
* (3) Integration analysis 
* Assigned epithelial single cells (G0G1 phase) from scRNA-seq into subclones of scDNA-seq
* Evaluated the subclone assignment by inferCNV package 
* Discovered phenotype biology of subclones  

## Data
The scDNA-seq and scRNA-seq datasets generated for this study are available in NCBI's dbGAP repositories, accession numbers phs001711 and phs001818. 

## Reference
[Single cell multi-omic mapping of subclonal architecture and pathway phenotype in primary gastric and metastatic colon cancers, bioRxiv](https://www.biorxiv.org/content/10.1101/2022.07.03.498616v1)

## Contact 
### xiangqi@stanford.edu 
