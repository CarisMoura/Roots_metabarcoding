## Root-fungal associations in transformed landscapes: Insights from network and indicator species analysis
# DNA metabarcoding of plant roots

## Description 
This repository contains scripts for the data processing and downstream analysis pipeline for roots metabarcoding sequences using a dataset comprising four land-use types (rainforest, jungle rubber, rubber, and oil palm plantations) sampled in Indonesia in the framework of the EFForTS project. This project incorporates bioinformatic, correlational analysis, network, and indicator species analysis from plants and fungi detected in a mixture of root samples. 

### Bioinformatic tools
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [Usearch 11](https://www.drive5.com/usearch/download.html)
* [Blast+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
### Usage for data analysis and visualization
```
library(phyloseq)
library(edgeR)
library(BiodiversityR)
library(vegan)
library(ggplot2)
library(indicspecies) ## Indicator species analysis
library(rnetcarto) ## Network analysis
```
 
