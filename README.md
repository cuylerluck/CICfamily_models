# CIC-family modelling 2025

The code in this repository accompanies the manuscript "First Generation Tools for the Modeling of Capicua - (CIC) Family Fusion Oncoprotein-Driven Cancers" by Luck et al (2025). 
These files are not plug-and-play but serve to provide transparency for how we analyzed the data. They may need to be modified for use in different environments or with different data.
With questions please contact the lead author at cuyler.luck@ucsf.edu and the corresponding author at ross.okimoto@ucsf.edu.
RNA-seq raw and processed (counts, log2cpm) data files are available at GEO, series [GSE295624](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295624), [GSE295623](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295623) and [GSE295625](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295625) for the April 2023, May 2024, and Feburary 2025 data, respectively.

There are two main groups of analyses that the code is divided into, with the purpose of each file described below:

## RNA-seq processing and analysis

***documentation_ATXN1_DUX4_Feb2025.txt*** : a text file describing processing of Feb2025 FASTQ files with STAR (and other tools) as performed locally on a MacBook Pro and on Wynton, UCSF's high-performance computing cluster.

***documentation_CICNUTM1_delE_May2024.txt*** : same as above, but for May2024 FASTQ files.

***documentation_CICNUTM1_April2023.txt*** : same as above, but for April2023 FASTQ files.

***293T_ATXN1_DUX4_Feb2025.R*** : R script detailing differential expression analysis of Feb2025 RNA-seq data and used for plotting figures.

***293T_delE_DE_May2024.R*** : R script detailing differential expression analysis of May2024 RNA-seq data and used for plotting figures.

***293T_CICNUTM1_DE_April2023.R*** : R script detailing differential expression analysis of April2023 RNA-seq data and used for plotting figures.


## Schematic generation in Fig 1

***fig1_plasmid_schema.R*** : an R script to generate the pieces of the schematic shown in Figure 1.
