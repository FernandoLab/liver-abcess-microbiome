# liver-abcess-microbiome

Data analysis methods for filtration and analysis of samples taken from cattle liver abscesses. Data are 16S rDNA amplicon libraries (V3-4 Hypervariable region). This repository contains the R commands run on the ASV count table produced from DADA2, sample metadata, phylogenetic tree of ASV relatedness, and taxonomic classification of ASV sequences. Analysis should be performed in R version 4.0.1 (as of 2020-06-06) in the following order:
  
  1.) la_microbiome_analysis.R
  
  2.) network_family.R
  
  3.) network_genus.R 

If the user would like to start from the raw fastq files the NCBI bioproject information and the commands to generate ASV count table, taxonomic classification, and phylogenetic tree are also included for transparency:

Raw fastq data can be accessed at: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA657040

The commands that were used for fastq processing are also included the LA_microbiome_read_processing.Rmd markdown file. The method for phylogenetic tree generation assumes that Mothur software is previously installed (install - https://mothur.org/wiki/installation/).

citations:
R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
  
Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41
