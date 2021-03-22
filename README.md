# liver-abcess-microbiome

Data analysis methods for filtration and analysis of samples taken from cattle liver abscesses. Data are 16S rDNA amplicon libraries (V3-4 Hypervariable region). The "sample_data" directory contains the ASV count table produced from DADA2, sample metadata, phylogenetic tree of ASV relatedness, and taxonomic classification of ASV sequences. Also provided are the R scripts used for the analysis and should be proccessed in the following order:
  1.) la_microbiome_analysis.R
  2.) network_family.R
  3.) network_genus.R 

All data should be imported and run in R: (at the time of analysis R version 4.0.1 (2020-06-06))
R Core Team (2018). R: A language and environment for statistical
  computing. R Foundation for Statistical Computing, Vienna, Austria.
  URL https://www.R-project.org/.
  
Raw fastq data can be accessed at: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA657040
Methods for fastq processing is also in the "sample_data" directory, and methods for phylogenetic tree generation assume the user has Mothur software previously installed (install - https://mothur.org/wiki/installation/).

