# liver-abcess-microbiome

Data analysis methods for filtration and analysis of samples taken from cattle liver abscesses. Data are 16S rDNA amplicon libraries (V3-4 Hypervariable region). This repository contains the R commands run on the ASV count table produced from DADA2, sample metadata, phylogenetic tree of ASV relatedness, and taxonomic classification of ASV sequences. Analysis should be performed in R version 4.0.1 (as of 2020-06-06) in the following order:
  
  1.) la_microbiome_analysis.R
  
  2.) network_family.R
  
  3.) network_genus.R 

  
Raw fastq data can be accessed at: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA657040

Methods for fastq processing is also in the "sample_data" directory, and methods for phylogenetic tree generation assume the user has Mothur software previously installed (install - https://mothur.org/wiki/installation/).

citation: R Core Team (2018). R: A language and environment for statistical
  computing. R Foundation for Statistical Computing, Vienna, Austria.
  URL https://www.R-project.org/.
