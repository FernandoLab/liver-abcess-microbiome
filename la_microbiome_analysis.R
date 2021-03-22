
#################################LIVER ABCSESS MICROBIOME ANALYSIS#######################

#Be sure to either set your work directory to where the data is stored, or change your
#file path to the lines where the images for the data generated.


#Load libraries and install anythin that you don't have:
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(QsRutils)
library(ggpubr)
library(knitr)
library(kableExtra)
library(ampvis2)
library(extrafont)
library(DESeq2)
library(Distanced)
library(RColorBrewer)
library(ampvis2)

##############################FUNCTIONS USED IN THE ANALYSIS########################
# A custom geometric mean function, with zero/NA tolerance.
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#A function to convert phyloseq object to ampvis object:
phyloseq_to_ampvis2 <- function(physeq) {
  #check object for class
  if(!any(class(physeq) %in% "phyloseq"))
    stop("physeq object must be of class \"phyloseq\"", call. = FALSE)
  
  #ampvis2 requires taxonomy and abundance table, phyloseq checks for the latter
  if(is.null(physeq@tax_table))
    stop("No taxonomy found in the phyloseq object and is required for ampvis2", call. = FALSE)
  
  #OTUs must be in rows, not columns
  if(phyloseq::taxa_are_rows(physeq))
    abund <- as.data.frame(phyloseq::otu_table(physeq)@.Data)
  else
    abund <- as.data.frame(t(phyloseq::otu_table(physeq)@.Data))
  
  #tax_table is assumed to have OTUs in rows too
  tax <- phyloseq::tax_table(physeq)@.Data
  
  #merge by rownames (OTUs)
  otutable <- merge(
    abund,
    tax,
    by = 0,
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
  )
  colnames(otutable)[1] <- "OTU"
  
  #extract sample_data (metadata)
  if(!is.null(physeq@sam_data)) {
    metadata <- data.frame(
      phyloseq::sample_data(physeq),
      row.names = phyloseq::sample_names(physeq), 
      stringsAsFactors = FALSE, 
      check.names = FALSE
    )
    
    #check if any columns match exactly with rownames
    #if none matched assume row names are sample identifiers
    samplesCol <- unlist(lapply(metadata, function(x) {
      identical(x, rownames(metadata))}))
    
    if(any(samplesCol)) {
      #error if a column matched and it's not the first
      if(!samplesCol[[1]])
        stop("Sample ID's must be in the first column in the sample metadata, please reorder", call. = FALSE)
    } else {
      #assume rownames are sample identifiers, merge at the end with name "SampleID"
      if(any(colnames(metadata) %in% "SampleID"))
        stop("A column in the sample metadata is already named \"SampleID\" but does not seem to contain sample ID's", call. = FALSE)
      metadata$SampleID <- rownames(metadata)
      
      #reorder columns so SampleID is the first
      metadata <- metadata[, c(which(colnames(metadata) %in% "SampleID"), 1:(ncol(metadata)-1L)), drop = FALSE]
    }
  } else
    metadata <- NULL
  
  #extract phylogenetic tree, assumed to be of class "phylo"
  if(!is.null(physeq@phy_tree)) {
    tree <- phyloseq::phy_tree(physeq)
  } else
    tree <- NULL
  
  #extract OTU DNA sequences, assumed to be of class "XStringSet"
  if(!is.null(physeq@refseq)) {
    #convert XStringSet to DNAbin using a temporary file (easiest)
    fastaTempFile <- tempfile(pattern = "ampvis2_", fileext = ".fa")
    Biostrings::writeXStringSet(physeq@refseq, filepath = fastaTempFile)
  } else
    fastaTempFile <- NULL
  
  #load as normally with amp_load
  ampvis2::amp_load(
    otutable = otutable,
    metadata = metadata,
    tree = tree,
    fasta = fastaTempFile
  )
}
#################################################################################

n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))



#load in the data with a modified sample data Tylan changed to Tylosin
asv_tab <- readRDS("/User_filepath_to_data/seqtab.nochim.rds")
which(row.names(asv_tab) == "S50")
which(row.names(asv_tab) == "S49")
seqtab.nochim.finaldf <- asv_tab[-c(44,46),]
rownames(seqtab.nochim.finaldf)

tax_tab <- readRDS("/User_filepath_to_data/taxa.rds")

asv_tab <- otu_table(seqtab.nochim.finaldf, taxa_are_rows = F) 
tax_tab <- tax_table(tax_tab)
samdat <- as.data.frame(read.delim("/User_filepath_to_data/metadata/sample_data.txt", header = T, row.names = 1, sep = "\t", strip.white = T, na.strings = "NA", stringsAsFactors = TRUE))
samdat <- sample_data(samdat)
tree <- read_tree(treefile = "/User_filepath_to_data/la_microbiome.phylip.tre")

la_ps <- phyloseq(asv_tab,tax_tab,samdat,tree)
la_ps

#Give ASVs an index so we don't have to use the whole sequence: 
la_ps_dna <- Biostrings::DNAStringSet(taxa_names(la_ps))
names(la_ps_dna) <- taxa_names(la_ps)
la_ps <- merge_phyloseq(la_ps, la_ps_dna)
taxa_names(la_ps) <- paste0("ASV", seq(ntaxa(la_ps)))
la_ps
head(refseq(la_ps))
#head(otu_table(la_ps))

##############################ALPHA DIVERSITY OF DATASET WITH PLOTS########################
#Alpha Diversity on raw dataset:
la_ps_trt_alpha <- plot_richness(laps_900up, "Subgroup", measures = c("Shannon","Simpson"), color = "Subgroup") +
  #facet_wrap(Subgroup~., scale="free_x") +
  ylab("Diversity Index") + 
  geom_boxplot(stat = "boxplot") +
  scale_colour_manual(values=col_vector) + 
  ggtitle("Alpha Diversity Measures by Treatment Groups") +
  theme(strip.text.x = element_text(angle=0),
        axis.text.x = element_text(angle=90),
        text = element_text(family="Times New Roman")) +
  stat_compare_means()
la_ps_trt_alpha + theme(legend.text = ggtext::element_markdown(family = 'Times'))
ggsave("/User_filepath_to_data/final_figures/trt_grps_alpha.tiff", width = 30, height = 18, units = "cm")

##############################FILTRATION BASED ON READ ABUNDANCE########################
#Rarefaction curves for samples less than 5,000 reads:
la_ps_less5000 <- prune_samples(rowSums(otu_table(la_ps)) < 5000, la_ps)
rarecurve(as.matrix(otu_table(la_ps_less5000)), step = 50, xlab = "Read Depth", ylab = "Species", label = TRUE)

#remove sample 6
list <- as(rownames(otu_table(la_ps)), "character")
list <- list[-45]                     # Remove list element with !=
list
#make a phyloseq object without sample 6:
la_ps <- prune_samples(samples = list, la_ps)
la_ps
#Make an ampvis object:
la_amp <- phyloseq_to_ampvis2(la_ps)
la_amp

#look at sample read depths
sort(sample_sums(la_ps), decreasing = T)
#there are 5 samples less than 1000 reads S27 is close wth 998 reads,
#but S9 S5 S4, S28 are all less #than 300 reads with S28 only having 85 reads.

# create a phyloseq object w/samples 900 reads or more:
laps_900up <- prune_samples(sample_sums(la_ps)>=900, la_ps)
sort(sample_sums(laps_900up), decreasing = T)
laps_900up <- prune_taxa(taxa_sums(laps_900up)>=1,laps_900up)
#create ampvis object:
laamp_900up<- phyloseq_to_ampvis2(laps_900up)
laamp_900up

##############################PREVALENCE FILTRATION AT 5%########################

#get rid of ASVs that are not classified at the phylum level:
la_trim_ps <- subset_taxa(laps_900up, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
#check if it worked
table(tax_table(la_trim_ps)[, "Phylum"], exclude = NULL)

#creatE a data frame for prevalence of Phyla:
theme_update(plot.title = element_text(hjust = 0.5))
# Compute prevalence of each feature, store as data.frame
prevdf_la = apply(X = otu_table(la_trim_ps),
                  MARGIN = ifelse(taxa_are_rows(la_trim_ps), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_la = data.frame(Prevalence = prevdf_la,
                       TotalAbundance = taxa_sums(la_trim_ps),
                       tax_table(la_trim_ps))
prevdf_la = subset(prevdf_la, Phylum %in% get_taxa_unique(la_trim_ps, "Phylum"))
#plot prevalence of phyla:
ggplot(prevdf_la, aes(TotalAbundance, Prevalence / nsamples(la_trim_ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")+ scale_colour_manual(values=col_vector) + ggtitle("Liver Abscess Bacterial Phyla Prevalence")
#create a prevalence threshold object of 5%
prevalenceThreshold = 0.05 * nsamples(la_trim_ps)
prevalenceThreshold
#filter taxa that are not present in at least 5% of the samples:
keepTaxa = rownames(prevdf_la)[(prevdf_la$Prevalence >= prevalenceThreshold)]
la_ps2 = prune_taxa(keepTaxa, la_trim_ps)
#summarize the number of reads per sample before prevalence filtration and after: 
sample_sums <- sample_sums(la_trim_ps)
filt_sample_sums <- sample_sums(la_ps2)
#check the amount of reads lost from filtration:
filt_read_loss <- cbind(sample_sums, filt_sample_sums)
colSums(filt_read_loss)

#Create an ampvis object of filtered data:
la_amp2 <- phyloseq_to_ampvis2(la_ps2)
la_amp2

##############################PRINCIPLE COORDINATE ANALYSIS AND PLOT########################
#log transform data to normalize (you can also look at relative abundance)
la_ps2.log <- transform_sample_counts(la_ps2, function(x) log(1+x))
#transform data to proportions as well:
#la_ps2.prop <- transform_sample_counts(la_ps2, function(x){x/sum(x)})
#perform ordinations based on phylogenetics or not, and weighted for abundance or not (unifrac weighted unifrac, bray)
la.wuf.log <- ordinate(la_ps2.log, method = "MDS", distance = "wunifrac")
la.uf.log <- ordinate(la_ps2.log, method = "MDS", distance = "unifrac")
la.bray.log <- ordinate(la_ps2.log, method = "MDS", distance = "bray")
#if you would prefer you can just use proportional transformation:
#la.wuf.prop <- ordinate(la_ps2.prop, method = "MDS", distance = "wunifrac")
#la.uf.prop <- ordinate(la_ps2.prop, method = "MDS", distance = "unifrac")
#la.bray.prop <- ordinate(la_ps2.prop, method = "MDS", distance = "bray")

#Adjust titles in pcoa plots to center:
theme(plot.title = element_text(hjust = 0.5))
#Plot principle coordinate analysis of weighted unifrac distance ordinations:
evals <- la.wuf.prop$values$Eigenvalues
plot_ordination(la_ps2.prop, la.wuf.prop, color = "Subgroup") +
  labs(col = "Subgroup") +  
  coord_fixed(sqrt(evals[2] / evals[1])) +
  scale_colour_manual(values=col_vector) +
  ggtitle("Principle Coordinate Analysis of Weighted Unifrac Distances")+
  theme(text = element_text(family="Times New Roman"))
ggsave("/User_filepath_to_data/final_figures/la_pcoa_vegan_wuf_log.tiff")

##############################PERMANOVA ANALYSIS########################
#calculate distances
la_log_wunifrac = phyloseq::distance(la_ps2.log, method = "wunifrac", type = "samples")
#la_log_unifrac = phyloseq::distance(la_ps2.log, method = "unifrac", type = "samples")
#la_log_bray = phyloseq::distance(la_ps2.log, method = "bray", type = "samples")

#extract metadata as a matrix
la_data <- as(sample_data(la_ps2.log), "matrix")
la_data <- as.data.frame(la_data)
class(la_data)

la_wunifrac_permanova <- adonis(la_log_wunifrac ~ Breed +Tylosin +Feedyard, data = la_data, permutations = 999)
la_wunifrac_permanova_table <- kable(la_wunifrac_permanova$aov.tab, "latex", booktabs = T, digits = 10, align = "c") %>%
  kable_styling(latex_options ="striped")%>%
  kable_styling(latex_options = "scale_down")%>%
  footnote(general = "Permutational multivariate analysis of variance based on weighted unifrac distances between liver abscess microbiome samples (999 permutations). The statistical model takes into acount breed, Tylosin treatment, and the interaction between Breed and Tylosin. Significance is determined at alpha < 0.05.",threeparttable = T, general_title = "Table 1. Weighted Unifrac Permanova ", title_format =c("bold"))
#save_kable(la_wunifrac_permanova_table, file ="/User_filepath_to_data/wunifrac_permanova.pdf")


##############################RELATIVE ABUNDANCE HEATMAP FIGURE########################
#create some titles for your treatments:
treatment_names <- list(
  'CB_No_Tylosin'="Crossbred \n No Tylosin",
  'CB_Tylosin'="Crossbred \n With Tylosin",
  'HF_No_Tylosin'="Holstein \n No Tylosin ",
  'HF_Tylosin'="Holstein \n With Tylosin"
)
treatment_labeller <- function(variable,value){
  return(treatment_names[value])
}


la_ps_ampvis <- phyloseq_to_ampvis2(la_ps2)
la_ps_ampvis
#Make a heatmap of the to 25 most abundant ASVs by their genus classification:
amp_heatmap(la_ps_ampvis,
            group_by = c("Subgroup","SampleID"),
            tax_aggregate = "Genus",
            tax_show = 25,
            color_vector = c("white", "darkred"),
            plot_colorscale = "sqrt",
            plot_values = FALSE) +
  facet_grid(.~Subgroup, scale="free_x", space = "free_x", labeller = treatment_labeller) +
  theme(axis.text.x = element_text(angle = 90, size=8, vjust = 1, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold.italic"),
        legend.position="right",
        strip.text.x = element_text(size = 12, face = "bold"),
        text = element_text(family="Times New Roman"))
ggsave("/User_filepath_to_data/final_figures/la_treatment_by_sample_heatmap.tiff", 
       height = 7, width = 11)

##############################DIFFERENTIAL ABUNDANCE HEATMAPS BY BREED########################
library(phylosmith)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(pheatmap)
library(extrafontdb)
library(extrafont)
library(DESeq2)
s <- set_sample_order(t(la_ps2), sort_on = c('Subgroup'))
head(otu_table(s))
r<-as(sample_data(s),'matrix')

dds <- phyloseq_to_deseq2(s, ~Breed)
geoMeans <- apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds,geoMeans=geoMeans)

dds = DESeq(dds, test = "Wald", fitType = "local")

res = results(dds)
resultsNames(dds)

res = res[order(res$padj, na.last=NA), ]
alpha = 0.05

sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(s)[rownames(sigtab), ], "matrix"))
View(sigtab)
#write.table(sigtab, file = "/User_filepath_to_data/ty_vs_noty_sigtab.tsv", sep = "\t")

df <- as.data.frame(colData(dds)[,c("Tylosin","Breed")])

vsd <- varianceStabilizingTransformation(dds)

#Get 25 top varying genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing=TRUE ),10)
order(rowVars(assay(vsd)))

#sig_tax_vector <- c("ASV19","ASV196","ASV1","ASV381","ASV261","ASV357","ASV251","ASV13","ASV348")
sig_tax_vector <- c("ASV392","ASV7","ASV19","ASV261", "ASV391","ASV44")
#make a subset of the log transformed counts for just the top 25 varying genes
deseq2_diff_Counts<-assay(vsd)[sig_tax_vector,]
#ddif2<- rel[sig_tax_vector,]
#Use pheatmap function to draw a heatmap
#INCLUDE NEXT LINE IF YOU WANT TO SAVE THE FIGURE IN A FILE
#pdf(file="gene.heatmap.pdf")

tax_vector <- sigtab[12]
row.names(tax_vector)<-row.names(sigtab)
tax_vector <- cbind(tax_vector,row.names(tax_vector))
tax_vector$col3 <- paste(tax_vector$`row.names(tax_vector)`,tax_vector$Genus, sep = "-")


rownames(deseq2_diff_Counts) <- tax_vector$col3
rownames(deseq2_diff_Counts)

colnames(df) <- c("Tylosin", "Cattle Type")

la_heatmap <- pheatmap::pheatmap(deseq2_diff_Counts,scale = "row", 
                                 annotation_col = df,
                                 #annotation_row = tax_vector$col3, 
                                 cluster_rows = T, 
                                 cluster_cols = F, 
                                 gaps_col = 20, 
                                 cutree_rows = 2 ,
                                 show_rownames=T,
                                 fontsize_col = 7,
                                 fontfamily='Times New Roman',
                                 fontface="italic",
                                 width = 7, 
                                 height = 5,
                                 clustering_distance_rows = "euclidean",
                                 fontface_row = "italic",
                                 filename = "/User_filepath_to_data/final_figures/diff_breed_heatmap_deseq2.pdf")

class(la_heatmap)
la_heatmap

##############################DIFFERENTIAL ABUNDANCE HEATMAPS BY TYLOSIN TREATMENT########################

s <- set_sample_order(t(la_ps2), sort_on = c('Subgroup'))
head(otu_table(s))
r<-as(sample_data(s),'matrix')

dds <- phyloseq_to_deseq2(s, ~Tylosin)
geoMeans <- apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds,geoMeans=geoMeans)

dds = DESeq(dds, test = "Wald", fitType = "local")

res = results(dds)
resultsNames(dds)

res = res[order(res$padj, na.last=NA), ]
alpha = 0.05

sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(s)[rownames(sigtab), ], "matrix"))
View(sigtab)
#write.table(sigtab, file = "/User_filepath_to_data/ty_vs_noty_sigtab.tsv", sep = "\t")

df <- as.data.frame(colData(dds)[,c("Tylosin","Breed")])

vsd <- varianceStabilizingTransformation(dds)

#Get 25 top varying genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing=TRUE ),10)
order(rowVars(assay(vsd)))

sig_tax_vector <- c("ASV19","ASV196","ASV1","ASV381","ASV261","ASV357","ASV251")
#make a subset of the log transformed counts for just the top 25 varying genes
deseq2_diff_Counts<-assay(vsd)[sig_tax_vector,]
#ddif2<- rel[sig_tax_vector,]
#Use pheatmap function to draw a heatmap
#INCLUDE NEXT LINE IF YOU WANT TO SAVE THE FIGURE IN A FILE
#pdf(file="gene.heatmap.pdf")

tax_vector <- sigtab[12]
row.names(tax_vector)<-row.names(sigtab)
tax_vector <- cbind(tax_vector,row.names(tax_vector))
tax_vector$col3 <- paste(tax_vector$`row.names(tax_vector)`,tax_vector$Genus, sep = "-")


rownames(deseq2_diff_Counts) <- tax_vector$col3
rownames(deseq2_diff_Counts)

colnames(df) <- c("Tylosin", "Cattle Type")

la_heatmap <- pheatmap::pheatmap(deseq2_diff_Counts,scale = "row", 
                                 annotation_col = df,
                                 #annotation_row = tax_vector$col3, 
                                 cluster_rows = T, 
                                 cluster_cols = F, 
                                 gaps_col = 20, 
                                 cutree_rows = 2 ,
                                 show_rownames=T,
                                 fontsize_col = 7,
                                 fontfamily='Times New Roman',
                                 fontface="italic",
                                 width = 7, 
                                 height = 5,
                                 fontface_row = "italic",
                                 filename = "/User_filepath_to_data/final_figures/diff_tylosin_heatmap_deseq2.pdf")

class(la_heatmap)

##############################NETWORKS SAVED IN DIFFERENT R SCRIPTS########################






