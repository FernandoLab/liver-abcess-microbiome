
# Install Required packages
#install.packages("igraph")
#install.packages("qgraph")
#install.packages("vegan")
#install.packages("MCL")
# Install SpiecEasi package
#install.packages("devtools")
#library(devtools)
#install_github("zdk123/SpiecEasi")

#install.packages("igraph", dependencies = TRUE)

setwd("/filepath_to_data/")
set.seed(123456)

library(igraph)
library(qgraph)
library(vegan)
library(MCL)
library(SpiecEasi)



# Load OTU table
#data("amgut1.filt")
# Change row and column names to a more readable format
#colnames(amgut1.filt) <- sapply(1:ncol(amgut1.filt),
# function(x) paste("OTU", x, sep = "_"))
#rownames(amgut1.filt) <- sapply(1:nrow--(amgut1.filt),
# function(x) paste("Sample", x, sep = "_"))


#user.table <- read.table(infile,
#                        header = T,
#                       row.names = 1,
#                      sep = ",")
la_family_ps <- tax_glom(la_ps2, taxrank = 'Family', NArm = T)
asv.table<- as(otu_table(la_family_ps),'matrix')
class(asv.table)

#change the asv consensus number to taxonomy name
colnames(asv.table) <- tax_table(la_family_ps)[,"Family"]

amgut1.filt <- asv.table

otu.relative <- amgut1.filt / rowSums(amgut1.filt)


# Create dissimilarity matrix
#distances <- vegdist(t(otu.relative),
#                     method = "bray")


# Convert distance object to a matrix
#diss.mat <- as.matrix(distances)
#diss.cutoff <- 0.6
diss.adj <- ifelse(diss.mat <= diss.cutoff, 1, 0)
# Construct microbiome network from adjacency matrix
diss.net <- graph.adjacency(diss.adj,
                            mode = "undirected",
                            diag = FALSE)





#correlation based network
#cor.matrix <- cor(otu.relative, method = "pearson")


# Convert correlation matrix to binary adjacency matrix
#cor.cutoff <- 0.3
#cor.adj <- ifelse(abs(cor.matrix) >= cor.cutoff, 1, 0)
# Construct microbiome network from adjacency matrix
#cor.net <- graph.adjacency(cor.adj,
#                           mode = "undirected",
#                           diag = FALSE)


# Execute this command after running Function 1
#cor.net.2 <- build.cor.net(amgut1.filt,
#                           method = 'pearson',
#                           num_perms = 100,
#                           sig_level = 0.01)
##### Function 1: Construct microbiome network using permutation
#build.cor.net <- function(MB, method, num_perms, sig_level) {
#  taxa <- dim(MB)[2]
#  MB.mat <- array(0, dim = c(taxa, taxa, num_perms + 1))
# Perform permutation
#  MBperm <- permatswap(MB, 'quasiswap', times = num_perms)
# Convert to relative abundance
#  MB.relative <- MB / rowSums(MB)
#  MB.mat[,,1] <- as.matrix(cor(MB.relative, method = method))
#  for(p in 2:num_perms) {
#    MBperm.relative <- MBperm$perm[[p-1]] / rowSums(MBperm$perm[[p-1]])
#    MB.mat[, , p] <- as.matrix(cor(MBperm.relative, method = method))
#  }
# Get p-values
#  pvals <- sapply(1:taxa,
#                  function(i) sapply(1:taxa, function(j)
#                    sum(MB.mat[i, j, 1] > MB.mat[i, j, 2:num_perms])))
#  pvals <- pvals / num_perms
# p-value correction
#  pvals_BH <- array(p.adjust(pvals, method = 'BH'),
#                    dim=c(nrow(pvals), ncol(pvals)))
# Build adjacency matrix
#  adj.mat <- ifelse(pvals_BH >= (1 - sig_level), 1, 0)
# Add names to rows & cols
#  rownames(adj.mat) <- colnames(MB)
#  colnames(adj.mat) <- colnames(MB)
# Build and return the network
#  graph <- graph.adjacency(adj.mat, mode = 'undirected', diag = FALSE)
#}


#Graphical Model Networks



# Compute (partial) correlations
#ebic.cor <- cor_auto(amgut1.filt)
# Identify graph with the best EBIC
#ebic.graph <- EBICglasso(ebic.cor, ncol(amgut1.filt), 0.5)
# Build the network
#ebic.qgnet <- qgraph(ebic.graph, DoNotPlot = TRUE)
# Convert to igraph network
#ebic.net <- as.igraph(ebic.qgnet, attributes = TRUE)



# Compute (partial) correlations
#fdr.cor <- cor_auto(amgut1.filt)
# Identify graphical model
#fdr.graph <- FDRnetwork(fdr.cor, cutoff = 0.01, method = "pval")
# Build the network
#fdr.qgnet <- qgraph(fdr.graph, DoNotPlot = TRUE)
# Convert to igraph network
#fdr.net <- as.igraph(fdr.qgnet, attributes = TRUE)


#SparCC and SPIEC-EASI Networks




# SparCC network
sparcc.matrix <- sparcc(amgut1.filt)
sparcc.cutoff <- 0.3
sparcc.adj <- ifelse(abs(sparcc.matrix$Cor) >= sparcc.cutoff, 1, 0)
# Add OTU names to rows and columns
rownames(sparcc.adj) <- colnames(amgut1.filt)
colnames(sparcc.adj) <- colnames(amgut1.filt)
# Build network from adjacency
sparcc.net <- graph.adjacency(sparcc.adj,
                              mode = "undirected",
                              diag = FALSE)



# SPIEC-EASI network
#SpiecEasi.matrix <- spiec.easi(amgut1.filt,
#                               method = 'glasso',
#                               lambda.min.ratio = 1e-2,
#                               nlambda = 20,
#                               icov.select.params = list(rep.num = 50))
# Add OTU names to rows and columns
#rownames(SpiecEasi.matrix$refit$stars) <- colnames(amgut1.filt)
#colnames(SpiecEasi.matrix$refit$stars) <- colnames(amgut1.filt)

#head(SpiecEasi.matrix$refit)
# Build network from adjacency
#SpiecEasi.net <- graph.adjacency(SpiecEasi.matrix$refit$stars,
#                                 mode = "undirected",
#                                 diag = FALSE)



#Hub Detection


# Use sparcc.net for the rest of the method
net <- sparcc.net
# Hub detection
net.cn <- closeness(net)
net.bn <- betweenness(net)
net.pr <- page_rank(net)$vector
net.hs <- hub_score(net)$vector


# Sort the species based on hubbiness score
net.hs.sort <- sort(net.hs, decreasing = TRUE)
# Choose the top 5 keystone species
net.hs.top5 <- head(net.hs.sort, n = 5)



#Cluster Detection


# Get clusters
wt <- walktrap.community(net)
ml <- multilevel.community(net)
# Get membership of walktrap clusters
membership(wt)
# Get clusters using MCL method
adj <- as_adjacency_matrix(net)
mc <- mcl(adj, addLoops = TRUE)


# Compare clusters detected by different methods
compare(membership(wt), membership(ml))
compare(membership(wt), mc$Cluster)
# Create customized membership for comparison
expected.cls <- sample(1:5, vcount(net), replace = T) %>% as_membership
compare(expected.cls, membership(wt))
# Plot clusters as dendrogram
plot_dendrogram(wt)


# Calculate modularity
modularity(net, membership(wt))




#Network Simulation


# Simulate networks
#num.nodes <- 50
#regular.net <- k.regular.game(num.nodes, k = 4)
#random.net <- erdos.renyi.game(num.nodes, p = 0.037)
#smallworld.net <- sample_smallworld(dim = 1, num.nodes, nei = 2, p = 0.2)
#scalefree.net <- barabasi.game(num.nodes)


#Network Features


# Network features
nodes <- V(net)
edges <- V(net)
node.names <- V(net)$name
num.nodes <- vcount(net)
num.edges <- ecount(net)
num.edges
num.nodes
clustering_coeff <- transitivity(net, type = "global")


# Obtain the neighbors of nodes 1 and 25
#asv1_neighbors <- neighbors(net, "ASV1")
#asv114_neighbors <- neighbors(net, "ASV114")
# Find neighbors shared by nodes 1 and 25
#intersection(asv1_neighbors, asv114_neighbors)


# Edges incident to OTU_1
#asv1.edges <- incident(net, "ASV1", mode = "all")
# Edges incident to OTU_1 and OTU_25
#otus.edges <- incident_edges(net, c("OTU_1", "OTU_25"), mode = "all")
# Extracting/printing the incident edges separately
#otus.edges$"OTU_1"
#otus.edges$"OTU_25"


net.knn <- knn(net, vids = V(net))
net.knn$knn

#sub.node1 <- subcomponent(net, v = "OTU_1", mode = "all")


clean.net <- delete.vertices(net, which(degree(net, mode = "all")
                                        == 0))



# Network components
net.comps <- components(net)
# Print components membership
net.comps$membership
# Print components sizes
net.comps$csize
# Print number of components
net.comps$no



# Largest component
#largest.comp <- V(net)[which.max(net.comps$csize) == net.comps$membership]
# Second component
#second.comp <- V(net)[net.comps$membership == 2]
# The component containing OTU_1
#otu1.comp <- V(net)[net.comps$membership ==
#which(names(net.comps$membership) == "OTU_1")]
# Largest component subnetwork
#largest.subnet <- induced_subgraph(net, largest.comp)
# Subnetwork for the component containing OTU_1
#otu1.subnet <- induced_subgraph(net, otu1.comp)


# Degrees
deg <- degree(net, mode = "all")
# Degree distribution
deg.dist <- degree_distribution(net, mode = "all", cumulative = T)
# Plot degree distribution
plot(deg.dist, xlab = "Nodes degree", ylab = "Probability")
lines(deg.dist)
# qgraph method
centralityPlot(net)


# Scalefreeness: Fit a power_law to the network
deg <- degree(net, mode = "in")
pl <- fit_power_law(deg, xmin = 10)
pl$KS.p
# Smallworldness
sw <- smallworldness(net, B = 10)
sw[1]



#intsect.edges <- intersection(net, otu1.subnet)
#union.edges <- union(net, otu1.subnet)

#node.similarity <- similarity(net, vids = V(net), mode = "all", method = "jaccard")


# Find articulation points
AP <- articulation.points(net)
AP




#Network Visualization


# Simple plotting
plot(net)

#pdf("/Volumes/Wes_External/LA_microbiome/sparcc_net_phylum.pdf", width = 15, height = 10 )
plot(wt, net)
#dev.off()


# Customized plotting
plot(net,
     main = "Microbiome Network",
     vertex.color = "white",
     vertex.size = 5,
     vertex.shape = "circle",
     vertex.frame.color = "green",
     Vertex.label.size = 1,
     Vertex.label.color = "black",
     edge.color = "grey",
     layout = layout.fruchterman.reingold)



# Function 2: Plot network with node size scaled to hubbiness
plot.net <- function(net, scores, outfile, title) {
  # Convert node label from names to numerical IDs.
  features <- V(net)$name
  col_ids <- seq(1, length(features))
  V(net)$name <- col_ids
  node.names <- features[V(net)]
  # Nodes' color.
  V(net)$color <- "white"
  # Define output image file.
  outfile <- paste(outfile, "tiff", sep=".")
  # Image properties.
  jpeg(outfile, width = 4800, height = 9200, res = 300, quality = 100)
  par(oma = c(4, 1, 1, 1))
  # Main plot function.
  plot(net, vertex.size = (scores*5)+4, vertex.label.cex = 1)
  title(title, cex.main = 4)
  # Plot legend containing OTU names.
  labels = paste(as.character(V(net)), node.names, sep = ") ")
  legend("bottom", legend = labels, xpd = TRUE, ncol = 5, cex = 1.2)
  dev.off()
}

# Execute this command after running Function 2
#plot.net(net, net.hs, outfile = "network_la_phylum", title = "SparCC Network: Phylum")




# Function 3: Plot network with clusters and node size scaled to hubbiness
plot.net.cls <- function(net, scores, cls, AP, outfile, title) {
  # Get size of clusters to find isolated nodes.
  cls_sizes <- sapply(groups(cls), length)
  # Randomly choosing node colors. Users can provide their own vector of colors.
  colors <- sample(colours(), length(cls))
  # Nodes in clusters will be color coded. Isolated nodes will be white.
  V(net)$color <- sapply(membership(cls),
                         function(x) {ifelse(cls_sizes[x]>1,
                                             colors[x], "white")})
  # Convert node label from names to numerical IDs.
  node.names <- V(net)$name
  col_ids <- seq(1, length(node.names))
  V(net)$name <- col_ids
  # To draw a halo around articulation points.
  #AP <- lapply(names(AP), function(x) x)
  #marks <- lapply(1:length(AP), function(x) which(node.names == AP[[x]]))
  # Define output image file.
  outfile <- paste(outfile, "jpg", sep=".")
  # Image properties.
  jpeg(outfile, width = 7200, height = 8700, res = 300, quality = 100)
  par(oma = c(4, 1, 1, 1))
  # Customized layout to avoid nodes overlapping.
  e <- get.edgelist(net)
  class(e) <- "numeric"
  l <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(net),
                                         area=8*(vcount(net)^2),
                                         repulse.rad=(vcount(net)^3.1))
  # Main plot function.
  plot(net, vertex.size = (scores*5)+4, vertex.label.cex=0.9,
       vertex.label.color = "black",
       mark.border="black",
       #mark.groups = marks,
       mark.col = "white",
       mark.expand = 10,
       mark.shape = 1,
       layout=l)
  title(title, cex.main=4)
  # Plot legend containing OTU names.
  labels = paste(as.character(V(net)), node.names, sep = ") ")
  legend("bottom", legend = labels, xpd = TRUE, ncol = 5, cex = 1.2)
  dev.off()
}
# Execute this command after running Function 3
#plot.net.cls(net, net.hs, wt, AP, outfile = "sparcc_la_network2_phylum", title = "SparCC Network: Phylum")


# Convert node label from names to numerical IDs.
features <- V(net)$name
col_ids <- seq(1, length(features))
V(net)$name <- col_ids
node.names <- features[V(net)]
wt <- walktrap.community(net)

pdf("/Users/westom/Desktop/sparcc_family_net2.pdf", width = 18, height = 25)
par(family = 'serif')
plot(wt,net,label.family='Times New Roman',vertex.label.cex=2.5)
title("SparCC Network: Family", cex.main=4)
labels = paste(as.character(V(net)), node.names, sep = ") ")
#legend("bottom", legend = labels, xpd = TRUE, ncol = 3, cex = 2, bty = 'n')
#title("SparCC Network: Family", cex.main=4)
dev.off()





