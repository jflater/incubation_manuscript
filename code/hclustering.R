library(phyloseq)
library(vegan)
library(tidyverse)
library(ape)
library(ggtree)
library(ggplot2)

tree <- read.tree("data/tree.nwk")
tree
inc_phy <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

# Add the tree file to the phyloseq object
inc_phy <- merge_phyloseq(inc_phy, tree)
inc_phy
inc.physeq <- subset_samples(inc_phy, day %in% c("7",
                                                        "14",
                                                        "21",
                                                        "35",
                                                        "49",
                                                        "97"))

#Rename treatments to more informative titles
data <- data.frame(sample_data(inc.physeq)) %>%
  mutate(treatment = recode(treatment,
                            'Control' = 'Reference',
                            'CompAlfa' = 'Mix')) %>%
  mutate(C_N = C_flash / N_flash, Inorganic_N = NH3 + NO3) %>%
  mutate(TreatmentAndDay = paste(treatment, day))
rownames(data) <- data$i_id
sample_data(inc.physeq) <- data
sample_data(inc.physeq)$day <- as.factor(sample_data(inc.physeq)$day)
sample_data(inc.physeq)

no.unclass <- subset_taxa(inc.physeq, !Phylum=="Bacteria_unclassified")
no.unclass <- subset_taxa(no.unclass, !Genus=="Gp6_unclassified")
rare6k.physeq <- rarefy_even_depth(physeq = no.unclass, sample.size = 6000, rngseed = 3242343, verbose = F)
rare6k.physeq

####Compost####
compost.physeq <- subset_samples(rare6k.physeq, treatment %in% c("Compost")) %>%
  filter_taxa(function(x) sum(x) >= 3, T) 
# get day and i_id for these samples, we will use this list to order the days on the tree
compost.rownames <- data.frame(sample_data(compost.physeq)) %>%
  select(day, i_id)
# List of each samples day label to use as rownames for the clustering plot
rownamesforcompost <- as.character(compost.rownames$day)
# Distance matrix, setting binary to T compustes a presence absence for the OTU table
compost.physeq.dist <- vegdist(t(data.frame(otu_table(compost.physeq))), method = "bray", binary = T)
# Hierachical clustering to check if days are grouped
compost.clustering <- hclust(compost.physeq.dist, method = "ward.D2")
compost.clustering$labels <- rownamesforcompost

treec <- as.phylo(compost.clustering)
colors <- c("7"="blue", "14"="blue", "21"="blue", "35"="red", "49"="red", "97"="red")
png("Figures/hclust_compost_test.png",height=5,width=6,units='in',res=300)
plot(treec, type = "unrooted", use.edge.length = TRUE, col = "gray80", cex = .4, font = 1, tip.color=colors[substr(treec$tip.label, 1, 6)])
dev.off()

####Alfalfa####
alfalfa.physeq <- subset_samples(rare6k.physeq, treatment %in% c("Alfalfa")) %>%
  filter_taxa(function(x) sum(x) >= 3, T)
# get day and i_id for these samples, we will use this list to order the days on the tree
alfalfa.rownames <- data.frame(sample_data(alfalfa.physeq)) %>%
  select(day, i_id)
# List of each samples day label to use as rownames for the clustering plot
rownamesforalfalfa <- as.character(alfalfa.rownames$day)
# Distance matrix
alfalfa.physeq.dist <- vegdist(t(data.frame(otu_table(alfalfa.physeq))), method = "bray", binary = T)
# Hierachical clustering to check if days are grouped
alfalfa.clustering <- hclust(alfalfa.physeq.dist, method = "ward.D2")
alfalfa.clustering$labels <- rownamesforalfalfa

treea <- as.phylo(alfalfa.clustering)
colors <- c("7"="blue", "14"="blue", "21"="blue", "35"="red", "49"="red", "97"="red")

png("Figures/hclust_alfalfa.png",height=5,width=6,units='in',res=300)
plot(treea, type = "unrooted", use.edge.length = TRUE, col = "gray80", cex = .4, font = 1, tip.color=colors[substr(treea$tip.label, 1, 6)])
dev.off()

####Reference####
reference.physeq <- subset_samples(rare6k.physeq, treatment %in% c("Reference")) %>%
  filter_taxa(function(x) sum(x) >= 3, T) 
# get day and i_id for these samples, we will use this list to order the days on the tree
reference.rownames <- data.frame(sample_data(reference.physeq)) %>%
  select(day, i_id)
# List of each samples day label to use as rownames for the clustering plot
rownamesforreference <- as.character(reference.rownames$day)
# Distance matrix
reference.physeq.dist <- vegdist(t(data.frame(otu_table(reference.physeq))), method = "bray", binary = T)
# Hierachical clustering to check if days are grouped
reference.clustering <- hclust(reference.physeq.dist, method = "ward.D2")
reference.clustering$labels <- rownamesforreference

treer <- as.phylo(reference.clustering)
colors <- c("7"="blue", "14"="blue", "21"="blue", "35"="red", "49"="red", "97"="red")

png("Figures/hclust_reference.png",height=5,width=6,units='in',res=300)
plot(treer, type = "unrooted", use.edge.length = TRUE, col = "gray80", cex = .4, font = 1, tip.color=colors[substr(treer$tip.label, 1, 6)])
dev.off()

####Mix####
mix.physeq <- subset_samples(rare6k.physeq, treatment %in% c("Mix")) %>%
  filter_taxa(function(x) sum(x) >= 3, T)
# get day and i_id for these samples, we will use this list to order the days on the tree
mix.rownames <- data.frame(sample_data(mix.physeq)) %>%
  select(day, i_id)
# List of each samples day label to use as rownames for the clustering plot
rownamesformix <- as.character(mix.rownames$day)
# Distance matrix
mix.physeq.dist <- vegdist(t(data.frame(otu_table(mix.physeq))), method = "bray", binary = T)
# Hierachical clustering to check if days are grouped
mix.clustering <- hclust(mix.physeq.dist, method = "ward.D2")
mix.clustering$labels <- rownamesformix

treem <- as.phylo(mix.clustering)
colors <- c("7"="blue", "14"="blue", "21"="blue", "35"="red", "49"="red", "97"="red")

png("Figures/hclust_mix.png",height=5,width=6,units='in',res=300)
plot(treem, type = "unrooted", use.edge.length = TRUE, col = "gray80", cex = .4, font = 1, tip.color=colors[substr(treem$tip.label, 1, 6)])
dev.off()



