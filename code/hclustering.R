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

# Normalization, depth cutoff based on rarefaction 
rare6k.physeq <- rarefy_even_depth(inc.physeq, sample.size = 6000,
                                   rngseed = 15879966) %>%
  filter_taxa(function(x) sum(x) >= 2, T) 

compost.physeq <- subset_samples(rare6k.physeq, treatment %in% c("Compost")) %>%
  filter_taxa(function(x) sum(x) >= 2, T) 
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

tiff("Figures/hclust_compost.tif",height=5,width=6,units='in',res=300)
plot(as.phylo(compost.clustering), type = "unrooted", use.edge.length = TRUE, col = "gray80")
dev.off()

alfalfa.physeq <- subset_samples(rare6k.physeq, treatment %in% c("Alfalfa")) %>%
  filter_taxa(function(x) sum(x) >= 2, T)
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

tiff("Figures/hclust_alfalfa.tif",height=5,width=6,units='in',res=300)
plot(as.phylo(alfalfa.clustering), type = "unrooted", use.edge.length = TRUE, col = "gray80")
dev.off()

reference.physeq <- subset_samples(rare6k.physeq, treatment %in% c("Reference")) %>%
  filter_taxa(function(x) sum(x) >= 2, T) 
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

tiff("Figures/hclust_reference.tif",height=5,width=6,units='in',res=300)
plot(as.phylo(reference.clustering), type = "unrooted", use.edge.length = TRUE, col = "gray80")
dev.off()

mix.physeq <- subset_samples(rare6k.physeq, treatment %in% c("Mix")) %>%
  filter_taxa(function(x) sum(x) >= 2, T)
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

tiff("Figures/hclust_mix.tif",height=5,width=6,units='in',res=300)
plot(as.phylo(mix.clustering), type = "unrooted", use.edge.length = TRUE, col = "gray80")
dev.off()

rare6k.physeq.data <- data.frame(sample_data(rare6k.physeq))
rare6k.physeq.data$response.group[rare6k.physeq.data$day == "0"] <- "baseline" 
rare6k.physeq.data$response.group[rare6k.physeq.data$day %in% c("7", "14", "21")] <- "early" 
rare6k.physeq.data$response.group[rare6k.physeq.data$day %in% c("35", "49", "97")] <- "late" 
sample_data(rare6k.physeq) <- rare6k.physeq.data
saveRDS(rare6k.physeq, "data/IncPhyseqRareClusteredTree")
