if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")



library(tidyverse)
library(phyloseq)
library(ggtree)
library(ggplot2)
library(scales)

tree <- read.tree("data/tree.nwk")
tree

inc_phy <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

# Add the tree file to the phyloseq object
inc_phy <- merge_phyloseq(inc_phy, tree)
inc_phy

#subset
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
head(sample_data(inc.physeq))

# Rarefy
rare6k.physeq <- rarefy_even_depth(inc.physeq, sample.size = 6000,
                                   rngseed = 15879966) %>%
  filter_taxa(function(x) sum(x) >= 1, T) 

firmc <- subset_taxa(rare6k.physeq, Phylum=="Firmicutes") %>%
  filter_taxa(function(x) sum(x) >= 1, T) %>%
  tax_glom("Genus")

p <- ggtree(firmc, ladderize = FALSE) + geom_text2(aes(subset=!isTip, label=label), size=4) +
  geom_tiplab(aes(label=Genus)) +
  theme(legend.position="right") + ggtitle("Firmicutes")
print(p)
