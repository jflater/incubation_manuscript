#Phylosmith
library(devtools)
install_github('schuyler-smith/phylosmith')
library(phylosmith)
library(phyloseq)
library(tidyverse)

# Load phyloseq object
inc.raw.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

# List of OTUs from each alien group
alf.aliens <- as.character(readRDS("data/alf.aliens.rds"))
comp.aliens <- readRDS("data/comp.aliens.rds")
mix.aliens <- readRDS("data/mix.aliens.rds")

# Remove day zero incubated samples, will also remove ingredient samples
physeq <- subset_samples(inc.raw.physeq, day %in% c("7",
                                                        "14",
                                                        "21",
                                                        "35",
                                                        "49",
                                                        "97"))

# prune out zero samples, though note here if work to generate the alf.aliens was on mintax 1, 5, or?
alf.a <- prune_taxa(alf.aliens, physeq) %>%
  filter_taxa(function(x) sum(x) > 1, T)

# This will make a list of samples in the correct day order
sample.order <- as.data.frame(sample_data(alf.a)) %>%
  arrange(day, replication) %>%
  select(i_id) %>%
  remove_rownames() 

abundance_heatmap_ggplot(alf.a, classification = 'Genus', treatment = c('treatment'), transformation = 'log')

plot_heatmap(alf.a, sample.label = "day", taxa.order = "Phylum", taxa.label = "Genus", trans = 'log',  
                              sample.order = as.character(sample.order$i_id), 
                              low = "yellow", high = "red", na.value = "grey") +
  facet_grid(. ~ treatment)

