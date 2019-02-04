library(phyloseq)
library(tidyverse)
library(gplots)

# Read in whole phyloseq object
incubation.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

# See what catgegories describe our treatments:
levels(sample_data(incubation.physeq)$treatment)

# Function to get a list of OTUs from a sample type
GetOTUs <- function(physeq, samples) {
    prune_samples(sample_data(physeq)$treatment %in% c(samples), physeq) %>%
    filter_taxa(function(x) sum(x) > 1, T) %>%
    tax_table() %>%
    row.names()
}

# Using function to generate list of OTUs from all starting soil, incubated soils and amendments
# Incubated samples:
Alfalfa.otus <- GetOTUs(incubation.physeq, c("Alfalfa"))
Compost.otus <- GetOTUs(incubation.physeq, c("Compost"))
CompAlfa.otus <- GetOTUs(incubation.physeq, c("CompAlfa"))
Control.otus <- GetOTUs(incubation.physeq, c("Control"))
# Amendment samples:
AlfalfaAmend.otus <- GetOTUs(incubation.physeq, c("AlfalfaAmend"))
CompostAmend.otus <- GetOTUs(incubation.physeq, c("CompostAmend"))
# Soil sample:
AlfalfaSoil.otus <- GetOTUs(incubation.physeq, c("AlfalfaSoil"))


GetAlienHeatMap <- function(physeq, control_otus, alien_otus, recieving_otus,samples){
  otus <- list(alien_otus, control_otus)
  venn <- venn(otus)
  alf.aliens <- attr(venn, "intersections")$A
  aliens <- list(alf.aliens, recieving_otus)
  aliens.venn <- venn(aliens)
  aliens.detected <- attr(aliens.venn,"intersections")$`A:B`
  alfalfa <- prune_samples(sample_data(physeq)$treatment %in% c(samples), physeq) %>%
    filter_taxa(function(x) sum(x) > 1, T) 
  alf.incubated.aliens <- prune_taxa(aliens.detected, alfalfa)
  alf.sample.order <- as.data.frame(sample_data(alf.incubated.aliens)) %>%
    arrange(day, replication) %>%
    select(i_id) %>%
    remove_rownames()   
  alf.alien.heatmap <- plot_heatmap(alf.incubated.aliens, sample.label = "day", taxa.order= "Phylum", taxa.label = "Genus",
                                  sample.order = as.character(alf.sample.order$i_id), 
                                  low = "#66CCFF", high = "#000033", na.value = "white")
  alf.alien.heatmap
}
GetAlienHeatMap(incubation.physeq, Control.otus, AlfalfaAmend.otus, Alfalfa.otus, c("Alfalfa"))
GetAlienHeatMap(incubation.physeq, Control.otus, CompostAmend.otus, Compost.otus, c("Compost"))

otus <- list(CompostAmend.otus, Control.otus)
otus
venn <- venn(otus)
venn
aliens <- attr(venn, "intersections")$A
aliens
alien.list <- list(aliens, Compost.otus)
alien.list
aliens.venn <- venn(alien.list)
aliens.venn
aliens.detected <- attr(aliens.venn,"intersections")$`A:B`
incubated <- prune_samples(sample_data(incubation.physeq)$treatment %in% c("Compost"), incubation.physeq) %>%
  filter_taxa(function(x) sum(x) > 5, T) %>%
  transform_sample_counts(function(x) x / sum(x))
incubated
incubated.aliens <- prune_taxa(aliens.detected, incubated) 
incubated.aliens
sample.order <- as.data.frame(sample_data(incubated.aliens)) %>%
  arrange(day, replication) %>%
  select(i_id) %>%
  remove_rownames() 
sample.order
alien.heatmap <- plot_heatmap(incubated.aliens, sample.label = "day", taxa.order= "Phylum", taxa.label = "Genus",
                                  sample.order = as.character(sample.order$i_id), 
                                  low = "#66CCFF", high = "#000033", na.value = "white")
alien.heatmap

otus <- list(AlfalfaAmend.otus, Control.otus)
venn <- venn(otus)
aliens <- attr(venn, "intersections")$A
alien.list <- list(aliens, Alfalfa.otus)
aliens.venn <- venn(alien.list)
aliens.detected <- attr(aliens.venn,"intersections")$`A:B`
incubated <- prune_samples(sample_data(incubation.physeq)$treatment %in% c("Alfalfa"), incubation.physeq) %>%
  filter_taxa(function(x) sum(x) > 5, T) %>%
  transform_sample_counts(function(x) x / sum(x))
incubated.aliens <- prune_taxa(aliens.detected, incubated) 
sample.order <- as.data.frame(sample_data(incubated.aliens)) %>%
  arrange(day, replication) %>%
  select(i_id) %>%
  remove_rownames() 
alien.heatmap <- plot_heatmap(incubated.aliens, sample.label = "day", taxa.order= "Phylum", taxa.label = "Genus",
                              sample.order = as.character(sample.order$i_id), 
                              low = "white", high = "red", na.value = "gray")
alien.heatmap

