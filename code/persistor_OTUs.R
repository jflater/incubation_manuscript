library(tidyverse)
library(gplots)
# Read in whole phyloseq object
incubation.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

GetOTUs <- function(physeq, samplestogetotusfrom) {
  subset_samples(physeq, treatment == "samplestogetotusfrom") %>%
    filter_taxa(function(x) sum(x) > 0, T) %>%
    tax_table() %>%
    row.names()
}

alfalfa.amendment.otus <- GetOTUs(incubation.physeq, "Alfalfa")

# Make new object with only samples from the inputs
inputs <- subset_samples(incubation.physeq, treatment %in% c("AlfalfaAmend", "AlfalfaSoil", "CompostAmend"))

# Get incubated alfalfa samples
alfala.physeq <- subset_samples(incubation.physeq, treatment == c("Alfalfa")) %>%
  filter_taxa(function(x) sum(x) > 0, T) 

# get OTUs from alfalfa samples
alfala.incubated.otus <- alfala.physeq %>% tax_table() %>%
  row.names()

alfalfa.amendment.otus <- subset_samples(inputs, treatment == "AlfalfaAmend") %>%
  filter_taxa(function(x) sum(x) > 0, T) %>%
  tax_table() %>%
  row.names()

compost.amendment.otus <- subset_samples(inputs, treatment == "CompostAmend") %>%
  filter_taxa(function(x) sum(x) > 0, T) %>%
  tax_table() %>%
  row.names()

soil.amendment.otus <- subset_samples(inputs, treatment == "AlfalfaSoil") %>%
  filter_taxa(function(x) sum(x) > 0, T) %>%
  tax_table() %>%
  row.names()

otus <- list(alfalfa.amendment.otus, compost.amendment.otus, soil.amendment.otus)
venn <- venn(otus)

common.otus <- attr(venn,"intersections")$`A:B:C`
only.alfa.input <- attr(venn, "intersections")$`A`
ex1 = prune_taxa(common.otus, inputs)
tax_table(ex1)
plot_heatmap(ex1)

# Alfalfa treated soil and alfalfa amendment
alfa.otus <- list(only.alfa.input, alfala.physeq.input)
alf.venn <- venn(alfa.otus)

ex2 = prune_taxa(attr(alf.venn, "intersections")$`A`, alfala.physeq) %>%
  filter_taxa(function(x) sum(x) > 0, T)


  
