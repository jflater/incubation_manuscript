library(tidyverse)
library(gplots)

# Read in whole phyloseq object
incubation.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

######
# Working on function to get a list of OTUs from a sample type
GetOTUs <- function(z, y) {
    subset_samples(z, treatment == as.character(y)) %>%
    filter_taxa(function(x) sum(x) > 0, T) %>%
    tax_table() %>%
    row.names()
}

# Testing the function
alfalfa.amendment.otus <- GetOTUs(incubation.physeq, "Alfalfa")
#####

# Make new object with only samples from the inputs
inputs <- subset_samples(incubation.physeq, treatment %in% c("AlfalfaAmend", "AlfalfaSoil", "CompostAmend"))

# Get incubated alfalfa samples
alfala.incubated.otus <- subset_samples(incubation.physeq, treatment == "Alfalfa") %>%
  filter_taxa(function(x) sum(x) > 0, T) %>%
  tax_table() %>%
  row.names()

# Get OTUs from alfalfa amendment
alfalfa.amendment.otus <- subset_samples(inputs, treatment == "AlfalfaAmend") %>%
  filter_taxa(function(x) sum(x) > 0, T) %>%
  tax_table() %>%
  row.names()

# Get OTUs from compost amendment
compost.amendment.otus <- subset_samples(inputs, treatment == "CompostAmend") %>%
  filter_taxa(function(x) sum(x) > 0, T) %>%
  tax_table() %>%
  row.names()

# Get OTUs from the starting soil
soil.amendment.otus <- subset_samples(inputs, treatment == "AlfalfaSoil") %>%
  filter_taxa(function(x) sum(x) > 0, T) %>%
  tax_table() %>%
  row.names()

# Make list of OTUs to use in venn()
otus <- list(alfalfa.amendment.otus, compost.amendment.otus, soil.amendment.otus)
venn <- venn(otus)

only.alfa.input <- attr(venn, "intersections")$`A`

# Alfalfa treated soil and alfalfa amendment
alfa.otus <- list(only.alfa.input, alfala.physeq.input)
alf.venn <- venn(alfa.otus)

ex2 = prune_taxa(attr(alf.venn, "intersections")$`A`, alfala.physeq) %>%
  filter_taxa(function(x) sum(x) > 0, T)


  
