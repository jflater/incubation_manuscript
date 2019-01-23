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


# Venn diagram for identifying OTUs unique to alfalfa amendment that were not detected in the incubated control soils
otus <- list(AlfalfaAmend.otus, Control.otus)
venn <- venn(otus)

# OTUs in alfalfa amendment but not detected in incubated control samples are candidates for "Alien" status, or thos otus
# that transfer from the amendment to the soil and were not present prior amendment of soil
alf.aliens <- attr(venn, "intersections")$A

# compare aliens to incubated alfalfa samples
aliens <- list(alf.aliens, Alfalfa.otus)
aliens.venn <- venn(aliens)

# aliens detected!
alf.aliens.detected <- attr(aliens.venn,"intersections")$`A:B`


  
