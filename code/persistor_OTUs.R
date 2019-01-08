library(tidyverse)
library(gplots)

# Read in whole phyloseq object
incubation.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")
levels(sample_data(incubation.physeq)$treatment)
######
# Working on function to get a list of OTUs from a sample type
GetOTUs <- function(physeq, samples) {
    prune_samples(sample_data(physeq)$treatment %in% c(samples), physeq) %>%
    filter_taxa(function(x) sum(x) > 0, T) %>%
    tax_table() %>%
    row.names()
}

incubated.Alfalfa.otus <- GetOTUs(incubation.physeq, c("Alfalfa"))
incubated.Compost.otus <- GetOTUs(incubation.physeq, c("Compost"))
incubated.CompAlfa.otus <- GetOTUs(incubation.physeq, c("CompAlfa"))
incubated.Control.otus <- GetOTUs(incubation.physeq, c("Control"))
incubated.AlfalfaAmend.otus <- GetOTUs(incubation.physeq, c("AlfalfaAmend"))
incubated.CompostAmend.otus <- GetOTUs(incubation.physeq, c("CompostAmend"))
incubated.AlfalfaSoil.otus <- GetOTUs(incubation.physeq, c("AlfalfaSoil"))

otus <- list(incubated.AlfalfaAmend.otus, incubated.CompostAmend.otus, incubated.AlfalfaSoil.otus)
venn <- venn(otus)
#####

# OTUs from the alfalfa amendment, not in soil
a <- list(incubated.Control.otus)
b <- list(incubated.AlfalfaSoil.otus)
c <- mapply(c, a, b, SIMPLIFY = F)
soilplusrefotus <- list(incubated.Control.otus, incubated.AlfalfaSoil.otus)

alfminussoil <- list()
only.alfa.input <- attr(venn, "intersections")$`A`

# Alfalfa treated soil and alfalfa amendment
alfa.otus <- list(only.alfa.input, alfala.physeq.input)
alf.venn <- venn(alfa.otus)

ex2 = prune_taxa(attr(alf.venn, "intersections")$`A`, alfala.physeq) %>%
  filter_taxa(function(x) sum(x) > 0, T)


  
