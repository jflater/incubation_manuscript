incubation.physeq <- readRDS("data/incubation_physeq_Aug18.RDS")

inputs <- subset_samples(incubation.physeq, treatment %in% c("AlfalfaAmend", "AlfalfaSoil", "CompostAmend"))

alfala.physeq <- subset_samples(incubation.physeq, treatment == c("Alfalfa")) %>%
  filter_taxa(function(x) sum(x) > 0, T) 

alfala.physeq.input <- alfala.physeq %>% tax_table() %>%
  row.names()

alfalfa.input <- subset_samples(inputs, treatment == "AlfalfaAmend") %>%
  filter_taxa(function(x) sum(x) > 0, T) %>%
  tax_table() %>%
  row.names()

compost.input <- subset_samples(inputs, treatment == "CompostAmend") %>%
  filter_taxa(function(x) sum(x) > 0, T) %>%
  tax_table() %>%
  row.names()

soil <- subset_samples(inputs, treatment == "AlfalfaSoil") %>%
  filter_taxa(function(x) sum(x) > 0, T) %>%
  tax_table() %>%
  row.names()

otus <- list(alfalfa.input, compost.input, soil)
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


  
