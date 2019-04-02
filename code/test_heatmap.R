
physeq <- incubation.physeq
otus <- list(AlfalfaAmend.otus, Control.otus) 
print("Looking for aliens between amendment and control soil")
venn <- venn(otus)
alf.aliens <- attr(venn, "intersections")$A
aliens <- list(alf.aliens, Alfalfa.otus)
print("Detecting aliens in amended soil")
aliens.venn <- venn(aliens)
aliens.detected <- attr(aliens.venn,"intersections")$`A:B`
rare.merged <- merge_samples(physeq, "TreatmentAndDay")

sample_data(rare.merged)$TreatmentAndDay <- levels(sample_data(physeq)$TreatmentAndDay)

incubated <- prune_samples(sample_data(rare.merged)$treatment %in% c("Alfalfa"), rare.merged) %>%
  filter_taxa(function(x) sum(x) > 0, T) %>%
  transform_sample_counts(function(x) x / sum(x)) 
incubated.aliens <- prune_taxa(aliens.detected, incubated) %>%
  tax_glom(taxrank = "Genus")
#Ask Schuyler about this, 
#test <- conglomerate_samples(incubated.aliens, treatment = "treatment", subset = "Alfalfa", merge_on = "day")
test <- incubated.aliens
sample.order <- as.data.frame(sample_data(test)) %>%
  arrange(day, replication) %>%
  select(i_id) %>%
  remove_rownames()   
alien.heatmap <- plot_heatmap(test, sample.order = "day", taxa.order= "Phylum", taxa.label = "Genus",
                              low = "yellow", high = "red", na.value = "gray") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) 
alien.heatmap

