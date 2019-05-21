library(phylosmith)
library(phyloseq)

inc.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

tree <- read.tree("data/tree.nwk")

inc.physeq <- merge_phyloseq(inc.physeq, tree)

rm(tree)
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

rm(data)

inc.physeq <- subset_samples(inc.physeq, day %in% c("7",
                                                 "14",
                                                 "21",
                                                 "35",
                                                 "49",
                                                 "97"))

no.unclass <- subset_taxa(inc.physeq, !Phylum=="Bacteria_unclassified")
no.unclass <- subset_taxa(no.unclass, !Genus=="Gp6_unclassified")

rare6k.physeq <- rarefy_even_depth(no.unclass, sample.size = 6000, rngseed = 3242343, verbose = F) %>%
  filter_taxa(function(x) sum(x) >= 3, T) %>%
  tax_glom(taxrank = "Genus")

rare6k.physeq

phylogeny_profile_ggplot(rare6k.physeq, classification = 'Phylum', 
                         treatment = c("treatment"), 
                         merge = F, 
                         relative_abundance = TRUE)

taxa_abundance_bars_ggplot(rare6k.physeq, classification = 'Phylum', 
                           treatment = c('treatment'), 
                           transformation = 'mean')

tsne_phyloseq_ggplot(rare6k.physeq, treatment = c('treatment'), perplexity = 3, circle = TRUE, labels = FALSE, colors = 'default')

nmds_phyloseq_ggplot(rare6k.physeq, c('treatment'), circle = TRUE, verbose = FALSE)
nmds_phyloseq_ggplot(rare6k.physeq, c('day'), circle = TRUE, verbose = FALSE)
