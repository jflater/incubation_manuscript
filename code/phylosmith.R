library(phylosmith)
library(phyloseq)
library(tidyverse)
library(ape)
library(viridis)
library(vegan)
library(kableExtra)
theme_set(theme_bw())
inc <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

#tree <- read.tree("data/tree.nwk")

#inc.physeq <- merge_phyloseq(inc.physeq, tree)

#rm(tree)
#Rename treatments to more informative titles
data <- data.frame(sample_data(inc)) %>% 
  mutate(treatment = recode(treatment,
                            'Control' = 'Reference',
                            'CompAlfa' = 'Mix')) %>% 
  mutate(C_N = C_flash / N_flash, Inorganic_N = NH3 + NO3) %>%
  mutate(TreatmentAndDay = paste(treatment, day))

rownames(data) <- data$i_id
sample_data(inc) <- data
sample_data(inc)$day <- as.factor(sample_data(inc)$day)

rm(data)
colnames(sample_data(inc))
levels(inc@sam_data$treatment)

inc.inputs <- subset_samples(inc, treatment %in% c("AlfalfaAmend", "CompostAmend", "AlfalfaSoil")) %>% rarefy_even_depth(rngseed = 3242343, verbose = F) %>%
  filter_taxa(function(x) sum(x) >= 3, T)

inc.physeq <- subset_samples(inc, day %in% c("7",
                                                 "14",
                                                 "21",
                                                 "35",
                                                 "49",
                                                 "97"))

no.unclass <- subset_taxa(inc.physeq, !Phylum=="Bacteria_unclassified")
no.unclass <- subset_taxa(no.unclass, !Genus=="Gp6_unclassified")

rare6k.physeq <- rarefy_even_depth(no.unclass, sample.size = 6000, rngseed = 3242343, verbose = F) %>%
  filter_taxa(function(x) sum(x) >= 3, T)

head(tax_table(rare6k.physeq))
rare6k.physeq

phylogeny_profile_ggplot(rare6k.physeq, classification = 'Phylum', 
                         treatment = c("treatment"), 
                         merge = F, 
                         relative_abundance = TRUE)

taxa_abundance_bars_ggplot(rare6k.physeq, classification = 'Phylum', 
                           treatment = c('treatment'), 
                           transformation = 'mean')

png("Figures/tsne_treatment.png",height=4,width=4,units='in',res=600)
tsne_phyloseq_ggplot(rare6k.physeq, treatment = c('treatment'), perplexity = 10, circle = TRUE, colors = 'default') +
  scale_fill_viridis(discrete = T, option = "viridis") + ggplot2::theme_bw() +
  guides(fill=guide_legend(title="Treatment"))
dev.off()

# It would be nice to be able to label points that fall outside the elipse with the sample ID for easier outlier analysis
nmds_phyloseq_ggplot(rare6k.physeq, c('treatment'), circle = TRUE, verbose = FALSE)
nmds_phyloseq_ggplot(rare6k.physeq, c('day'), circle = TRUE, verbose = FALSE)


# Input ordinations
tsne_phyloseq_ggplot(inc.inputs, treatment = c('treatment'), perplexity = 10, circle = TRUE, colors = 'default') +
  scale_fill_viridis(discrete = T, option = "viridis") + ggplot2::theme_bw() +
  guides(fill=guide_legend(title="Treatment"))

####

PCoA.phy <- ordinate(inc.inputs, method = "PCoA", distance = "bray")
PCoA.treatment <- plot_ordination(inc.inputs, PCoA.phy, type = "samples", color = "treatment") +
  scale_fill_viridis_d(option = "viridis", aesthetics = "colour") +
  labs(color = "Input material") +
  ggtitle("PCoA ordination of the\nbray-curtis dissimilarity")
png("Figures/input_PCOA.png",height=4,width=4,units='in',res=600)
PCoA.treatment
dev.off()
d <- vegdist(t(data.frame(otu_table(inc.inputs))))
d.table <- adonis2(formula = d ~ treatment, data = data.frame(sample_data(inc.inputs)))
d.table
c_html=htmlTable(round(d.table, digits = 3),escape.html=FALSE)
c_html
png("Figures/input_PCOA_stats.png",height=4,width=4,units='in',res=600)
c_html
dev.off()
