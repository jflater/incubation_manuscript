library(ggpubr)
# melt and glom for ggplot

tree <- read.tree("data/tree.nwk")
tree
inc_phy <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

# Add the tree file to the phyloseq object
inc_phy <- merge_phyloseq(inc_phy, tree)
inc_phy
inc.physeq <- subset_samples(inc_phy, day %in% c("7",
                                                 "14",
                                                 "21",
                                                 "35",
                                                 "49",
                                                 "97"))

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
sample_data(inc.physeq)

no.unclass <- subset_taxa(inc.physeq, !Phylum=="Bacteria_unclassified")
no.unclass <- subset_taxa(no.unclass, !Genus=="Gp6_unclassified")
rare6k.physeq <- rarefy_even_depth(physeq = no.unclass, sample.size = 6000, rngseed = 3242343, verbose = F)
rare6k.physeq

PhylumRelativeAbundanceDf <- rare6k.physeq %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {
    x/sum(x)
  })

phyla <- c("Verrucomicrobia", "Proteobacteria", "Acidobacteria")

# # get top 10 taxa
# phylum.control <- subset_samples(PhylumRelativeAbundanceDf, treatment %in% c("Reference")) %>%
#   prune_taxa(phyla)
# phylum.alfalfa <- subset_samples(PhylumRelativeAbundanceDf, treatment %in% c("Alfalfa"))
# phylum.compalfa <- subset_samples(PhylumRelativeAbundanceDf, treatment %in% c("Mix"))
# phylum.compost <- subset_samples(PhylumRelativeAbundanceDf, treatment %in% c("Compost")) 
# 
# 
# 
# ContOTUs <- names(sort(taxa_sums(phylum.control), TRUE)[1:5])
# AlfOTUs<- names(sort(taxa_sums(phylum.alfalfa), TRUE)[1:5])
# CompAlfOTUs <- names(sort(taxa_sums(phylum.compalfa), TRUE)[1:5])
# CompOTUs <- names(sort(taxa_sums(phylum.compost), TRUE)[1:5])
# 
# phylumotus <- c(AlfOTUs, CompAlfOTUs, CompAlfOTUs, CompOTUs)

inc10phylum <- subset_taxa(PhylumRelativeAbundanceDf, Phylum %in% c("Verrucomicrobia", "Proteobacteria", "Acidobacteria")) %>%
  psmelt()

# plot
#ggplot(inc10phylum, aes(x = treatment, y = Abundance), color = treatment) +
#  facet_grid(Phylum~day, scales = "free") +
#  geom_boxplot(aes(color = treatment), position = "dodge") 

boxplot <- ggboxplot(data = inc10phylum, x = "treatment", y = "Abundance", color = "treatment", legend = "none", ggtheme = theme_bw()) +
  rotate_x_text(angle = 45) +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "Reference") +
  facet_grid(Phylum ~ day, scales = "free")
# Use ggboxplot
boxplot + scale_color_viridis(discrete = T, option = "viridis") + ggplot2::theme_bw()
png("Figures/relaboxplot.png",height=8,width=10,units='in',res=600)
boxplot + scale_color_viridis(discrete = T, option = "viridis") + ggplot2::theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) + 
  ylab("Relative abundance") +
  xlab("Amended microcosms") +
  theme(legend.position = "none")
dev.off()
png("Figures/relaboxplot.png",height=8,width=10,units='in',res=600)
dev.off()