library(ggpubr)
# melt and glom for ggplot
PhylumRelativeAbundanceDf <- rare6k.physeq %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {
    x/sum(x)
  })

# get top 10 taxa
phylum.control <- subset_samples(PhylumRelativeAbundanceDf, treatment %in% c("Reference"))
phylum.alfalfa <- subset_samples(PhylumRelativeAbundanceDf, treatment %in% c("Alfalfa"))
phylum.compalfa <- subset_samples(PhylumRelativeAbundanceDf, treatment %in% c("Mix"))
phylum.compost <- subset_samples(PhylumRelativeAbundanceDf, treatment %in% c("Compost")) 

ContOTUs <- names(sort(taxa_sums(phylum.control), TRUE)[1:5])
AlfOTUs<- names(sort(taxa_sums(phylum.alfalfa), TRUE)[1:5])
CompAlfOTUs <- names(sort(taxa_sums(phylum.compalfa), TRUE)[1:5])
CompOTUs <- names(sort(taxa_sums(phylum.compost), TRUE)[1:5])

phylumotus <- c(AlfOTUs, CompAlfOTUs, CompAlfOTUs, CompOTUs)

inc10phylum <- prune_taxa(phylumotus, PhylumRelativeAbundanceDf) %>%
  psmelt()

# plot
ggplot(inc10phylum, aes(x = treatment, y = Abundance), color = treatment) +
  facet_grid(Phylum~day, scales = "free") +
  geom_boxplot(aes(color = treatment), position = "dodge") 

boxplot <- ggboxplot(data = inc10phylum, x = "treatment", y = "Abundance", color = "treatment", legend = "none") +
  rotate_x_text(angle = 45) +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "Reference") +
  facet_grid(Phylum ~ day, scales = "free")
# Use ggboxplot

png("Figures/relaboxplot.png",height=9,width=9,units='in',res=600)
boxplot
dev.off()