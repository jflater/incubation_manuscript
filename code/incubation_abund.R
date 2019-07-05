library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)

inc <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

inc.physeq <- subset_samples(inc, day %in% c("7",
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

# Normalization, depth cutoff based on rarefaction 
rare6k.physeq <- rarefy_even_depth(inc.physeq, sample.size = 6000,
                                   rngseed = 15879966) %>%
  filter_taxa(function(x) sum(x) >= 2, T) 

RelativeAbundanceDf <- function(physeq) {
  physeq %>% tax_glom(taxrank = "Phylum") %>% transform_sample_counts(function(x) {
    x/sum(x)
  }) %>% psmelt() %>% filter(Abundance > 0.02) %>% arrange(Phylum)
}

# Function to plot relative abundance
PlotRelativeAbundance <- function(df) {
  ggplot(df, aes(x = as.factor(treatment), y = Abundance, fill = Phylum)) + 
    geom_bar(stat = "identity") +
    theme(axis.title.x = element_blank()) + 
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) + 
    ylab("Phylogenetic distribution\nof OTUs with abundance > 2%") 
}

treatment_names <- c(
  `1` = "Alfalfa",
  `2` = "Mix",
  `3` = "Compost",
  `4` = "Reference"
)
day_names <- c(
  `1` = "7",
  `2` = "14",
  `3` = "21",
  `4` = "35",
  `5` = "49",
  `6` = "97"
)

rare.merged <- merge_samples(rare6k.physeq, "TreatmentAndDay")

sample_data(rare.merged)$TreatmentAndDay <- levels(sample_data(rare6k.physeq)$TreatmentAndDay)

relllll <- PlotRelativeAbundance(RelativeAbundanceDf(rare.merged)) +
  facet_grid(~ day, labeller = labeller(day = as_labeller(day_names))) +
  scale_x_discrete(labels = treatment_names) +
  rotate_x_text(angle = 45)
relllll 

png("Figures/rela_abund.png",height=3,width=4,units='in',res=600)
relllll +
  scale_fill_viridis_d(option = "viridis", aesthetics = "fill")
dev.off()

