library(phyloseq)
library(vegan)
library(tidyverse)
library(viridis)
# This script used to generate relative abundance bar plots for the incubation inputs
# tps is variable name for the temporary physeq object
# Read in the phyloseq object, raw object from mothur, see code file mothur_to_phyloseq.R
inc.raw.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

# We will be characterizing the ammendments and starting soil, exclude all else
inc.inputs <- subset_samples(inc.raw.physeq, treatment %in% c("AlfalfaAmend", 
                                                              "AlfalfaSoil",
                                                              "CompostAmend"))
no.unclass <- subset_taxa(inc.inputs, !Phylum=="Bacteria_unclassified")
no.unclass <- subset_taxa(no.unclass, !Genus=="Gp6_unclassified")
tps = no.unclass

# remove unused data
rm(inc.inputs)
rm(inc.raw.physeq)
df <- as.data.frame(sample_data(tps)) # Put sample_data into a ggplot-friendly data.frame

# We need to add data for inorganic nitrogen concentration and C:N ratio. 
# Inorganic N can be calculated by the sumation of `NH3` and `NO3` and C:N is `C_flash` / `N_flash`:
# Pull out sample data and make a data.frame, pipe that into mutate() to create the new columns
df <- df %>%
  mutate(C_N = C_flash / N_flash, Inorganic_N = NH3 + NO3) 

rownames(df) <- df$i_id
sample_data(tps) <- df
# Sample names have start with a 2, we need to change these
names <- gsub(pattern = "2mm\\.", replacement = "", x = sample_names(tps))
sample_names(tps) <- names
sample_names(tps)

# Make table of input characteristics
colnames(df)
df
td <- df %>%
  select(N_flash, C_flash, NH3, NO3, Inorganic_N, C_N) 
td

d <- td[c("2mm.dry.alfalfa.soil", "2mm.dry.compost7", "2mm.dry.alfalfa.plant3"),] 
d
row.names(d) <- c("Soil", "Compost", "Alfalfa")

#td.table <- kable(d)
# table of characteristics
#td.table

#Let's check some of the of taxa counts from theses samples
# Minimum sample size
min(taxa_sums(otu_table(tps)))

# Remove singletons and 0 counts
tps <- filter_taxa(tps, function(x) sum(x) >= 2, T)
min(taxa_sums(otu_table(tps)))

tps <- rarefy_even_depth(tps, rngseed = 423423423)
physeq.dist <- vegdist(t(data.frame(otu_table(tps))), method = "bray")

PCoA_ord <- ordinate(
physeq = tps, 
method = "PCoA",
distance = physeq.dist
)

input.PCoA <- plot_ordination(physeq = tps, ordination = PCoA_ord,
                              color = "treatment",
                              title = "Input PCoA") + 
geom_point(size = 1.5) +
theme_bw() +
theme(legend.title = element_blank())

input.PCoA 

# Stats for input.PCoA
#input.PCoA.stats <- adonis(physeq.dist ~ treatment, data = data.frame(sample_data(tps)))
#png("Figures/input_PCoA_stats.png",height=2,width=9,units='in',res=300)
#kable(input.PCoA.stats$aov.tab)
#dev.off()


# Relative Abundance in inputs
# Put phyloseq object into a df with .02% phylum (glomed at phylum level)
RelativeAbundanceDf <- function(physeq) {
  physeq %>% tax_glom(taxrank = "Phylum") %>% transform_sample_counts(function(x) {
    x/sum(x)
  }) %>% psmelt() %>% arrange(Phylum)
}

treatment_names <- c(
  "AlfalfaAmend" = "Alfalfa amendment",
  "AlfalfaSoil" = "Starting soil",
  "CompostAmend" = "Compost amendment")

# Function to plot relative abundance
PlotRelativeAbundance <- function(df) {
  ggplot(df, aes(x = as.factor(treatment), y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", color = "black") +
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) + 
  ylab("Relative abundance of phyla") +
  xlab("Incubation inputs")
}

tps.merged <- merge_samples(tps, "treatment")
sample_data(tps.merged)$treatment <- levels(sample_data(tps)$treatment)
ggtheme = theme_bw()
plot <- PlotRelativeAbundance(RelativeAbundanceDf(tps.merged)) + 
  scale_x_discrete(labels = treatment_names) 
png("Figures/rela_abund_input.png",height=4,width=6,units='in',res=300)
plot + scale_fill_viridis(discrete = T, option = "viridis") + ggplot2::theme_bw()
dev.off()
plot + scale_fill_viridis(discrete = T, option = "viridis") + ggplot2::theme_bw()
