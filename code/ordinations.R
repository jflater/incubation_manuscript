if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("microbiome")

library(phyloseq)
library(vegan)
library(microbiome)
library(tidyverse)
inc.raw.physeq <- readRDS("Data/incubation_physeq_Aug18.RDS")

inc.physeq <- subset_samples(inc.raw.physeq, day %in% c("7",
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

# Normalization, depth cutoff based on ward.D2faction 
rare6k.physeq <- rarefy_even_depth(no.unclass, sample.size = 6000,
                                   rngseed = 15879966) %>%
  filter_taxa(function(x) sum(x) >= 3, T) 



phylum.sum = tapply(taxa_sums(rare6k.physeq), tax_table(rare6k.physeq)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
rare6k.physeq.top5 = prune_taxa((tax_table(rare6k.physeq)[, "Phylum"] %in% top5phyla), rare6k.physeq)

####NMDS
#######################
tmp.ps <- rare6k.physeq
d = vegdist(t(data.frame(otu_table(rare6k.physeq))),
            method = "bray")

NMDS.phy <- metaMDS(t(data.frame(otu_table(tmp.ps))), distance = "bray", k = 3, autotransform = F, trymax = 100)
                   
NMDS.treatment<- plot_ordination(tmp.ps, NMDS.phy, type = "samples", color = "day")

png("Figures/NMDS_incubation_axes12.png",height=4,width=9,units='in',res=600)
NMDS.treatment + stat_ellipse(geom = "polygon", type = "norm",
                   alpha = 0.2, aes(fill = day)) + 
  ggtitle("Axes 1 and 2, stress = 0.081")
dev.off()

NMDS.treatment.13 <- plot_ordination(tmp.ps, NMDS.phy, type = "samples", color = "day", axes = c(1, 3))

png("Figures/NMDS_incubation_axes13.png",height=4,width=9,units='in',res=600)
NMDS.treatment.13 + stat_ellipse(geom = "polygon", type = "norm",
                   alpha = 0.2, aes(fill = day)) +
  ggtitle("Axes 1 and 3, stress = 0.081")
dev.off()


NMDS.treatment<- plot_ordination(tmp.ps, NMDS.phy, type = "samples", color = "treatment")

png("Figures/NMDS_incubation_treatmentaxes12.png",height=4,width=9,units='in',res=600)
NMDS.treatment + stat_ellipse(geom = "polygon", type = "norm",
                              alpha = 0.2, aes(fill = treatment)) + 
  ggtitle("Axes 1 and 2, stress = 0.081")
dev.off()

NMDS.treatment.13 <- plot_ordination(tmp.ps, NMDS.phy, type = "samples", color = "treatment", axes = c(1, 3))

png("Figures/NMDS_incubation_treatmentaxes13.png",height=4,width=9,units='in',res=600)
NMDS.treatment.13 + stat_ellipse(geom = "polygon", type = "norm",
                                 alpha = 0.2, aes(fill = treatment)) +
  ggtitle("Axes 1 and 3, stress = 0.081")
dev.off()

PCoA.phy <- ordinate(tmp.ps, method = "PCoA",
                    distance = d)
PCoA.treatment <- plot_ordination(tmp.ps, PCoA.phy, type = "samples", color = "treatment")
PCoA.treatment <- PCoA.treatment + stat_ellipse(geom = "polygon", type = "norm", alpha = 0.3, aes(fill = treatment))

png("Figures/PCoA.treatment.png",height=4,width=9,units='in',res=600)
PCoA.treatment
dev.off()

PCoA.day <- plot_ordination(tmp.ps, PCoA.phy, type = "samples", color = "day")
PCoA.day <- PCoA.day + stat_ellipse(geom = "polygon", type = "norm", alpha = 0.3, aes(fill = day))

png("Figures/PCoA.day.png",height=4,width=9,units='in',res=600)
PCoA.day
dev.off()

PCoA.day.34 <- plot_ordination(tmp.ps, PCoA.phy, type = "samples", color = "day", axes = c(3, 4))
PCoA.day.34 <- PCoA.day.34 + stat_ellipse(geom = "polygon", type = "norm", alpha = 0.3, aes(fill = day))

a.tnd <- adonis(formula = d ~ TreatmentAndDay, data = data.frame(sample_data(tmp.ps)))
adonis(formula = d ~ day, data = data.frame(sample_data(tmp.ps)))

png("Figures/adonis_incubation_TreatmentAndDay.png",height=4,width=9,units='in',res=600)
grid.table(a.tnd$aov.tab)
dev.off()

a.dt <- adonis(formula = d ~ day * treatment, data = data.frame(sample_data(tmp.ps)))
adonis(formula = d ~ treatment * day, data = data.frame(sample_data(tmp.ps)))

png("Figures/adonis_incubation_day*trt.png",height=2,width=9,units='in',res=600)
grid.table(a.dt$aov.tab)
dev.off()



########################
#### Alfalfa
inc <- subset_samples(rare6k.physeq, treatment %in% c("Alfalfa"))
inc
min(taxa_sums(inc))
inc <- filter_taxa(inc, function(x) sum(x) >= 2, T)

inc.dist <- vegdist(t(data.frame(otu_table(inc))), method = "bray")
cap_ord <- ordinate(
  physeq = inc, 
  method = "CAP",
  distance = inc.dist,
  formula = ~ MBC_mg.kg_per_dry_wt_soil + C_N + Inorganic_N
)
x <- anova.cca(cap_ord, by = "margin")
# CAP plot
cap_plot <- plot_ordination(
  physeq = inc, 
  ordination = cap_ord, 
  color = "day", 
  axes = c(1,2)
) + 
  geom_point(aes(colour = day), alpha = 0.8, size = 4) + 
  geom_point(colour = "grey90", size = 1.5) 

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = .8 * CAP1, 
                 y = .8 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
Alfalfa <- cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) +
  ggtitle("Alfalfa CAP")


########################
#### Compost
inc <- subset_samples(rare6k.physeq, treatment %in% c("Compost"))
inc
min(taxa_sums(inc))
inc <- filter_taxa(inc, function(x) sum(x) >= 2, T)

inc.dist <- vegdist(t(data.frame(otu_table(inc))), method = "bray")
cap_ord <- ordinate(
  physeq = inc, 
  method = "CAP",
  distance = inc.dist,
  formula = ~ MBC_mg.kg_per_dry_wt_soil + C_N + Inorganic_N
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = inc, 
  ordination = cap_ord, 
  color = "day", 
  axes = c(1,2)
) + 
  geom_point(aes(colour = day), alpha = 0.8, size = 4) + 
  geom_point(colour = "grey90", size = 1.5) 

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = .8 * CAP1, 
                 y = .8 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
Compost <- cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) +
  ggtitle("Compost CAP")

########################
#### Mix
inc <- subset_samples(rare6k.physeq, treatment %in% c("Mix"))
inc
min(taxa_sums(inc))
inc <- filter_taxa(inc, function(x) sum(x) >= 2, T)

inc.dist <- vegdist(t(data.frame(otu_table(inc))), method = "bray")
cap_ord <- ordinate(
  physeq = inc, 
  method = "CAP",
  distance = inc.dist,
  formula = ~ MBC_mg.kg_per_dry_wt_soil + C_N + Inorganic_N
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = inc, 
  ordination = cap_ord, 
  color = "day", 
  axes = c(1,2)
) + 
  geom_point(aes(colour = day), alpha = 0.8, size = 4) + 
  geom_point(colour = "grey90", size = 1.5) 

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = .8 * CAP1, 
                 y = .8 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
Mix <- cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) +
  ggtitle("Mix CAP")

########################
#### Control
inc <- subset_samples(rare6k.physeq, treatment %in% c("Reference"))
inc
min(taxa_sums(inc))
inc <- filter_taxa(inc, function(x) sum(x) >= 2, T)

inc.dist <- vegdist(t(data.frame(otu_table(inc))), method = "bray")
cap_ord <- ordinate(
  physeq = inc, 
  method = "CAP",
  distance = inc.dist,
  formula = ~ MBC_mg.kg_per_dry_wt_soil + C_N + Inorganic_N
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = inc, 
  ordination = cap_ord, 
  color = "day", 
  axes = c(1,2)
) + 
  geom_point(aes(colour = day), alpha = 0.8, size = 4) + 
  geom_point(colour = "grey90", size = 1.5) 

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = .8 * CAP1, 
                 y = .8 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
Control <- cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) +
  ggtitle("Reference CAP")
tiff("Figures/CAP_Control.tif",height=4,width=9,units='in',res=300)
Control
dev.off()  

tiff("Figures/CAP_Alfalfa.tif",height=4,width=9,units='in',res=300)
Alfalfa
dev.off()  

tiff("Figures/CAP_Compost.tif",height=4,width=9,units='in',res=300)
Compost
dev.off()

tiff("Figures/CAP_Mix.tif",height=4,width=9,units='in',res=300)
Mix
dev.off()  