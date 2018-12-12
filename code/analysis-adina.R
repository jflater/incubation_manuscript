setwd("~/Box Sync/Papers/args-sc/")
library(phyloseq)
library(plyr)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(tidyr)

abundance_limit = 0.0

arg_phy <- readRDS("arg_phy_updated.RDS")
#int_phy <- readRDS("int-phyloseq.RDS")
#ann_card <- read.delim(sep="\t", file="ann2.txt", row.names=1, header=FALSE)
#ann_new <- tax_table(as.matrix(ann_card))
#arg_phy2 <- phyloseq(sample_data(arg_phy), otu_table(arg_phy), ann_new)

GP.ord <- ordinate(arg_phy, "NMDS", "bray")
p1 = plot_ordination(arg_phy, GP.ord, color="Matrix", shape="Treatment", size=12)
p1 + theme_bw()

# RESULT 1 - ARGS IN MANURE AND SOIL BACKGROUNDS
# ARG_PHY HAS 355 GENES AND ARG2 (CARD DATABASE SUBSET) 84 GENES, SO LET'S GO WITH ARG_PHY
arg_manure <- subset_samples(arg_phy, Matrix == "Manure")
arg_manure2 <- prune_taxa(taxa_sums(arg_manure) > abundance_limit, arg_manure)  #170 genes in manure
arg_manure_melt <- psmelt(arg_manure2)
arg_manure_melt$norm <- arg_manure_melt$Abundance/arg_manure_melt$hk_count
arg_manure_melt.all <- ddply(arg_manure_melt, .(Sample, Day, OTU, Treatment), summarise, MEAN = mean(norm), SE = sd(norm)/sqrt(length(norm)))
arg_manure_melt.all$presence.absence <- ifelse(arg_manure_melt.all$MEAN > 0, 1, 0)
arg_manure_melt.all2 <- ddply(arg_manure_melt.all, .(OTU), summarise, TOT_SAMPLES = sum(presence.absence))
all_manure_samples <- subset(arg_manure_melt.all2, TOT_SAMPLES == 4) #100 genes
#all_manure_samples <- subset(arg_manure_melt.all2, TOT_SAMPLES >= 1) #170 genes
otus.in.all.manures <- all_manure_samples$OTU

manure_otus <- prune_taxa(otus.in.all.manures, arg_phy)
manure_otus <- subset_samples(manure_otus, Matrix == "Manure") 
# 100 genes in all 4 manure samples
manure_otus <- prune_taxa(taxa_sums(manure_otus) > abundance_limit, manure_otus) # changed for 0
manure_otus.melt <- psmelt(manure_otus)
manure_otus.melt$norm <- manure_otus.melt$Abundance/manure_otus.melt$hk_count
manure_dist <- ddply(manure_otus.melt, .(OTU, V2), summarise, MEAN = mean(norm), SE = sd(norm)/sqrt(length(norm)))
rank_order_otus <- manure_dist[order(-manure_dist$MEAN),]$OTU
manure_dist$OTU = factor(manure_dist$OTU, levels = rank_order_otus)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
manure_dist2 <- subset(manure_dist, MEAN >= 2.44)
p <-ggplot(manure_dist2,aes(OTU, MEAN, color=V2))+geom_point(stat = "identity") + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=1, size=10)) + scale_x_discrete(breaks=manure_dist$OTU, labels=manure_dist$V2)
p+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + xlab("Gene Name") + ylab("Average Gene Count per recA gene")

# Soil background
arg_soil <- subset_samples(arg_phy, Matrix == "Soil"& Treatment == "Control")
arg_soil2 <- prune_taxa(taxa_sums(arg_soil) > abundance_limit, arg_soil)  #93 genes in soil
arg_soil_melt <- psmelt(arg_soil2)
arg_soil_melt$norm <- arg_soil_melt$Abundance/arg_soil_melt$hk_count
arg_soil_melt.all <- ddply(arg_soil_melt, .(ID, Day, OTU, Treatment), summarise, MEAN = mean(norm), SE = sd(norm)/sqrt(length(norm)))
arg_soil_melt.all$presence.absence <- ifelse(arg_soil_melt.all$MEAN > 0, 1, 0)
arg_soil_melt.all2 <- ddply(arg_soil_melt.all, .(OTU), summarise, TOT_SAMPLES = sum(presence.absence))
all_soil_samples <- subset(arg_soil_melt.all2, TOT_SAMPLES == 4)  
#all_soil_samples <- subset(arg_soil_melt.all2, TOT_SAMPLES >= 1) 
otus.in.all.soils <- all_soil_samples$OTU
soil_otus <- prune_taxa(otus.in.all.soils, arg_phy)
# 24 genes in all 4 soil samples
soil_otus <- prune_taxa(taxa_sums(soil_otus) > abundance_limit, soil_otus) #24 OTUs
soil_otus.melt <- psmelt(soil_otus)
soil_otus.melt$norm <- soil_otus.melt$Abundance/soil_otus.melt$hk_count
soil_dist <- ddply(soil_otus.melt, .(OTU, V2), summarise, MEAN = mean(norm), SE = sd(norm)/sqrt(length(norm)))
rank_order_otus <- soil_dist[order(-soil_dist$MEAN),]$OTU
soil_dist$OTU = factor(soil_dist$OTU, levels = rank_order_otus)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
soil_dist2 <- subset(soil_dist, MEAN >= abundance_limit)
p <-ggplot(soil_dist2,aes(OTU, MEAN, color=V2))+geom_point(stat = "identity") + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=0, size=10)) + scale_x_discrete(breaks=soil_dist$OTU, labels=soil_dist$V2)
p+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + xlab("Gene Name") + ylab("Average Gene Count per recA gene")


# soil and manure ARGS together
same_otus <- intersect(taxa_names(soil_otus), taxa_names(manure_otus))  #6 when in only in one sample
manure_otus_in_environment <- prune_taxa(taxa_names(manure_otus), arg_phy)
manure_otus_in_environment2 <- subset_samples(manure_otus_in_environment, Treatment == "Control") 
manure_otus_in_environment2 <- subset_samples(manure_otus_in_environment2, Matrix != "Effluent") 
manure_dist <- psmelt(manure_otus_in_environment2)
manure_dist$norm <- manure_dist$Abundance/manure_dist$hk_count
manure_dist2 <- ddply(manure_dist, .(OTU, V2, Matrix), summarise, MEAN = mean(norm), SE = sd(norm)/sqrt(length(norm)))
manure_dist2$shared <- ifelse(manure_dist2$OTU %in% soil_dist$OTU, 1, 0)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
temp <- subset(manure_dist2, Matrix == "Manure")
rank_order_otus <- temp[order(-temp$MEAN),]$OTU
manure_dist2$OTU = factor(manure_dist2$OTU, levels = rank_order_otus)
p <-ggplot(manure_dist2,aes(OTU, MEAN, color=shared))+geom_point(stat = "identity") + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=0, size=10)) + scale_x_discrete(breaks=manure_dist$OTU, labels=manure_dist$V2)
p 


arg_phy2 <- subset_samples(arg_phy, Treatment == "Manure")
GP.ord <- ordinate(arg_phy2, "NMDS", "bray")
p1 = plot_ordination(arg_phy, GP.ord, shape="Matrix", color="Day")
p1 + theme_bw()


# manure otus, requires all OTUS to be in all 4 manure biological replicates
meta <- sample_data(arg_phy)
arg_manure <- subset_samples(arg_phy, Matrix == "Manure")
#int_manure <- subset_samples(int_phy, Matrix == "Manure")
arg_manure2 <- prune_taxa(taxa_sums(arg_manure) > abundance_limit, arg_manure)
#int_manure2 <- prune_taxa(taxa_sums(int_manure) > 0, int_manure)
arg_manure_melt <- psmelt(arg_manure2)
arg_manure_melt$norm <- arg_manure_melt$Abundance/arg_manure_melt$hk_count
arg_manure_melt.all <- ddply(arg_manure_melt, .(ID, Day, OTU, Treatment), summarise, MEAN = mean(norm), SE = sd(norm)/sqrt(length(norm)))
arg_manure_melt.all$presence.absence <- ifelse(arg_manure_melt.all$MEAN > 0, 1, 0)
arg_manure_melt.all2 <- ddply(arg_manure_melt.all, .(OTU), summarise, TOT_SAMPLES = sum(presence.absence))
all_manure_samples <- subset(arg_manure_melt.all2, TOT_SAMPLES == 4)
otus.in.all.manures <- all_manure_samples$OTU
manure_otus <- prune_taxa(otus.in.all.manures, arg_phy) #49 taxa in CARD, 100 taxa in total
#manure_otus <- prune_taxa(otus.in.all.manures, arg_phy2) #49 taxa in CARD, 100 taxa in total

manure_otus <- subset_samples(manure_otus, Matrix == "Manure")
manure_otus <- prune_taxa(taxa_sums(manure_otus) > abundance_limit, manure_otus)


# soil otus overlapping with 49 manure OTUs
#arg_soil <- subset_samples(arg_phy, Matrix == "Soil" & Treatment == "Manure")
#int_soil <- subset_samples(int_phy, Matrix == "Soil" & Treatment == "Manure")
#arg_soil2 <- prune_taxa(taxa_sums(arg_soil) > 0, arg_soil)
#int_soil2 <- prune_taxa(taxa_sums(int_soil) > 0, int_soil)
# arg_soil_melt <- psmelt(arg_soil2)
# arg_soil_melt.all <- ddply(arg_soil_melt, .(ID, Day, OTU, Treatment), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
# arg_soil_melt.all$presence.absence <- ifelse(arg_soil_melt.all$MEAN > 0, 1, 0)
# arg_soil_melt.all2 <- ddply(arg_soil_melt.all, .(OTU), summarise, TOT_SAMPLES = sum(presence.absence))
# all_soil_samples <- subset(arg_soil_melt.all2, TOT_SAMPLES == 4)
# otus.in.all.soils <- all_soil_samples$OTU


soil_otus <- prune_taxa(rownames(otu_table(manure_otus)), arg_phy)  #manure OTUs only
soil_otus_manured <- subset_samples(soil_otus, Treatment == "Control" & Matrix == "Soil") #select only control soil
soil_otus_manured2 <- prune_taxa(taxa_sums(soil_otus_manured) == 0, soil_otus_manured) # removes 61 OTUs that are NOT control soils
soil_manure_otus <- taxa_names(soil_otus_manured2) # OTUs in manure but not in soil
effluent_otus_manured <- subset_samples(soil_otus, Treatment == "Control" & Matrix == "Effluent")
effluent_otus_manured2 <- prune_taxa(taxa_sums(effluent_otus_manured) == 0, effluent_otus_manured)
effluent_manure_otus <- taxa_names(effluent_otus_manured2)  # OTUs in manure but not in effluent
manure_only <- intersect(soil_manure_otus, effluent_manure_otus)
soil_only <- setdiff(soil_manure_otus, effluent_manure_otus) #11 only in soil manured
eff_only <- setdiff(effluent_manure_otus, soil_manure_otus) #13 only in effluent manured (need to remove things that are in effluent control)

sw_otus <- prune_taxa(manure_only, arg_phy)
sw_otus.melt <- psmelt(sw_otus)
foo <- sw_otus.melt
foo$norm <- foo$Abundance/foo$hk_count
# Heatmap
foo$Matrix <- factor(foo$Matrix, levels = c("Manure", "Soil", "Effluent"))
temp2 <- subset(foo, OTU == "gb|FR734406|+|0-906|ARO:3003106|Erm(42)")
sample_order <- temp2[order(temp2$Matrix),]$Sample
foo$Sample <- factor(foo$Sample, levels = sample_order)

foo_temp <- subset(foo, Treatment == "Control" & Matrix == "Manure")
foo_temp2 <- ddply(foo_temp, .(OTU, V2), summarize, MEAN = mean(norm))
order <- foo_temp2[order(foo_temp2$V2, foo_temp2$MEAN),]$OTU
foo$Day = factor(foo$Day, levels=c("0", "10","24", "38", "59","80","108"))
foo$OTU = factor(foo$OTU, levels = order)
p <-ggplot(foo,aes(Sample,OTU,fill=norm))+
  geom_tile(color = "white", size = 0.25) 
p <-p + facet_grid(Treatment ~ Day,  scale = "free", space = "free")+scale_y_discrete(breaks=foo$OTU, labels=foo$V2) +
  scale_x_discrete(breaks = foo$Sample, labels = foo$Matrix)
p +scale_fill_gradient(low = "white", high = "red", trans = "log") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Heatmap Without controls
sw_otus <- prune_taxa(manure_only, arg_phy)
sw_otus <- subset_samples(sw_otus, Treatment != "Control" | Matrix == "Manure")
foo <- psmelt(sw_otus)
foo$norm <- foo$Abundance/foo$hk_count
foo$norm <- ifelse(foo$norm > abundance_limit, foo$norm, 0)
foo$Matrix <- factor(foo$Matrix, levels = c("Manure", "Soil", "Effluent"))
temp2 <- subset(foo, OTU == "gb|FR734406|+|0-906|ARO:3003106|Erm(42)")
sample_order <- temp2[order(temp2$Matrix, as.numeric(temp2$Day)),]$Sample
foo$Sample <- factor(foo$Sample, levels = sample_order)
foo_temp <- subset(foo, Treatment == "Control" & Matrix == "Manure")
foo_temp2 <- ddply(foo_temp, .(OTU, V2), summarize, MEAN = mean(norm))
order <- foo_temp2[order(foo_temp2$MEAN),]$OTU
foo$Matrix = factor(foo$Matrix, levels = c("Manure", "Soil", "Effluent"))
foo$Day = factor(foo$Day, levels=c("0", "10","24", "38", "59","80","108"))
foo$OTU = factor(foo$OTU, levels = order)

p <-ggplot(foo,aes(Sample,OTU,fill=norm))+
  geom_tile(color = "white", size = 0.25) 
p <-p + facet_grid(~ Matrix,  scale = "free", space = "free")+scale_y_discrete(breaks=foo$OTU, labels=foo$V2) +
  scale_x_discrete(breaks = foo$Sample, labels = foo$Day)
p +scale_fill_gradient(low = "white", high = "red", trans = "log") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


# SOIL ONLY

sw_otus <- prune_taxa(eff_only, arg_phy)
sw_otus.melt <- psmelt(sw_otus)
foo <- sw_otus.melt
foo$norm <- foo$Abundance/foo$hk_count
# Heatmap
foo$Matrix <- factor(foo$Matrix, levels = c("Manure", "Soil", "Effluent"))
temp2 <- subset(foo, OTU == "gb|M18896.2|+|206-2126|ARO:3000190|tetO")
sample_order <- temp2[order(temp2$Matrix),]$Sample
foo$Sample <- factor(foo$Sample, levels = sample_order)

foo_temp <- subset(foo, Treatment == "Control" & Matrix == "Manure")
foo_temp2 <- ddply(foo_temp, .(OTU, V2), summarize, MEAN = mean(norm))
order <- foo_temp2[order(foo_temp2$V2, foo_temp2$MEAN),]$OTU
foo$Day = factor(foo$Day, levels=c("0", "10","24", "38", "59","80","108"))
foo$OTU = factor(foo$OTU, levels = order)
p <-ggplot(foo,aes(Sample,OTU,fill=norm))+
  geom_tile(color = "white", size = 0.25) 
p <-p + facet_grid(Treatment ~ Day,  scale = "free", space = "free")+scale_y_discrete(breaks=foo$OTU, labels=foo$V2) +
  scale_x_discrete(breaks = foo$Sample, labels = foo$Matrix)
p +scale_fill_gradient(low = "white", high = "red", trans = "log") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Heatmap Without controls
sw_otus <- prune_taxa(manure_only, arg_phy)
sw_otus <- subset_samples(sw_otus, Treatment != "Control" | Matrix == "Manure")
foo <- psmelt(sw_otus)
foo$norm <- foo$Abundance/foo$hk_count
foo$norm <- ifelse(foo$norm > abundance_limit, foo$norm, 0)
foo$Matrix <- factor(foo$Matrix, levels = c("Manure", "Soil", "Effluent"))
temp2 <- subset(foo, OTU == "gb|FR734406|+|0-906|ARO:3003106|Erm(42)")
sample_order <- temp2[order(temp2$Matrix, as.numeric(temp2$Day)),]$Sample
foo$Sample <- factor(foo$Sample, levels = sample_order)
foo_temp <- subset(foo, Treatment == "Control" & Matrix == "Manure")
foo_temp2 <- ddply(foo_temp, .(OTU, V2), summarize, MEAN = mean(norm))
order <- foo_temp2[order(foo_temp2$MEAN),]$OTU
foo$Matrix = factor(foo$Matrix, levels = c("Manure", "Soil", "Effluent"))
foo$Day = factor(foo$Day, levels=c("0", "10","24", "38", "59","80","108"))
foo$OTU = factor(foo$OTU, levels = order)

p <-ggplot(foo,aes(Sample,OTU,fill=norm))+
  geom_tile(color = "white", size = 0.25) 
p <-p + facet_grid(~ Matrix,  scale = "free", space = "free")+scale_y_discrete(breaks=foo$OTU, labels=foo$V2) +
  scale_x_discrete(breaks = foo$Sample, labels = foo$Day)
p +scale_fill_gradient(low = "white", high = "red", trans = "log") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Bar graph
sw_otus <- prune_taxa(manure_only, arg_phy)
sw_otus <- subset_samples(sw_otus, Treatment != "Control" | Matrix == "Manure")
foo <- psmelt(sw_otus)
foo$norm <- foo$Abundance/foo$hk_count
#foo$norm <- ifelse(foo$norm > norm_limit, foo$norm, 0)
#foo$Matrix <- factor(foo$Matrix, levels = c("Manure", "Soil", "Effluent"))
#temp2 <- subset(foo, OTU == "gb|FR734406|+|0-906|ARO:3003106|Erm(42)")
#sample_order <- temp2[order(temp2$Matrix),]$Sample
#foo$Sample <- factor(foo$Sample, levels = sample_order)
#foo_temp <- subset(foo, Treatment == "Control" & Matrix == "Manure")
foo2 <- ddply(foo, .(OTU, V2, Matrix, Day), summarize, MEAN = mean(norm))
foo2$Presence.absence <- ifelse(foo2$MEAN > 0, 1, 0)
order <- foo_temp2[order(foo_temp2$MEAN),]$OTU
foo2 <- subset(foo2, Day == "24" | Day == "59" | Day == "108")
foo2$Day = factor(foo2$Day, levels=c("24", "59","108"))
foo2$Matrix = factor(foo2$Matrix, levels=c("Soil", "Effluent"))
foo$OTU = factor(foo$OTU, levels = order)
p <- ggplot(foo2,aes(x = Matrix, y = Presence.absence, fill = OTU))+geom_bar(stat = "identity") + facet_grid(~Day, scale = "free", space="free")
p + theme(legend.position = 'none')

# Heatmap - abundance curves
ordering <- foo[order(-foo$norm),]
ordering$rank_order <- seq.int(nrow(ordering))
p <-ggplot(ordering,aes(rank_order, log(norm)))+geom_point(stat = "identity")# + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=0, size=10)) 
p

# Making Table
arg_manure <- subset_samples(arg_phy, Matrix == "Manure")
arg_manure2 <- prune_taxa(taxa_sums(arg_manure) > abundance_limit, arg_manure)  #163 genes in manure
arg_manure_melt <- psmelt(arg_manure2)
arg_manure_melt.all <- ddply(arg_manure_melt, .(ID, Day, OTU, Treatment), summarise, MEAN = mean(norm), SE = sd(norm)/sqrt(length(norm)))
arg_manure_melt.all$presence.absence <- ifelse(arg_manure_melt.all$MEAN > 0, 1, 0)
arg_manure_melt.all2 <- ddply(arg_manure_melt.all, .(OTU), summarise, TOT_SAMPLES = sum(presence.absence))
all_manure_samples <- subset(arg_manure_melt.all2, TOT_SAMPLES == 4) #100 genes
all_manure_samples <- subset(arg_manure_melt.all2, TOT_SAMPLES >= 1) #163 genes
otus.in.all.manures <- all_manure_samples$OTU

manure_otus <- prune_taxa(otus.in.all.manures, arg_phy)
manure_otus_day24 <- subset_samples(manure_otus, Day == "24" & Matrix != "Manure")
manure_otus_day24 <- subset_samples(manure_otus_day24, Treatment == "Manure")
soil_24 <- subset_samples(manure_otus_day24, Matrix == "Soil")
soil_24 <- prune_taxa(taxa_sums(soil_24) > abundance_limit, soil_24) #130 OTUs 
soil_24 
eff_24 <- subset_samples(manure_otus_day24, Matrix == "Effluent")
eff_24 <- prune_taxa(taxa_sums(eff_24) > abundance_limit, eff_24) #88
eff_24
same_otus <- intersect(taxa_names(soil_24), taxa_names(eff_24)) #87
length(same_otus)

manure_otus <- prune_taxa(otus.in.all.manures, arg_phy)
manure_otus_day59 <- subset_samples(manure_otus, Day == "59" & Matrix != "Manure")
manure_otus_day59 <- subset_samples(manure_otus_day59, Treatment == "Manure")
soil_59 <- subset_samples(manure_otus_day59, Matrix == "Soil")
soil_59 <- prune_taxa(taxa_sums(soil_59) > abundance_limit, soil_59) #131
soil_59 
eff_59 <- subset_samples(manure_otus_day59, Matrix == "Effluent")
eff_59 <- prune_taxa(taxa_sums(eff_59) > abundance_limit, eff_59)  #92
eff_59
same_otus <- intersect(taxa_names(soil_59), taxa_names(eff_59)) #89
length(same_otus) 

manure_otus <- prune_taxa(otus.in.all.manures, arg_phy)
manure_otus_day108 <- subset_samples(manure_otus, Day == "108" & Matrix != "Manure")
manure_otus_day108 <- subset_samples(manure_otus_day108, Treatment == "Manure")
soil_108 <- subset_samples(manure_otus_day108, Matrix == "Soil")
soil_108 <- prune_taxa(taxa_sums(soil_108) > abundance_limit, soil_108) #126
soil_108 
eff_108 <- subset_samples(manure_otus_day108, Matrix == "Effluent")
eff_108 <- prune_taxa(taxa_sums(eff_108) > abundance_limit, eff_108) #119
eff_108
same_otus <- intersect(taxa_names(soil_108), taxa_names(eff_108)) #113 shared; 
length(same_otus)


### Without control otus
manure_otus <- prune_taxa(otus.in.all.manures, arg_phy)
manure_otus <- prune_taxa(manure_only, manure_otus)
manure_otus_day24 <- subset_samples(manure_otus, Day == "24" & Matrix != "Manure")
manure_otus_day24 <- subset_samples(manure_otus_day24, Treatment == "Manure")
soil_24 <- subset_samples(manure_otus_day24, Matrix == "Soil")
soil_24 <- prune_taxa(taxa_sums(soil_24) > abundance_limit, soil_24) #130 OTUs 
soil_24 
eff_24 <- subset_samples(manure_otus_day24, Matrix == "Effluent")
eff_24 <- prune_taxa(taxa_sums(eff_24) > abundance_limit, eff_24) #88
eff_24
same_otus <- intersect(taxa_names(soil_24), taxa_names(eff_24)) #87
length(same_otus)

manure_otus <- prune_taxa(otus.in.all.manures, arg_phy)
manure_otus <- prune_taxa(manure_only, manure_otus)
manure_otus_day59 <- subset_samples(manure_otus, Day == "59" & Matrix != "Manure")
manure_otus_day59 <- subset_samples(manure_otus_day59, Treatment == "Manure")
soil_59 <- subset_samples(manure_otus_day59, Matrix == "Soil")
soil_59 <- prune_taxa(taxa_sums(soil_59) > abundance_limit, soil_59) #131
soil_59 
eff_59 <- subset_samples(manure_otus_day59, Matrix == "Effluent")
eff_59 <- prune_taxa(taxa_sums(eff_59) > abundance_limit, eff_59)  #92
eff_59
same_otus <- intersect(taxa_names(soil_59), taxa_names(eff_59)) #89
length(same_otus) 

manure_otus <- prune_taxa(otus.in.all.manures, arg_phy)
manure_otus <- prune_taxa(manure_only, manure_otus)
manure_otus_day108 <- subset_samples(manure_otus, Day == "108" & Matrix != "Manure")
manure_otus_day108 <- subset_samples(manure_otus_day108, Treatment == "Manure")
soil_108 <- subset_samples(manure_otus_day108, Matrix == "Soil")
soil_108 <- prune_taxa(taxa_sums(soil_108) > abundance_limit, soil_108) #126
soil_108 
eff_108 <- subset_samples(manure_otus_day108, Matrix == "Effluent")
eff_108 <- prune_taxa(taxa_sums(eff_108) > abundance_limit, eff_108) #119
eff_108
same_otus <- intersect(taxa_names(soil_108), taxa_names(eff_108)) #113 shared; 
length(same_otus)

#checking differences
diff1 <- setdiff(taxa_names(soil_108), taxa_names(eff_108))
diff2 <- setdiff(taxa_names(eff_108), taxa_names(soil_108))
diff <- prune_taxa(diff2, arg_phy)
diff.melt <- psmelt(diff)
diff_dist <- ddply(diff.melt, .(OTU, Matrix, Treatment, Day), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
diff_dist <- subset(diff_dist, Matrix != "Manure")
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p <-ggplot(diff_dist,aes(OTU, MEAN, color=Treatment))+geom_point(stat = "identity") + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=0, size=10)) 
p + facet_grid(Day~Matrix)

#abundance curves
soil_manure <- subset_samples(arg_phy, Matrix == "Soil" & Treatment == "Manure")
soil_manure.melt <- psmelt(soil_manure)
soil_manure <- ddply(soil_manure.melt, .(OTU, Day), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
soil_manure <- subset(soil_manure, MEAN > 0)
ordering <- soil_manure[order(-soil_manure$MEAN),]
ordering$rank_order <- seq.int(nrow(ordering))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p <-ggplot(ordering,aes(rank_order, log(MEAN)))+geom_point(stat = "identity")# + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=0, size=10)) 
p


effluent_manure <- subset_samples(arg_phy, Matrix == "Effluent" & Treatment == "Manure")
effluent_manure.melt <- psmelt(effluent_manure)
effluent_manure <- ddply(effluent_manure.melt, .(OTU, Day), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
effluent_manure <- subset(effluent_manure, MEAN > 0)
ordering <- effluent_manure[order(-effluent_manure$MEAN),]
ordering$rank_order <- seq.int(nrow(ordering))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p <-ggplot(ordering,aes(rank_order, log10(MEAN)))+geom_point(stat = "identity")# + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=0, size=10)) 
p
# foo3 <- soil_otus.melt
# foo <- foo3 %>% mutate(new.day = ifelse(Treatment == "Control", "0-C", Day))
# foo <- foo %>% mutate(new.day2 = ifelse(Matrix == "Manure", "0-C-M", Day))
# soil_plus_controls <- foo

# 
# 
# # This is old code for circular heatmaps
# soil_dist <- ddply(soil_otus.melt, .(OTU, Day, V6), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
# soil_dist_day24 <- subset(soil_dist, Day == "24")
# soil_dist_day59<- subset(soil_dist, Day == "59")
# soil_dist_day108<- subset(soil_dist, Day == "108")
# foo <- merge(manure_dist, soil_dist_day24, by = "OTU")
# foo2 <- merge(foo,soil_dist_day59, by = "OTU" )
# foo3 <- merge(foo2, soil_dist_day108, by = "OTU")
# new_dist <- foo3[,c(1,2,3,7,11,15)]
# colnames(new_dist) <- c("OTU", "V6", "MEAN_manure", "MEAN_day24", "MEAN_day59", "MEAN_day108")
# 



#foo3 <- effluent_otus.melt
#foo <- foo3 %>% mutate(new.day = ifelse(Treatment == "Control", paste(Day, "-C", sep=''), Day))
#foo <- foo %>% mutate(new.day2 = ifelse(Matrix == "Manure", "0-C-M", new.day))
#effluent_plus_controls <- foo


# arg_effluent_melt <- psmelt(arg_effluent2)
# arg_effluent_melt.all <- ddply(arg_effluent_melt, .(ID, Day, OTU, Treatment), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
# arg_effluent_melt.all$presence.absence <- ifelse(arg_effluent_melt.all$MEAN > 0, 1, 0)
# arg_effluent_melt.all2 <- ddply(arg_effluent_melt.all, .(OTU), summarise, TOT_SAMPLES = sum(presence.absence))
# all_effluent_samples <- subset(arg_effluent_melt.all2, TOT_SAMPLES == 4)
# otus.in.all.effluents <- all_effluent_samples$OTU

# # # This is old code for circular heatmaps
# # 
# sw_otus <- prune_taxa(manure_only, arg_phy)
# sw_otus <- subset_samples(sw_otus, Treatment != "Control" | Matrix == "Manure")
# foo <- psmelt(sw_otus)
# foo2 <- ddply(foo, .(OTU, Matrix, Day), summarize, MEAN = mean(Abundance))
# #foo2$Presence.absence <- ifelse(foo2$MEAN > 0, 1, 0)
# order <- foo_temp2[order(foo_temp2$MEAN),]$OTU
# foo2 <- subset(foo2, Day == "24" | Day == "59" | Day == "108")
# foo_24_s <- subset(foo2, Day == "24" & Matrix == "Soil")
# foo_59_s <- subset(foo2, Day == "59" & Matrix == "Soil")
# foo_108_s <- subset(foo2, Day == "108" & Matrix == "Soil")
# foo_24_e <- subset(foo2, Day == "24" & Matrix == "Effluent")
# foo_59_e <- subset(foo2, Day == "59" & Matrix == "Effluent")
# foo_108_e <- subset(foo2, Day == "108" & Matrix == "Effluent")
# merged_foo <- Reduce(function(x, y) merge(x, y, by = "OTU"), list(foo_24_s, foo_59_s, foo_108_s, foo_24_e, foo_59_e, foo_108_e))
# colnames(merged_foo) <- c(paste("V", 1:19, sep=""))

# Old circular plots
# 
library(tidyr) 
library(circlize)
circos.clear()
test = merged_foo
test <- test %>% arrange(-V4)
test$row_id <- seq.int(nrow(test))
factors = test$V2
factors = factor(factors, levels = c("Soil"))
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 5)
circos.initialize(factors, xlim = cbind(c(0, 0), table(factors)))



# circos.track(factors = factors, x = test$row_id, y = test$MEAN_day24, ylim = c(0,1),
#              panel.fun = function(x, y) {
#                sector.index = CELL_META$sector.index
#                circos.points(x, y, pch = 16, cex = 1)
#                #circos.axis(labels.cex = 0.6)
#              })
# circos.track(factors = factors, x = test$row_id, y = test$MEAN_day59, ylim = c(0,1),
#              panel.fun = function(x, y) {
#                sector.index = CELL_META$sector.index
#                circos.points(x, y, pch = 16, cex = 1)
#                #circos.axis(labels.cex = 0.6)
#              })
# circos.track(factors = factors, x = test$row_id, y = test$MEAN_day108, ylim = c(0,1),
#              panel.fun = function(x, y) {
#                sector.index = CELL_META$sector.index
#                circos.points(x, y, pch = 16, cex = 1)
#                #circos.axis(labels.cex = 0.6)
#              })
# #circos.trackHist(factors = factors, x = test$row_id, col = "#999999", border = "#999999")
# circos.clear()
# col_fun = colorRamp2(c(-2, 0, 2), c("green", "black", "red"))
# circos.initializeWithIdeogram(plotType = NULL)
# circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
#   chr = get.cell.meta.data("sector.index")
#   xlim = get.cell.meta.data("xlim")
#   ylim = get.cell.meta.data("ylim")
#   circos.rect(xlim[1], 0, xlim[2], 0.5,
#               col = col_fun(y))
#   #circos.text(mean(xlim), 0.9, chr, cex = 0.5, facing = "clockwise", niceFacing = TRUE)
# }, bg.border = NA)
# 
# 
# test

# Integrall Section -------------------------------------------------------

int_phy <- readRDS("int-phyloseq.RDS")

# RESULT 1 - ARGS IN MANURE AND SOIL BACKGROUNDS
# ARG_PHY HAS 355 GENES AND ARG2 (CARD DATABASE SUBSET) 84 GENES, SO LET'S GO WITH ARG_PHY
abundance_limit = 0 
int_manure <- subset_samples(int_phy, Matrix == "Manure")
int_manure2 <- prune_taxa(taxa_sums(int_manure) > abundance_limit, int_manure)   
int_manure_melt <- psmelt(int_manure2)
int_manure_melt.all <- ddply(int_manure_melt, .(ID, Day, OTU, Treatment), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
int_manure_melt.all$presence.absence <- ifelse(int_manure_melt.all$MEAN > 0, 1, 0)
int_manure_melt.all2 <- ddply(int_manure_melt.all, .(OTU), summarise, TOT_SAMPLES = sum(presence.absence))
all_manure_samples <- subset(int_manure_melt.all2, TOT_SAMPLES == 4) #75 genes
otus.in.all.manures <- all_manure_samples$OTU

manure_otus <- prune_taxa(otus.in.all.manures, int_phy)
manure_otus <- subset_samples(manure_otus, Matrix == "Manure") 
manure_otus <- prune_taxa(taxa_sums(manure_otus) > abundance_limit, manure_otus) 
manure_otus.melt <- psmelt(manure_otus)
manure_dist <- ddply(manure_otus.melt, .(OTU, Annotation), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
rank_order_otus <- manure_dist[order(-manure_dist$MEAN),]$OTU
manure_dist$OTU = factor(manure_dist$OTU, levels = rank_order_otus)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p <-ggplot(manure_dist,aes(OTU, MEAN, color=Annotation))+geom_point(stat = "identity") + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=0, size=10)) + scale_x_discrete(breaks=manure_dist$OTU, labels=manure_dist$Annotation)
p 

# Soil background
int_soil <- subset_samples(int_phy, Matrix == "Soil"& Treatment == "Control")
int_soil2 <- prune_taxa(taxa_sums(int_soil) > abundance_limit, int_soil)  #93 genes in soil
int_soil_melt <- psmelt(int_soil2)
int_soil_melt.all <- ddply(int_soil_melt, .(ID, Day, OTU, Treatment), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
int_soil_melt.all$presence.absence <- ifelse(int_soil_melt.all$MEAN > 0, 1, 0)
int_soil_melt.all2 <- ddply(int_soil_melt.all, .(OTU), summarise, TOT_SAMPLES = sum(presence.absence))
all_soil_samples <- subset(int_soil_melt.all2, TOT_SAMPLES == 4)  
#all_soil_samples <- subset(int_soil_melt.all2, TOT_SAMPLES >= 1) 
otus.in.all.soils <- all_soil_samples$OTU
soil_otus <- prune_taxa(otus.in.all.soils, int_phy)
# 24 genes in all 4 soil samples
soil_otus <- prune_taxa(taxa_sums(soil_otus) > abundance_limit, soil_otus) #24 OTUs
soil_otus.melt <- psmelt(soil_otus)
soil_dist <- ddply(soil_otus.melt, .(OTU, Annotation), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
rank_order_otus <- soil_dist[order(-soil_dist$MEAN),]$OTU
soil_dist$OTU = factor(soil_dist$OTU, levels = rank_order_otus)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p <-ggplot(soil_dist,aes(OTU, MEAN, color=Annotation))+geom_point(stat = "identity") + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=0, size=10)) + scale_x_discrete(breaks=soil_dist$OTU, labels=soil_dist$Annotation)
p

# soil and manure intS together
same_otus <- intersect(taxa_names(soil_otus), taxa_names(manure_otus))  #6 when in only in one sample
manure_otus_in_environment <- prune_taxa(taxa_names(manure_otus), int_phy)
manure_otus_in_environment2 <- subset_samples(manure_otus_in_environment, Treatment == "Control") 
manure_otus_in_environment2 <- subset_samples(manure_otus_in_environment2, Matrix != "Effluent") 
manure_dist <- psmelt(manure_otus_in_environment2)
manure_dist2 <- ddply(manure_dist, .(OTU, Annotation, Matrix), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
manure_dist2$shared <- ifelse(manure_dist2$OTU %in% soil_dist$OTU, 1, 0)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
temp <- subset(manure_dist2, Matrix == "Manure")
rank_order_otus <- temp[order(-temp$MEAN),]$OTU
manure_dist2$OTU = factor(manure_dist2$OTU, levels = rank_order_otus)
p <-ggplot(manure_dist2,aes(OTU, MEAN, color=shared))+geom_point(stat = "identity") + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=0, size=10)) + scale_x_discrete(breaks=manure_dist$OTU, labels=manure_dist$Annotation)
p 

# manure otus, requires all OTUS to be in all 4 manure biological replicates
meta <- sample_data(int_phy)
int_manure <- subset_samples(int_phy, Matrix == "Manure")
#int_manure <- subset_samples(int_phy, Matrix == "Manure")
int_manure2 <- prune_taxa(taxa_sums(int_manure) > 0, int_manure)
#int_manure2 <- prune_taxa(taxa_sums(int_manure) > 0, int_manure)
int_manure_melt <- psmelt(int_manure2)
int_manure_melt.all <- ddply(int_manure_melt, .(ID, Day, OTU, Treatment), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
int_manure_melt.all$presence.absence <- ifelse(int_manure_melt.all$MEAN > 0, 1, 0)
int_manure_melt.all2 <- ddply(int_manure_melt.all, .(OTU), summarise, TOT_SAMPLES = sum(presence.absence))
all_manure_samples <- subset(int_manure_melt.all2, TOT_SAMPLES == 4)
otus.in.all.manures <- all_manure_samples$OTU
manure_otus <- prune_taxa(otus.in.all.manures, int_phy) #49 taxa in CARD, 100 taxa in total
#manure_otus <- prune_taxa(otus.in.all.manures, int_phy2) #49 taxa in CARD, 100 taxa in total

manure_otus <- subset_samples(manure_otus, Matrix == "Manure")
manure_otus <- prune_taxa(taxa_sums(manure_otus) > 0, manure_otus)


# soil otus overlapping with 49 manure OTUs
#int_soil <- subset_samples(int_phy, Matrix == "Soil" & Treatment == "Manure")
#int_soil <- subset_samples(int_phy, Matrix == "Soil" & Treatment == "Manure")
#int_soil2 <- prune_taxa(taxa_sums(int_soil) > 0, int_soil)
#int_soil2 <- prune_taxa(taxa_sums(int_soil) > 0, int_soil)
# int_soil_melt <- psmelt(int_soil2)
# int_soil_melt.all <- ddply(int_soil_melt, .(ID, Day, OTU, Treatment), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
# int_soil_melt.all$presence.absence <- ifelse(int_soil_melt.all$MEAN > 0, 1, 0)
# int_soil_melt.all2 <- ddply(int_soil_melt.all, .(OTU), summarise, TOT_SAMPLES = sum(presence.absence))
# all_soil_samples <- subset(int_soil_melt.all2, TOT_SAMPLES == 4)
# otus.in.all.soils <- all_soil_samples$OTU


soil_otus <- prune_taxa(rownames(otu_table(manure_otus)), int_phy)
#soil_otus <- prune_taxa(rownames(otu_table(manure_otus)), int_phy2)
soil_otus_manured <- subset_samples(soil_otus, Treatment == "Control" & Matrix == "Soil")
soil_otus_manured2 <- prune_taxa(taxa_sums(soil_otus_manured) == 0, soil_otus_manured)
soil_manure_otus <- taxa_names(soil_otus_manured2)
effluent_otus_manured <- subset_samples(soil_otus, Treatment == "Control" & Matrix == "Effluent")
effluent_otus_manured2 <- prune_taxa(taxa_sums(effluent_otus_manured) == 0, effluent_otus_manured)
effluent_manure_otus <- taxa_names(effluent_otus_manured2)
manure_only <- intersect(soil_manure_otus, effluent_manure_otus)
sw_otus <- prune_taxa(manure_only, int_phy)
sw_otus.melt <- psmelt(sw_otus)
foo <- sw_otus.melt

# Heatmap
foo$Matrix <- factor(foo$Matrix, levels = c("Manure", "Soil", "Effluent"))
temp2 <- subset(foo, OTU == "gnl|AF458080|idbc3968")
sample_order <- temp2[order(temp2$Matrix),]$Sample
foo$Sample <- factor(foo$Sample, levels = sample_order)

foo_temp <- subset(foo, Treatment == "Control" & Matrix == "Manure")
foo_temp2 <- ddply(foo_temp, .(OTU, Annotation), summarize, MEAN = mean(Abundance))
order <- foo_temp2[order(foo_temp2$MEAN),]$OTU
foo$Day = factor(foo$Day, levels=c("0", "10","24", "38", "59","80","108"))
foo$OTU = factor(foo$OTU, levels = order)
p <-ggplot(foo,aes(Sample,OTU,fill=Abundance))+
  geom_tile(color = "white", size = 0.25) 
p <-p + facet_grid(Treatment ~ Day,  scale = "free", space = "free")+scale_y_discrete(breaks=foo$OTU, labels=foo$Annotation) +
  scale_x_discrete(breaks = foo$Sample, labels = foo$Matrix)
p +scale_fill_gradient(low = "white", high = "red", trans = "log") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Heatmap Without controls
sw_otus <- prune_taxa(manure_only, int_phy)
sw_otus <- subset_samples(sw_otus, Treatment != "Control" | Matrix == "Manure")
foo <- psmelt(sw_otus)
foo$Abundance <- ifelse(foo$Abundance > abundance_limit, foo$Abundance, 0)
foo$Matrix <- factor(foo$Matrix, levels = c("Manure", "Soil", "Effluent"))
temp2 <- subset(foo, OTU == "gnl|AF458080|idbc3968")
sample_order <- temp2[order(temp2$Matrix),]$Sample
foo$Sample <- factor(foo$Sample, levels = sample_order)
foo_temp <- subset(foo, Treatment == "Control" & Matrix == "Manure")
foo_temp2 <- ddply(foo_temp, .(OTU, Annotation), summarize, MEAN = mean(Abundance))
order <- foo_temp2[order(foo_temp2$MEAN),]$OTU
foo$Day = factor(foo$Day, levels=c("0", "10","24", "38", "59","80","108"))
foo$OTU = factor(foo$OTU, levels = order)
p <-ggplot(foo,aes(Sample,OTU,fill=Abundance))+
  geom_tile(color = "white", size = 0.25) 
p <-p + facet_grid( ~ Day,  scale = "free", space = "free")+scale_y_discrete(breaks=foo$OTU, labels=foo$Annotation) +
  scale_x_discrete(breaks = foo$Sample, labels = foo$Matrix)
p +scale_fill_gradient(low = "white", high = "red", trans = "log") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Heatmap - abundance curves
ordering <- foo[order(-foo$Abundance),]
ordering$rank_order <- seq.int(nrow(ordering))
p <-ggplot(ordering,aes(rank_order, log(Abundance)))+geom_point(stat = "identity")# + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=0, size=10)) 
p

# Clustering
sw_otus
distance_sw <- distance(sw_otus, method="euclidean")
sw.hclust <- hclust(distance_sw, method = 'complete')
plot(sw.hclust, hang= -1)
hcd = as.dendrogram(sw.hclust)
plot(cut(hcd, h = 0.4)$upper, main = "Upper tree of cut at h=0.4")

library(dendextend)
cbPalette <- c("red", "green", "blue")
colorCode <- c(Manure=cbPalette[1], Effluent = cbPalette[2], Soil = cbPalette[3])
labels_colors(hcd) <- colorCode[sample_data(sw_otus)$Matrix][order.dendrogram(hcd)]
labels(hcd) <- sample_data(sw_otus)$Day[order.dendrogram(hcd)]
plot(hcd)

# Making Table
int_manure <- subset_samples(int_phy, Matrix == "Manure")
int_manure2 <- prune_taxa(taxa_sums(int_manure) > abundance_limit, int_manure)  #163 genes in manure
int_manure_melt <- psmelt(int_manure2)
int_manure_melt.all <- ddply(int_manure_melt, .(ID, Day, OTU, Treatment), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
int_manure_melt.all$presence.absence <- ifelse(int_manure_melt.all$MEAN > 0, 1, 0)
int_manure_melt.all2 <- ddply(int_manure_melt.all, .(OTU), summarise, TOT_SAMPLES = sum(presence.absence))
all_manure_samples <- subset(int_manure_melt.all2, TOT_SAMPLES == 4) #100 genes
all_manure_samples <- subset(int_manure_melt.all2, TOT_SAMPLES >= 1) #163 genes
otus.in.all.manures <- all_manure_samples$OTU

manure_otus <- prune_taxa(otus.in.all.manures, int_phy)
manure_otus_day24 <- subset_samples(manure_otus, Day == "24" & Matrix != "Manure")
manure_otus_day24 <- subset_samples(manure_otus_day24, Treatment == "Manure")
soil_24 <- subset_samples(manure_otus_day24, Matrix == "Soil")
soil_24 <- prune_taxa(taxa_sums(soil_24) > abundance_limit, soil_24) #130 OTUs 
soil_24 
eff_24 <- subset_samples(manure_otus_day24, Matrix == "Effluent")
eff_24 <- prune_taxa(taxa_sums(eff_24) > abundance_limit, eff_24) #88
eff_24
same_otus <- intersect(taxa_names(soil_24), taxa_names(eff_24)) #87
length(same_otus)

manure_otus <- prune_taxa(otus.in.all.manures, int_phy)
manure_otus_day59 <- subset_samples(manure_otus, Day == "59" & Matrix != "Manure")
manure_otus_day59 <- subset_samples(manure_otus_day59, Treatment == "Manure")
soil_59 <- subset_samples(manure_otus_day59, Matrix == "Soil")
soil_59 <- prune_taxa(taxa_sums(soil_59) > abundance_limit, soil_59) #131
soil_59 
eff_59 <- subset_samples(manure_otus_day59, Matrix == "Effluent")
eff_59 <- prune_taxa(taxa_sums(eff_59) > abundance_limit, eff_59)  #92
eff_59
same_otus <- intersect(taxa_names(soil_59), taxa_names(eff_59)) #89
length(same_otus) 

manure_otus <- prune_taxa(otus.in.all.manures, int_phy)
manure_otus_day108 <- subset_samples(manure_otus, Day == "108" & Matrix != "Manure")
manure_otus_day108 <- subset_samples(manure_otus_day108, Treatment == "Manure")
soil_108 <- subset_samples(manure_otus_day108, Matrix == "Soil")
soil_108 <- prune_taxa(taxa_sums(soil_108) > abundance_limit, soil_108) #126
soil_108 
eff_108 <- subset_samples(manure_otus_day108, Matrix == "Effluent")
eff_108 <- prune_taxa(taxa_sums(eff_108) > abundance_limit, eff_108) #119
eff_108
same_otus <- intersect(taxa_names(soil_108), taxa_names(eff_108)) #113 shared; 
length(same_otus)

#checking differences
diff1 <- setdiff(taxa_names(soil_108), taxa_names(eff_108))
diff2 <- setdiff(taxa_names(eff_108), taxa_names(soil_108))
diff <- prune_taxa(diff2, int_phy)
diff.melt <- psmelt(diff)
diff_dist <- ddply(diff.melt, .(OTU, Matrix, Treatment, Day), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
diff_dist <- subset(diff_dist, Matrix != "Manure")
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p <-ggplot(diff_dist,aes(OTU, MEAN, color=Treatment))+geom_point(stat = "identity") + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=0, size=10)) 
p + facet_grid(Day~Matrix)

#abundance curves
soil_manure <- subset_samples(int_phy, Matrix == "Soil" & Treatment == "Manure")
soil_manure.melt <- psmelt(soil_manure)
soil_manure <- ddply(soil_manure.melt, .(OTU, Day), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
soil_manure <- subset(soil_manure, MEAN > 0)
ordering <- soil_manure[order(-soil_manure$MEAN),]
ordering$rank_order <- seq.int(nrow(ordering))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p <-ggplot(ordering,aes(rank_order, log(MEAN)))+geom_point(stat = "identity")# + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=0, size=10)) 
p


effluent_manure <- subset_samples(int_phy, Matrix == "Effluent" & Treatment == "Manure")
effluent_manure.melt <- psmelt(effluent_manure)
effluent_manure <- ddply(effluent_manure.melt, .(OTU, Day), summarise, MEAN = mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
effluent_manure <- subset(effluent_manure, MEAN > 0)
ordering <- effluent_manure[order(-effluent_manure$MEAN),]
ordering$rank_order <- seq.int(nrow(ordering))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p <-ggplot(ordering,aes(rank_order, log10(MEAN)))+geom_point(stat = "identity")# + geom_errorbar(limits, width=0) + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=0, size=10)) 
p

# Occupancy Curve

foo <- sw_otus.melt
foo$Presence.absence <- ifelse(foo$Abundance > abundance_limit, 1, 0)
foo_matrix <- subset(foo, Treatment == "Manure" & Matrix == "Soil")
foo2 <- ddply(foo_matrix, .(OTU, V2), summarise, SUM = sum(Presence.absence))
ordering <- foo2[order(foo2$SUM),]
ordering$rank <- seq.int(nrow(ordering))
ordering$percent <- ordering$SUM/length(unique(foo_matrix$Sample))
p <-ggplot(ordering,aes(rank, percent))+geom_point(stat = "identity")+scale_x_discrete(breaks = ordering$rank, labels = ordering$rank)
p
