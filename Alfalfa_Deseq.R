library(tidyverse)
library(phyloseq)
library(DESeq2)
library(naniar)
library(ggrepel)
library(viridis)
library(knitr)
library(kableExtra)
########
####### MAKE THE OTHER AMENDMENT SCRIPTS LOOK LIKE THIS!!! NEED TO EDIT THEM AS OF JUNE 24 2019
######
# Load phyloseq object as inc.physeq, not rarefied but unclassified OTUs removed if at phylum level
inc.physeq <- readRDS("data/RDS/not.rare.nounclass")

# LFC calculation function
who_diff_day <- function(DDS, choice1, choice2, phy.object){
  res = results(DDS, contrast = c("Treatment_Response", choice1, choice2), cooksCutoff = FALSE)
  #plotCounts(AlfalfaDDS, gene="OTU_311", intgroup="day")
  #Use above line to check if an OTU is increasing or decreasing depending on order of contrast
  alpha = 0.01
  #alpha = 0.1
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy.object)[rownames(sigtab), ], "matrix"))
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  # Phylum order
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  # Genus order
  x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
  #ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=phylum)) + geom_point(size=2) + 
  #  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=1.0)) +
  #  ggtitle("Day 0 to Day 7")
  return(sigtab)
}
# function plot log2FoldChange 
log_plot <- function(sigtab,t1){
  sigtab <- sigtab %>%
    rownames_to_column(var = "OTU") %>%
    filter(log2FoldChange >= 2) 
  
  ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + 
    coord_flip() +
    ggtitle(t1)
} 
#####             
# Alfalfa Resopnders  

## Early   
alf.early <- subset_samples(inc.physeq, Treatment_Response %in% c("Alfalfa_early", "Reference_early")) %>%
  filter_taxa(function(x) sum(x) >= 3, T) 
# Be very careful of the design formula in the who_diff_day() function
# This function also selects only LFC >= 2 and alpha 0.01 for significant and increasing otus to be returned
log.plot.early.alf <- alf.early %>%
  phyloseq_to_deseq2( ~ Treatment_Response) %>%
  DESeq(test = "Wald", fitType = "local") %>%
  who_diff_day("Alfalfa_early", "Reference_early", alf.early) %>%
  log_plot("Alfalfa OTUS in early group that are significantly changing compared to reference early")
# Save a data frame of these results
log.plot.early.alf.data <- log.plot.early.alf$data %>%
  mutate(trt = c("Alfalfa_early"))
log.plot.early.alf + scale_colour_viridis_d(option = "plasma") +
  theme_dark()

## Late 
alf.late <- subset_samples(inc.physeq, Treatment_Response %in% c("Alfalfa_late", "Reference_late")) %>%
  filter_taxa(function(x) sum(x) >= 3, T)
sample_data(alf.late)
# Make deseq and plot as above
log.plot.late.alf <- alf.late %>%
  phyloseq_to_deseq2( ~ Treatment_Response) %>%
  DESeq(test = "Wald", fitType = "local") %>%
  who_diff_day("Alfalfa_late", "Reference_late", alf.late) %>%
  log_plot("Alfalfa OTUS in late group that are significantly changing compared to reference late")
# Save a data frame of these results
log.plot.late.alf.data <- log.plot.late.alf$data %>%
  mutate(trt = c("Alfalfa_late"))
log.plot.late.alf + scale_colour_viridis_d(option = "plasma") +
  theme_dark()

## How many OTUs for the early and late groups resopnding with LFC > 2?
length(log.plot.early.alf.data$OTU)
length(log.plot.late.alf.data$OTU)

otustokeep <- intersect(log.plot.early.alf.data$OTU, log.plot.late.alf.data$OTU)

all_alf1 <- rbind.data.frame(log.plot.early.alf.data, log.plot.late.alf.data)

ggplot(all_alf1, aes(x=Genus, y=log2FoldChange, color=Phylum, shape = trt)) + geom_point(size=2) + 
    coord_flip() 

common_alf <- log.plot.early.alf.data %>%
  select(OTU, Phylum, Genus, Alfalfa_early_log2FoldChange = log2FoldChange) %>%
  full_join(log.plot.late.alf.data) %>%
  select(OTU, Phylum, Genus, Alfalfa_early_log2FoldChange, Alfalfa_late_log2FoldChange = log2FoldChange) %>%
  filter(OTU %in% otustokeep)

early_alf_OTUS <- log.plot.early.alf.data %>%
  select(OTU, Phylum, Class, Order, Family, Genus, Alfalfa_early_log2FoldChange = log2FoldChange) 

late_alf_OTUS <- log.plot.late.alf.data %>%
  select(OTU, Phylum, Class, Order, Family, Genus, Alfalfa_late_log2FoldChange = log2FoldChange)

all_alf <- full_join(early_alf_OTUS, late_alf_OTUS)

head(all_alf)

#ggplot(common_alf, aes(x = Alfalfa_early_log2FoldChange, y = Alfalfa_late_log2FoldChange, color = Phylum, label = Genus)) + geom_point(size = 1) +
#  coord_flip() +
#  geom_text(check_overlap = T, hjust = 0, nudge_x = 0.05) 


p <- ggplot(all_alf,
       aes(x = Alfalfa_early_log2FoldChange, 
           y = Alfalfa_late_log2FoldChange,
           color = Phylum,
           label = OTU)) +
  geom_miss_point() +
  geom_text_repel(aes(label=ifelse(Alfalfa_early_log2FoldChange>3 & Alfalfa_late_log2FoldChange>3,as.character(OTU),'')),hjust=0,vjust=0) 
p
p + scale_colour_viridis_d(option = "plasma") +
  theme_dark()
saveRDS(all_alf, file = "data/RDS/LFC_alf_OTUs_June19.RDS")

OTUs <- all_alf %>%
  filter(Alfalfa_early_log2FoldChange>3 & Alfalfa_late_log2FoldChange>3)

kable(OTUs) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)