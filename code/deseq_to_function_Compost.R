library(phyloseq)
library(vegan)
library(tidyverse)
library(gplots)
library(DESeq2)
library(ggtree)


inc.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

tree <- read.tree("data/tree.nwk")

inc.physeq <- merge_phyloseq(inc.physeq, tree)

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

inc.physeq.data <- data.frame(sample_data(inc.physeq))
inc.physeq.data$response.group[inc.physeq.data$day == "0"] <-
  "baseline"
inc.physeq.data$response.group[inc.physeq.data$day %in% c("7", "14", "21")] <-
  "early"
inc.physeq.data$response.group[inc.physeq.data$day %in% c("35", "49", "97")] <-
  "late"

inc.physeq.data <- inc.physeq.data %>% 
  mutate(Treatment_Response = paste(treatment, response.group, sep = '_'))

rownames(inc.physeq.data) <- data$i_id
sample_data(inc.physeq) <- inc.physeq.data

no.unclass <- subset_taxa(inc.physeq, !Phylum=="Bacteria_unclassified")
no.unclass <- subset_taxa(no.unclass, !Genus=="Gp6_unclassified")
inc.physeq <- no.unclass

who_diff_day <- function(DDS, choice1, choice2, phy.object){
  res = results(DDS, contrast = c("response.group", choice1, choice2), cooksCutoff = FALSE)
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

# Use inc.physeq, not rarefied as DESeq does this
# Use treatment and response group columns to capture early alfalfa and baseline alfalfa and reference to compare
#######
comp.physeq <- subset_samples(inc.physeq, Treatment_Response %in% c("Compost_early", "Compost_baseline", "Reference_baseline")) %>%
  filter_taxa(function(x) sum(x) >= 3, T) %>%
  tax_glom(taxrank = "Genus")

log.plot.early.comp <- comp.physeq %>%
  phyloseq_to_deseq2( ~ response.group) %>%
  DESeq(test = "Wald", fitType = "local") %>%
  who_diff_day("early", "baseline", comp.physeq) %>%
  log_plot("Compost OTUS in early group that are significantly changing compared to day 0")

log.plot.early.comp

png("Figures/log.plot.early.comp",height=5,width=6,units='in',res=300)
plot(plot)
dev.off()

comp.late.physeq <- subset_samples(inc.physeq, Treatment_Response %in% c("Compost_late", "Compost_early", "Reference_early")) %>%
  filter_taxa(function(x) sum(x) >= 3, T) %>%
  tax_glom(taxrank = "Genus")

log.plot.late.comp <- comp.late.physeq %>%
  phyloseq_to_deseq2( ~ response.group) %>%
  DESeq(test = "Wald", fitType = "local") %>%
  who_diff_day("late", "early", comp.late.physeq) %>%
  log_plot("Compost OTUS in late group that are significantly changing compared to early group")

log.plot.late.comp

early.comp.otus <- log.plot.early.comp$data %>%
  rownames_to_column() %>%
  select(OTU, Phylum, Genus, log2FoldChange) %>%
  filter(log2FoldChange >= 2) %>%
  arrange(desc(log2FoldChange))

early.comp.otus 

saveRDS(early.comp.otus, "data/early.comp.otus.rds")

late.comp.otus <- log.plot.late.comp$data %>%
  rownames_to_column() %>%
  select(OTU, Phylum, Genus, log2FoldChange) %>%
  filter(log2FoldChange >= 2)

late.comp.otus 

saveRDS(late.comp.otus, "data/late.comp.otus.rds")