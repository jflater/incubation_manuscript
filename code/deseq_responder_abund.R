library(phyloseq)
library(vegan)
library(tidyverse)
library(gplots)
library(DESeq2)
inc.raw.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

inc.physeq <- subset_samples(inc.raw.physeq, day %in% c("0",
                                                        "7",
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

inc.physeq.data <- data.frame(sample_data(inc.physeq))
inc.physeq.data$response.group[inc.physeq.data$day == "0"] <- "baseline" 
inc.physeq.data$response.group[inc.physeq.data$day %in% c("7", "14", "21")] <- "early" 
inc.physeq.data$response.group[inc.physeq.data$day %in% c("35", "49", "97")] <- "late" 
inc.physeq.data <- inc.physeq.data %>%
  mutate(Treatment_Response = paste(treatment, response.group, sep = '_'))
rownames(inc.physeq.data) <- data$i_id
sample_data(inc.physeq) <- inc.physeq.data

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
alf.physeq <- subset_samples(inc.physeq, treatment %in% c("Alfalfa")) %>%
  filter_taxa(function(x) sum(x) >= 2, T)

log.plot.early.alf <- alf.physeq %>%
  phyloseq_to_deseq2( ~ response.group) %>%
  DESeq(test = "Wald", fitType = "local") %>%
  who_diff_day("early", "baseline", alf.physeq) %>%
  log_plot("Alfalfa OTUS in early group that are significantly changing compared to day 0")
log.plot.early.alf

log.plot.late.alf <- alf.physeq %>%
  phyloseq_to_deseq2( ~ response.group) %>%
  DESeq(test = "Wald", fitType = "local") %>%
  who_diff_day("late", "early", alf.physeq) %>%
  log_plot("Alfalfa OTUS in late group that are significantly changing compared to early group")
log.plot.late.alf
###
comp.physeq <- subset_samples(inc.physeq, treatment %in% c("Compost")) %>%
  filter_taxa(function(x) sum(x) >= 2, T)

log.plot.early.comp <- comp.physeq %>%
  phyloseq_to_deseq2( ~ response.group) %>%
  DESeq(test = "Wald", fitType = "local") %>%
  who_diff_day("early", "baseline", comp.physeq) %>%
  log_plot("Compost OTUS in early group that are significantly changing compared to day 0")
log.plot.early.comp

log.plot.late.comp <- comp.physeq %>%
  phyloseq_to_deseq2( ~ response.group) %>%
  DESeq(test = "Wald", fitType = "local") %>%
  who_diff_day("late", "early", comp.physeq) %>%
  log_plot("Compost OTUS in late group that are significantly changing compared to early group")
log.plot.late.comp
###
mix.physeq <- subset_samples(inc.physeq, treatment %in% c("Mix")) %>%
  filter_taxa(function(x) sum(x) >= 2, T)

log.plot.early.mix <- mix.physeq %>%
  phyloseq_to_deseq2( ~ response.group) %>%
  DESeq(test = "Wald", fitType = "local") %>%
  who_diff_day("early", "baseline", mix.physeq) %>%
  log_plot("Mix OTUS in early group that are significantly changing compared to day 0")
log.plot.early.mix

log.plot.late.mix <- mix.physeq %>%
  phyloseq_to_deseq2( ~ response.group) %>%
  DESeq(test = "Wald", fitType = "local") %>%
  who_diff_day("late", "early", mix.physeq) %>%
  log_plot("Mix OTUS in late group that are significantly changing compared to early group")
log.plot.late.mix
###
ref.physeq <- subset_samples(inc.physeq, treatment %in% c("Reference")) %>%
  filter_taxa(function(x) sum(x) >= 2, T)

log.plot.early.ref <- ref.physeq %>%
  phyloseq_to_deseq2( ~ response.group) %>%
  DESeq(test = "Wald", fitType = "local") %>%
  who_diff_day("early", "baseline", ref.physeq) %>%
  log_plot("Reference OTUS in early group that are significantly changing compared to day 0")
log.plot.early.ref

log.plot.late.ref <- ref.physeq %>%
  phyloseq_to_deseq2( ~ response.group) %>%
  DESeq(test = "Wald", fitType = "local") %>%
  who_diff_day("late", "early", ref.physeq) %>%
  log_plot("Reference OTUS in late group that are significantly changing compared to early group")
log.plot.late.ref

early.alf.otus <- log.plot.early.alf$data %>%
  rownames_to_column() %>%
  select(OTU, Phylum, Genus, log2FoldChange) %>%
  filter(log2FoldChange >= 2) %>%
  arrange(desc(log2FoldChange))
early.alf.otus 
saveRDS(early.alf.otus, "data/early.alf.otus.rds")
###
early.comp.otus <- log.plot.early.comp$data %>%
  rownames_to_column() %>%
  select(OTU, Phylum, Genus, log2FoldChange) %>%
  filter(log2FoldChange >= 2)
early.comp.otus
saveRDS(early.comp.otus, "data/early.comp.otus.rds")
###
early.mix.otus <- log.plot.early.mix$data %>%
  rownames_to_column() %>%
  select(OTU, Phylum, Genus, log2FoldChange) %>%
  filter(log2FoldChange >= 2)
early.mix.otus
saveRDS(early.mix.otus, "data/early.mix.otus.rds")
###
early.ref.otus <- log.plot.early.ref$data %>%
  rownames_to_column() %>%
  select(OTU, Phylum, Genus, log2FoldChange) %>%
  filter(log2FoldChange >= 2)
early.ref.otus
saveRDS(early.ref.otus, "data/early.ref.otus.rds")

venn <- venn(list(Alfalfa = early.alf.otus$OTU, Reference = early.ref.otus$OTU,
                  Mix = early.mix.otus$OTU, Compost = early.comp.otus$OTU))

late.alf.otus <- log.plot.late.alf$data %>%
  rownames_to_column() %>%
  select(OTU, Phylum, Genus, log2FoldChange) %>%
  filter(log2FoldChange >= 2)
late.alf.otus 
saveRDS(late.alf.otus, "data/late.alf.otus.rds")

###
late.comp.otus <- log.plot.late.comp$data %>%
  rownames_to_column() %>%
  select(OTU, Phylum, Genus, log2FoldChange) %>%
  filter(log2FoldChange >= 2)
late.comp.otus
saveRDS(late.comp.otus, "data/late.comp.otus.rds")

###
late.mix.otus <- log.plot.late.mix$data %>%
  rownames_to_column() %>%
  select(OTU, Phylum, Genus, log2FoldChange) %>%
  filter(log2FoldChange >= 2)
late.mix.otus
saveRDS(late.mix.otus, "data/late.mix.otus.rds")

###
late.ref.otus <- log.plot.late.ref$data %>%
  rownames_to_column() %>%
  select(OTU, Phylum, Genus, log2FoldChange) %>%
  filter(log2FoldChange >= 2)
late.ref.otus
saveRDS(late.ref.otus, "data/late.ref.otus.rds")

venn2 <- venn(list(Alfalfa = late.alf.otus$OTU,
                   Reference = late.ref.otus$OTU,
                   Mix = late.mix.otus$OTU,
                   Compost = late.comp.otus$OTU))

alf.early.otus <- attr(venn, "intersections")$Alfalfa

rare6k.physeq <- rarefy_even_depth(inc.physeq, sample.size = 6000,
                                   rngseed = 15879966) %>%
  filter_taxa(function(x) sum(x) >= 2, T) 

# rare6k.physeq.data <- data.frame(sample_data(rare6k.physeq))
# rare6k.physeq.data$response.group[rare6k.physeq.data$day == "0"] <- "baseline" 
# rare6k.physeq.data$response.group[rare6k.physeq.data$day %in% c("7", "14", "21")] <- "early" 
# rare6k.physeq.data$response.group[rare6k.physeq.data$day %in% c("35", "49", "97")] <- "late" 
# sample_data(rare6k.physeq) <- rare6k.physeq.data

RelativeAbundanceDf <- function(physeq) {
  physeq %>% tax_glom(taxrank = "Phylum") %>% transform_sample_counts(function(x) {
    x/sum(x)
  }) %>% psmelt() %>% filter(Abundance > 0.02) %>% arrange(Phylum)
}

PlotRelativeAbundance <- function(df) {
  ggplot(df, aes(x = as.factor(Treatment_Response), y = Abundance, fill = Phylum)) + 
    geom_bar(stat = "identity") +
    theme(axis.title.x = element_blank()) + 
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) + 
    ylab("Phylogenetic distribution of OTUs with abundance > 2%") 
}


early.alf <- subset_samples(rare6k.physeq, treatment %in% c("Alfalfa")) %>%
  subset_samples(day %in% c("0", "7", "14", "21")) 
only.early.alf <- subset_taxa(early.alf, rownames(tax_table(early.alf)) %in% alf.early.otus)

comp.early.otus <- attr(venn, "intersections")$Compost
early.comp <- subset_samples(rare6k.physeq, treatment %in% c("Compost")) %>%
  subset_samples(day %in% c("0", "7", "14", "21")) 
only.early.comp <- subset_taxa(early.comp, rownames(tax_table(early.comp)) %in% comp.early.otus)

lotu <- paste(c(comp.early.otus, alf.early.otus))
compalf <- subset_samples(rare6k.physeq, treatment %in% c("Compost", "Alfalfa")) %>%
  subset_samples(day %in% c("7", "14", "21", "35", "49", "97")) 

only.early.compalf <- subset_taxa(compalf, rownames(tax_table(compalf)) %in% lotu) 
# plot abundance
plot_bar(only.early.compalf, x = "response.group", fill = "Phylum", facet_grid = ~treatment) +
  geom_bar(aes(color= Phylum, fill= Phylum), stat = "identity", position = "stack")

alf.late.otus <- attr(venn2, "intersections")$Alfalfa
late.alf <- subset_samples(rare6k.physeq, treatment %in% c("Alfalfa")) %>%
  subset_samples(day %in% c("7", "14", "21", "35", "49", "97")) 
only.late.alf <- subset_taxa(late.alf, rownames(tax_table(late.alf)) %in% alf.late.otus)

comp.late.otus <- attr(venn2, "intersections")$Compost
late.comp <- subset_samples(rare6k.physeq, treatment %in% c("Compost")) %>%
  subset_samples(day %in% c("7", "14", "21", "35", "49", "97")) 
only.late.comp <- subset_taxa(late.comp, rownames(tax_table(late.comp)) %in% comp.late.otus)

lotu2 <- paste(c(comp.late.otus))
compalf.late <- subset_samples(rare6k.physeq, treatment %in% c("Compost", "Alfalfa", "Reference", "Mix")) %>%
  subset_samples(day %in% c("7", "14", "21", "35", "49", "97")) 

only.late.compalf <- subset_taxa(compalf.late, rownames(tax_table(compalf.late)) %in% lotu2)
# plot abundance
plot <- plot_bar(only.late.compalf, x = "response.group", fill = "Phylum", facet_grid = ~treatment) +
  geom_bar(aes(color= Phylum, fill= Phylum), stat = "identity", position = "stack")
tiff("Figures/late_comp_abund.tif",height=4,width=6,units='in',res=400)
plot
dev.off()
tiff("Figures/earlyalflog.tif",height=8,width=6,units='in',res=400)
log.plot.early.alf
dev.off()
tiff("Figures/latealflog.tif",height=8,width=6,units='in',res=400)
log.plot.late.alf
dev.off()