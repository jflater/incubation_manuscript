library(tidyverse)
library(phyloseq)
library(DESeq2)

# Load phyloseq object as inc.physeq
inc.physeq <- readRDS("data/RDS/not.rare.nounclass")

# LFC calculation and visualization function
plot_deseq <- function(phyloseq_object){
  # Subset phyloseq object 
  physeq <- phyloseq_object %>%
    filter_taxa(function(x) sum(x) >= 3, T)
  # Convert to desey and perfrom differential expression analysis
  diagdds = phyloseq_to_deseq2(physeq, ~ treatment)
  diagdds = DESeq(diagdds, test="Wald", fitType="local")
  res = results(diagdds, cooksCutoff = FALSE)
  alpha = 0.01
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  ####
  #theme_set(theme_bw())
  #scale_fill_discrete <- function(palname = "Set1", ...) {
  #  scale_fill_brewer(palette = palname, ...)
  #}
  # Phylum order
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  # Genus order
  x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
  plot <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
  return(plot)
}

Alfalfa_early = c("Alfalfa_early", "Reference_early")
Alfalfa_early <- subset_samples(inc.physeq, Treatment_Response %in% Alfalfa_early)
Alfalfa_early <- plot_deseq(Alfalfa_early)
#Alfalfa_early
Alfalfa_early_data <- Alfalfa_early$data %>%
  rownames_to_column(var = "OTU") %>%
  filter(log2FoldChange >= 0) 
Alfalfa_early_data$trt <- c("Alfalfa_early")

Alfalfa_late = c("Alfalfa_late", "Reference_late")
Alfalfa_late <- subset_samples(inc.physeq, Treatment_Response %in% Alfalfa_late)
Alfalfa_late <- plot_deseq(Alfalfa_late)
#Alfalfa_late
Alfalfa_late_data <- Alfalfa_late$data %>%
  rownames_to_column(var = "OTU") %>%
  filter(log2FoldChange >= 0) 
Alfalfa_late_data$trt <- c("Alfalfa_late")

otustokeep <- intersect(Alfalfa_early_data$OTU, Alfalfa_late_data$OTU)
test <- rbind.data.frame(Alfalfa_early_data, Alfalfa_late_data) %>%
  filter(OTU %in% otustokeep)

ggplot(test, aes(x=Genus, y=log2FoldChange, color=Phylum, shape = trt)) + geom_point(size=2) + 
    coord_flip() 

all_alf <- Alfalfa_early_data %>%
  select(OTU, Phylum, Genus, Alfalfa_early_log2FoldChange = log2FoldChange) %>%
  filter(OTU %in% otustokeep) %>%
  full_join(Alfalfa_late_data) %>%
  select(OTU, Phylum, Genus, Alfalfa_early_log2FoldChange, Alfalfa_late_log2FoldChange = log2FoldChange)

ggplot(all_alf, aes(x = Alfalfa_early_log2FoldChange, y = Alfalfa_late_log2FoldChange, color = Phylum, label = Genus)) + geom_point(size = 1) +
  coord_flip() +
  geom_text(check_overlap = T, hjust = 0, nudge_x = 0.05) 

alf_OTUS <- Alfalfa_early_data %>%
  select(OTU, Phylum, Genus, Alfalfa_early_log2FoldChange = log2FoldChange) %>%
  full_join(Alfalfa_late_data) %>%
  select(OTU, Phylum, Genus, Alfalfa_early_log2FoldChange, Alfalfa_late_log2FoldChange = log2FoldChange)

head(alf_OTUS)

saveRDS(alf_OTUS, file = "data/RDS/LFC_alf_OTUs_June19.RDS")
