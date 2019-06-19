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

Compost_early = c("Compost_early", "Reference_early")
Compost_early <- subset_samples(inc.physeq, Treatment_Response %in% Compost_early)
Compost_early <- plot_deseq(Compost_early)
#Compost_early
Compost_early_data <- Compost_early$data %>%
  rownames_to_column(var = "OTU") %>%
  filter(log2FoldChange >= 0) 
Compost_early_data$trt <- c("Compost_early")

Compost_late = c("Compost_late", "Reference_late")
Compost_late <- subset_samples(inc.physeq, Treatment_Response %in% Compost_late)
Compost_late <- plot_deseq(Compost_late)
#Compost_late
Compost_late_data <- Compost_late$data %>%
  rownames_to_column(var = "OTU") %>%
  filter(log2FoldChange >= 0) 
Compost_late_data$trt <- c("Compost_late")

otustokeep <- intersect(Compost_early_data$OTU, Compost_late_data$OTU)
test <- rbind.data.frame(Compost_early_data, Compost_late_data) %>%
  filter(OTU %in% otustokeep)

ggplot(test, aes(x=Genus, y=log2FoldChange, color=Phylum, shape = trt)) + geom_point(size=2) + 
  coord_flip() 

all_comp <- Compost_early_data %>%
  select(OTU, Phylum, Genus, Compost_early_log2FoldChange = log2FoldChange) %>%
  filter(OTU %in% otustokeep) %>%
  full_join(Compost_late_data) %>%
  select(OTU, Phylum, Genus, Compost_early_log2FoldChange, Compost_late_log2FoldChange = log2FoldChange)

ggplot(all_comp, aes(x = Compost_early_log2FoldChange, y = Compost_late_log2FoldChange, color = Phylum, label = Genus)) + geom_point(size = 1) +
  coord_flip() +
  geom_text(check_overlap = T, hjust = 0, nudge_x = 0.05) 

# Save list of OTUs that have LFC > 0 from early and late data
comp_OTUS <- Compost_early_data %>%
  select(OTU, Phylum, Genus, Compost_early_log2FoldChange = log2FoldChange) %>%
  full_join(Compost_late_data) %>%
  select(OTU, Phylum, Genus, Compost_early_log2FoldChange, Compost_late_log2FoldChange = log2FoldChange)

head(comp_OTUS)

saveRDS(comp_OTUS, file = "data/RDS/LFC_comp_OTUs_June19.RDS")

