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

Mix_early = c("Mix_early", "Reference_early")
Mix_early <- subset_samples(inc.physeq, Treatment_Response %in% Mix_early)
Mix_early <- plot_deseq(Mix_early)
#Mix_early
Mix_early_data <- Mix_early$data %>%
  rownames_to_column(var = "OTU") %>%
  filter(log2FoldChange >= 0) 
Mix_early_data$trt <- c("Mix_early")

Mix_late = c("Mix_late", "Reference_late")
Mix_late <- subset_samples(inc.physeq, Treatment_Response %in% Mix_late)
Mix_late <- plot_deseq(Mix_late)
#Mix_late
Mix_late_data <- Mix_late$data %>%
  rownames_to_column(var = "OTU") %>%
  filter(log2FoldChange >= 0) 
Mix_late_data$trt <- c("Mix_late")

otustokeep <- intersect(Mix_early_data$OTU, Mix_late_data$OTU)
test <- rbind.data.frame(Mix_early_data, Mix_late_data) %>%
  filter(OTU %in% otustokeep)

ggplot(test, aes(x=Genus, y=log2FoldChange, color=Phylum, shape = trt)) + geom_point(size=2) + 
  coord_flip() 

all_mix <- Mix_early_data %>%
  select(OTU, Phylum, Genus, Mix_early_log2FoldChange = log2FoldChange) %>%
  filter(OTU %in% otustokeep) %>%
  full_join(Mix_late_data) %>%
  select(OTU, Phylum, Genus, Mix_early_log2FoldChange, Mix_late_log2FoldChange = log2FoldChange)

ggplot(all_mix, aes(x = Mix_early_log2FoldChange, y = Mix_late_log2FoldChange, color = Phylum, label = Genus)) + geom_point(size = 1) +
  coord_flip() +
  geom_text(check_overlap = T, hjust = 0, nudge_x = 0.05) 

mix_OTUS <- Mix_early_data %>%
  select(OTU, Phylum, Genus, Mix_early_log2FoldChange = log2FoldChange) %>%
  full_join(Mix_late_data) %>%
  select(OTU, Phylum, Genus, Mix_early_log2FoldChange, Mix_late_log2FoldChange = log2FoldChange)

head(mix_OTUS)

saveRDS(mix_OTUS, file = "data/RDS/LFC_mix_OTUs_June19.RDS")
