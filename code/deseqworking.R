physeq <- subset_samples(inc.physeq, Treatment_Response %in% c("Alfalfa_early", "Reference_early")) %>%
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
plot <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
plot
