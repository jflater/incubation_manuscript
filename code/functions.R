# A file containing commonly used functions created or modified by Jared Flater.
# July 10 2019

## Functions
#
#
# Get a list of OTUs from a phyloseq object, requires subseting to desired samples first, this will remove OTUs observed less than 5 times in the group of samples
# Note the column Treatment_Response is unique to incubation microcosm data
GetOTUs <- function(physeq, samples) {
    prune_samples(sample_data(physeq)$Treatment_Response %in% c(samples), physeq) %>%
    filter_taxa(function(x) sum(x) >= 5, T) %>%
    tax_table() %>%
    row.names()
}
#
#
# A function that takes a list of OTUs and returns the average relative abundance for those OTUs by treatment group
RelaOTUs <- function(physeq, samples, OTUs) {
  prune_samples(sample_data(physeq)$Treatment_Response %in% c(samples), physeq) %>% # run physeq
    subset_taxa(rownames(tax_table(physeq)) %in% c(OTUs)) %>%
  #  filter_taxa(function(x) sum(x) >= 1, T) %>% # Remove observed less than five times
    transform_sample_counts(function(x) x / sum(x)) %>% # Transform to relative abundance
    psmelt() %>%
    select(OTU, Abundance) %>%
    group_by(OTU) %>%
    summarise_at(.vars = names(.)[2],.funs = c(mean="mean"))
}
#
#
# LFC calculation function
who_diff_day <- function(DDS, choice1, choice2, phy.object){
  res = results(DDS, contrast = c("Treatment_Response", choice1, choice2), cooksCutoff = FALSE)
  #plotCounts(DDS, gene="OTU_311", intgroup="day")
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
#
#
# function plot log2FoldChange 
log_plot <- function(sigtab,t1){
  sigtab <- sigtab %>%
    rownames_to_column(var = "OTU") %>%
    filter(log2FoldChange >= 2) 
  
  ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + 
    coord_flip() +
    ggtitle(t1)
} 