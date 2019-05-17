library(ggtree)

# Example from ggtree data
beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)

genotype_file <- system.file("examples/Genotype.txt", package="ggtree")
genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)

p <- ggtree(beast_tree, layout='circular', color="#4DAF4A", size=2, branch.length='none', right=T) +
  annotate('text', x=0, y=40, label='ggtree', family='mono', size=16)

p2 <- gheatmap(p, genotype, width=0.2, hjust='left', colnames_angle=-10, font.size=1.5)  +
  scale_fill_manual(values=c("#E41A1C","#377EB8","#FC8D59")) + theme_tree()
 
open_tree(p2, 80) %>% rotate_tree(80)

# Incubation data
tree <- read.tree("data/tree.nwk")
inc.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")
inc <- merge_phyloseq(inc.physeq, tree)
inc <- subset_samples(inc, day %in% c("7", "14", "21", "35", "49", "97"))

# too big, need smaller tree, use a phyloseq object with tree added and then only use one phyla
TopNOTUs = names(sort(taxa_sums(inc), TRUE)[1:100])
inc10 = prune_taxa(TopNOTUs, inc) %>%
  filter_taxa(function(x) sum(x) >= 3, T)

library(tidyverse)
t1 <- phy_tree(inc10)
data <- psmelt(inc10) 
#### working here, need a data.frame with 10 rows for each OTU and then columns with more info

pp <- ggtree(t1, layout='circular', color="#4DAF4A", size=2, branch.length='none', right=T) 
pp

ph <- gheatmap(pp, cast.data, colnames = T, width=0.2, hjust='left', colnames_angle=-10, font.size=1.5)  +
   theme_tree() 
ph

