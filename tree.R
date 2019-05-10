library(ggtree)

beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)

genotype_file <- system.file("examples/Genotype.txt", package="ggtree")
genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)

p <- ggtree(beast_tree, layout='circular', color="#4DAF4A", size=2, branch.length='none', right=T) +
  annotate('text', x=0, y=40, label='ggtree', family='mono', size=16)

tree <- read.tree("data/tree.nwk")

# too big, need smaller tree, use a phyloseq object with tree added and then only use one phyla

pp <- ggtree(tree, layout='circular', color="#4DAF4A", size=2, branch.length='none', right=T) +
  annotate('text', x=0, y=40, label='ggtree', family='mono', size=16)

p2 <- gheatmap(p, genotype, width=0.2, hjust='left', colnames_angle=-10, font.size=1.5)  +
  scale_fill_manual(values=c("#E41A1C","#377EB8","#FC8D59")) + theme_tree()

open_tree(p2, 80) %>% rotate_tree(80)