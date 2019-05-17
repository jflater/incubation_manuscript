library(phyloseq)
library(biomformat)
# download source code for microfiltR, or use copied script below
source("code/microfiltR_source_code.R")
inc_phy <- readRDS("data/RDS/IncPhyseqRareClusteredTree")
inc_phy

write.dataset.biom <- function(ps, filePATH, filePREFIX, writeFASTA=TRUE, rename=FALSE, useREFSEQ=FALSE){
  
  #pull seqs from refseq slot or extract from ASV ID for fasta format
  if (isTRUE(useREFSEQ)){
    #from phyloseq refseq slot
    f.onames <- phyloseq::refseq(ps)
  } else {
    f.onames <- phyloseq::taxa_names(ps)
  }
  
  if (isTRUE(rename)){
    phyloseq::taxa_names(ps) <- paste("ASV", 1:length(phyloseq::taxa_names(ps)), sep = "")
    names(f.onames) <- paste0(">", phyloseq::taxa_names(ps))
  } else {
    names(f.onames) <- paste0(">", phyloseq::taxa_names(ps))
  }
  
  #generate biom file
  suppressWarnings(ps.b <- biomformat::make_biom(
    data = format.ASV.tab(ps),
    sample_metadata = as.data.frame(phyloseq::sample_data(ps)),
    observation_metadata = as.data.frame(phyloseq::tax_table(ps)), 
    matrix_element_type = "int"
  )
  )
  
  #create output string
  if (isTRUE(writeFASTA)){
    fa <- print(paste0(filePATH, filePREFIX, "_ASVs.fasta"))
  }
  bo <- print(paste0(filePATH, filePREFIX, "_ASV_table.biom"))
  
  #write output
  if (isTRUE(writeFASTA)){
    write.table(x = f.onames, file = fa, quote = FALSE, sep = "\n", col.names = FALSE)
  }
  #biom export
  biomformat::write_biom(x = ps.b, biom_file = bo)
  
  #return phyloseq object with taxa renamed to ASV1, etc.
  return(ps)
}
biom <- write.dataset.biom(inc_phy, filePATH = "data/", writeFASTA = F, filePREFIX = F, rename = F, useREFSEQ = F)
test <- make_biom(data = inc_phy@otu_table, sample_metadata = inc_phy@sam_data, observation_metadata = inc_phy@tax_table)
biom_inc <- write_biom(x = test, biom_file = "data/inc_biom.biom")
## Now use this file and the export to graphlan python script to get into graphlan format

OTU10 = names(sort(taxa_sums(), TRUE)[1:10])
b10=prune_taxa(OTU10, b1)
