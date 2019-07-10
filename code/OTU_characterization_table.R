# Generate a table of OTU numbers for a set of phyloseq samples. In this case there will be: Total OTUs observed (totus), soil-enriched OTUs (sotus), amendment enriched (aotus), and responding OTUs (rotus)
# source the file containing functions to use
source("code/functions.R")

# Load libraries
if (!require(phyloseq)) install.packages('phyloseq')
library(phyloseq)
if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)
if (!require(phylosmith)) install.packages('phylosmith')
library(phylosmith)

# Read in phyloseq object for the incubation microcosm study
inc.physeq <- readRDS("../incubation_manuscript/data/RDS/incubation_physeq_Aug18.RDS")

# Suggest looking at sample_data(phy)
# Rename treatments to more informative titles, response.group determined with hclustering.R
data <- data.frame(sample_data(inc.physeq)) %>% 
  mutate(treatment = recode(treatment,
                            'Control' = 'Reference',
                            'CompAlfa' = 'Mix')) %>% 
  mutate(C_N = C_flash / N_flash, Inorganic_N = NH3 + NO3) %>%
  mutate(TreatmentAndDay = paste(treatment, day))

rownames(data) <- data$i_id
sample_data(inc.physeq) <- data
sample_data(inc.physeq)$day <- as.factor(sample_data(inc.physeq)$day)
sample_data(inc.physeq)$treatment <- as.character(sample_data(inc.physeq)$treatment)

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

# If you want to remove unclassified bacteria from your data set
no.unclass <- subset_taxa(inc.physeq, !Phylum=="Bacteria_unclassified")
no.unclass <- subset_taxa(no.unclass, !Genus=="Gp6_unclassified")

# Remove some stuff for now
rm(data)
rm(inc.physeq.data)

# Reference and amendment OTUs
refotus <- prune_samples(sample_data(inc.physeq)$treatment %in% c("Reference"), inc.physeq) %>%
  filter_taxa(function(x) sum(x) >= 5, T) %>%
  tax_table() %>%
  row.names()
# Alfalfa amendment OTUs
amendalfotus <- prune_samples(sample_data(inc.physeq)$treatment %in% c("AlfalfaAmend"), inc.physeq) %>%
  filter_taxa(function(x) sum(x) >= 5, T) %>%
  tax_table() %>%
  row.names()
amendcompotus <- prune_samples(sample_data(inc.physeq)$treatment %in% c("CompostAmend"), inc.physeq) %>%
  filter_taxa(function(x) sum(x) >= 5, T) %>%
  tax_table() %>%
  row.names()

# Now lets determine the number of OTUs in the early samples
early_alf_totus <- GetOTUs(inc.physeq, c("Alfalfa_early")) 
early_comp_totus <- GetOTUs(inc.physeq, c("Compost_early")) 
early_mix_totus <- GetOTUs(inc.physeq, c("Mix_early")) 
early_ref_totus <- GetOTUs(inc.physeq, c("Reference_early"))

# OTUs in early group orginating from soil
early_alf_sotus <- Reduce(intersect, list(early_alf_totus, refotus))
early_comp_sotus <- Reduce(intersect, list(early_comp_totus, refotus))
early_mix_sotus <- Reduce(intersect, list(early_mix_totus, refotus))

# OTUs in early group orginating from amendment, remember to change amendments as we go and we first need to combine alfala and compost because we don't have a sample of just those mixed then sequenced, only seperately  
amendmixotus <- unique(c(amendalfotus, amendcompotus))
early_alf_aotus <- Reduce(intersect, list(early_alf_totus, amendalfotus))
early_comp_aotus <- Reduce(intersect, list(early_comp_totus, amendcompotus))
early_mix_aotus <- Reduce(intersect, list(early_mix_totus, amendmixotus))

table <- data.frame("OTU.Group" = c("totus", "sotus", "aotus"),
                    "Early.Alfalfa" = c(length(early_alf_totus), length(early_alf_sotus), length(early_alf_aotus)),
                    "Early.Compost" = c(length(early_comp_totus), length(early_comp_sotus), length(early_comp_aotus)),
                    "Early.Mix" = c(length(early_mix_totus), length(early_mix_sotus), length(early_mix_aotus)))

view(table)
gridExtra::grid.table(table)
png("Figures/OTU_char_table_early.png", height = 1, width = 5.3, units = 'in', res = 600)
gridExtra::grid.table(table)
dev.off()

# Now for late
late_alf_totus <- GetOTUs(inc.physeq, c("Alfalfa_late")) 
late_comp_totus <- GetOTUs(inc.physeq, c("Compost_late")) 
late_mix_totus <- GetOTUs(inc.physeq, c("Mix_late")) 
late_ref_totus <- GetOTUs(inc.physeq, c("Reference_late"))

# OTUs in late group orginating from soil
late_alf_sotus <- Reduce(intersect, list(late_alf_totus, refotus))
late_comp_sotus <- Reduce(intersect, list(late_comp_totus, refotus))
late_mix_sotus <- Reduce(intersect, list(late_mix_totus, refotus))

# OTUs in late group orginating from amendment, remember to change amendments as we go and we first need to combine alfala and compost because we don't have a sample of just those mixed then sequenced, only seperately  
amendmixotus <- unique(c(amendalfotus, amendcompotus))
late_alf_aotus <- Reduce(intersect, list(late_alf_totus, amendalfotus))
late_comp_aotus <- Reduce(intersect, list(late_comp_totus, amendcompotus))
late_mix_aotus <- Reduce(intersect, list(late_mix_totus, amendmixotus))

table.2 <- data.frame("OTU.Group" = c("totus", "sotus", "aotus"),
                    "Late.Alfalfa" = c(length(late_alf_totus), length(late_alf_sotus), length(late_alf_aotus)),
                    "Late.Compost" = c(length(late_comp_totus), length(late_comp_sotus), length(late_comp_aotus)),
                    "Late.Mix" = c(length(late_mix_totus), length(late_mix_sotus), length(late_mix_aotus)))

view(table.2)
gridExtra::grid.table(table.2)
png("Figures/OTU_char_table_late.png", height = 1, width = 5.3, units = 'in', res = 600)
gridExtra::grid.table(table.2)
dev.off()
