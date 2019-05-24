inc.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

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

no.unclass <- subset_taxa(inc.physeq, !Phylum=="Bacteria_unclassified")
no.unclass <- subset_taxa(no.unclass, !Genus=="Gp6_unclassified")
inc.physeq <- no.unclass

alf.physeq <- subset_samples(inc.physeq, Treatment_Response %in% c("Alfalfa_early", "Reference_early")) %>%
  filter_taxa(function(x) sum(x) >= 3, T)

diagdds = phyloseq_to_deseq2(alf.physeq, ~ treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="local")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(alf.physeq)[rownames(sigtab), ], "matrix"))
head(sigtab)

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
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
