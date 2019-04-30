#####
inc.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

tree <- read.tree("data/tree.nwk")

inc.physeq <- merge_phyloseq(inc.physeq, tree)

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
#####
## Select early responding OTUs with greater than 2 lfc, use list generated in deseq_to_function
## two lists below are from deseq_to_function.R
early.alf.otus <- readRDS("data/early.alf.otus.rds")
late.alf.otus <- readRDS("data/late.alf.otus.rds")

library(gplots)

alf.venn <- venn(list(Early = early.alf.otus$OTU, Late = late.alf.otus$OTU))
attr(alf.venn, "intersections")$`Early:Late`
310 + 5 + 153
length(early.alf.otus$OTU)
length(late.alf.otus$OTU)
length(early.alf.otus$OTU) + length(late.alf.otus$OTU)

alf.list <- c(early.alf.otus$OTU, late.alf.otus$OTU)
alf.list

length(alf.list)
length(base::unique(alf.list))

alf.melt <- prune_taxa(alf.list, inc.physeq) %>%
  subset_samples(treatment %in% c("Alfalfa")) %>%
  filter_taxa(function(x) sum(x) >= 1, T) %>%
  psmelt()

both <- attr(alf.venn, "intersections")$`Early:Late`
early <- attr(alf.venn, "intersections")$Early
late <- attr(alf.venn, "intersections")$Late

alf.melt$group[alf.melt$OTU %in% both] <- "Early:Late"
alf.melt$group[alf.melt$OTU %in% early] <- "Early"
alf.melt$group[alf.melt$OTU %in% late] <- "Late"

p = ggplot(alf.melt, aes(x=day, y=Abundance, fill=group))
p = p + geom_bar(color="black", stat="identity", position="dodge") +
  facet_grid(. ~ Phylum)
print(p)

RelativeAbundanceDf <- function(physeq) {
  physeq %>% tax_glom(taxrank = "Phylum") %>% transform_sample_counts(function(x) {
    x/sum(x)
  }) %>% psmelt() %>% filter(Abundance > 0.02) %>% arrange(Phylum)
}

# Function to plot relative abundance
PlotRelativeAbundance <- function(df) {
  ggplot(df, aes(x = as.factor(treatment), y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity") +
    theme(axis.title.x = element_blank()) +
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
    ylab("Phylogenetic distribution of OTUs with abundance > 2%")
}

treatment_names <- c(
  `1` = "Alfalfa",
  `2` = "Mix",
  `3` = "Compost",
  `4` = "Reference"
)
day_names <- c(
  `1` = "7",
  `2` = "14",
  `3` = "21",
  `4` = "35",
  `5` = "49",
  `6` = "97"
)

rare.merged <- merge_samples(rare6k.physeq, "TreatmentAndDay")

sample_data(rare.merged)$TreatmentAndDay <- levels(sample_data(rare6k.physeq)$TreatmentAndDay)

relllll <- PlotRelativeAbundance(RelativeAbundanceDf(rare.merged)) +
  facet_grid(~ day, labeller = labeller(day = as_labeller(day_names))) +
  scale_x_discrete(labels = treatment_names) +
  rotate_x_text(angle = 45)
relllll



