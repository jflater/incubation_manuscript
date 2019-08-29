library(phyloseq)
library(vegan)
library(tidyverse)
library(nlme)
library(emmeans)
library(ggpubr)
library(agricolae)
library(broom)
inc.raw.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

inc.physeq <- subset_samples(inc.raw.physeq, day %in% c("7",
                                                        "14",
                                                        "21",
                                                        "35",
                                                        "49",
                                                        "97"))

#Rename treatments to more informative titles
data <- data.frame(sample_data(inc.physeq)) %>%
  mutate(treatment = recode(treatment,
                            'Control' = 'Reference',
                            'CompAlfa' = 'Mix')) %>%
  mutate(C_N = C_flash / N_flash, Inorganic_N = NH3 + NO3) %>%
  mutate(TreatmentAndDay = paste(treatment, day))
data$treatment <- relevel(data$treatment, ref = "Reference") 
data$day <- as.factor(data$day)
rownames(data) <- data$i_id
sample_data(inc.physeq) <- data
sample_data(inc.physeq)$day <- as.factor(sample_data(inc.physeq)$day)
sample_data(inc.physeq)

inc.model.data <- lme(MBC_mg.kg_per_dry_wt_soil~treatment * day, random=~1|replication 
                      , data = data
                      , weights = varIdent(form= ~1|day*treatment)
                      , control = lmeControl(opt = "optim", msVerbose = TRUE))

em <- emmeans(inc.model.data, c("day", "treatment"), data = data)
plot(em)

#my.eff <- Effect(c("treatment", "day"), summary(inc.model.data))
#plot(my.eff)

sum_em <- summary(em)
plot(sum_em)

theme_set(theme_bw())
p <- ggplot(data = data, aes(x = day, y = MBC_mg.kg_per_dry_wt_soil     )) +
  geom_point(aes(colour = treatment), size = 4) +
  stat_summary(aes(group = treatment), fun.y = mean,  geom = "line", size = 2, colour = "steelblue") +
  geom_pointrange(size = 1, pch = 1, data = sum_em, aes(x = day, y = emmean, ymin = lower.CL, ymax = upper.CL, group = treatment)) +
  xlab("Day") +
  ylab("Microbial Biomass") +
  facet_wrap(~treatment) +
  theme(
    axis.title.y=element_text(colour = "black", size = 17, hjust = 0.5, margin=margin(0,12,0,0)),
    axis.title.x=element_text(colour = "black", size = 17),
    axis.text.x=element_text(colour = "black", size=15),
    axis.text.y=element_text(colour = "black", size=15),
    legend.position="none",
    legend.text=element_text(size=12.5),
    legend.key=element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text.x=element_text(size=15)
  )
png("Figures/MBC_mg.kg_per_dry_wt_soil.png",height=5,width=5,units='in',res=300)
p
dev.off()
p


em2 <- emmeans(inc.model.data, c("treatment", "day"), data = data)

lambdas <- list(
  "Alfalfa - Reference" = c(-1, 1, rep(0, 24 - 2))
  , "Mix - Reference" = c(-1, 0, 1, rep(0, 24 - 3))
  , "Compost - Reference" = c(-1, 0, 0, 1, rep(0, 24 - 4))
  , "Alfalfa - Reference" = c(rep(0, 4), -1, 1, rep(0, 24 - 6))
  , "Mix - Reference" = c(rep(0, 4), -1, 0, 1, rep(0, 24 - 7))
  , "Compost - Reference" = c(rep(0, 4), -1, 0, 0, 1, rep(0, 24 - 8))
  , "Alfalfa - Reference" = c(rep(0, 8), -1, 1, rep(0, 24 - 10))
  , "Mix - Reference" = c(rep(0, 8), -1, 0, 1, rep(0, 24 - 11))
  , "Compost - Reference" = c(rep(0, 8), -1, 0, 0, 1, rep(0, 24 - 12))
  , "Alfalfa - Reference" = c(rep(0, 12), -1, 1, rep(0, 24 - 14))
  , "Mix - Reference" = c(rep(0, 12), -1, 0, 1, rep(0, 24 - 15))
  , "Compost - Reference" = c(rep(0, 12), -1, 0, 0, 1, rep(0, 24 - 16))
  , "Alfalfa - Reference" = c(rep(0, 16), -1, 1, rep(0, 24 - 18))
  , "Mix - Reference" = c(rep(0, 16), -1, 0, 1, rep(0, 24 - 19))
  , "Compost - Reference" = c(rep(0, 16), -1, 0, 0, 1, rep(0, 24 - 20))
  , "Alfalfa - Reference" = c(rep(0, 20), -1, 1, rep(0, 24 - 22))
  , "Mix - Reference" = c(rep(0, 20), -1, 0, 1, rep(0, 24 - 23))
  , "Compost - Reference" = c(rep(0, 20), -1, 0, 0, 1)
)

sum_em2 <- summary(contrast(em2, lambdas), infer = c(TRUE, TRUE), adjust = "mvt")
sum_em2

sum_em2$day <- factor(rep(c(7, 14, 21, 35, 49, 97), each = 3))

theme_set(theme_bw())
p <- ggplot(data = sum_em2, aes(x = day, y = estimate)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL, group = contrast, colour = contrast), size = 0.7) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  xlab("Day") +
  ylab(paste('Microbial biomass \n difference from reference')) +
  scale_y_continuous(breaks = seq(-10, 800, 25)) +
  facet_wrap(~contrast) +
  theme(
    axis.title.y=element_text(colour = "black", size = 17, hjust = 0.5, margin=margin(0,12,0,0)),
    axis.title.x=element_text(colour = "black", size = 17),
    axis.text.x=element_text(colour = "black", size=15),
    axis.text.y=element_text(colour = "black", size=15),
    legend.position="none",
    legend.text=element_text(size=12.5),
    legend.key=element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text.x=element_text(size=15)
  )
png("Figures/MBC_mg.kg_per_dry_wt_soil_plot_diff.png",height=5,width=8,units='in',res=300)
p
dev.off()
p
