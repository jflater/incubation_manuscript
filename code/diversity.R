library(phyloseq)
library(ggplot2)
# This object has day 0 and amendement samples removed and is rarefied to 6000
# I think Chao and observed need singletons included
physeq <- readRDS("data/")
plot_richness(physeq)
shannon <- plot_richness(physeq, measures = "Shannon")
shannon.df <- shannon$data
# load summary function
source("Functions/summarySE.R")
summary.shannon.df <- summarySE(shannon.df, measurevar = "value", groupvars = c("treatment", "day"))
pd <- position_dodge(0.2) # move them .05 to the left and right
shannon.diversity <- ggplot(summary.shannon.df, aes(x=day
                               , y=value
                               , colour=treatment
                               , group = treatment)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se)
                , width=.1, position=pd) +
  geom_point() +
  geom_line() + 
  ggtitle("Shannon diversity")

png("Figures/shannon.diversity.png",height=4,width=9,units='in',res=300)
shannon.diversity + 
  ggtitle("Shannon diversity for incubated microcosm, rarefied to 6k")
dev.off()
#####
Chao1 <- plot_richness(physeq, measures = "Chao1")
Chao1.df <- Chao1$data
# load summary function
source("Functions/summarySE.R")
summary.Chao1.df <- summarySE(Chao1.df, measurevar = "value", groupvars = c("treatment", "day"))
pd <- position_dodge(0.2) # move them .05 to the left and right
Chao1.diversity <- ggplot(summary.Chao1.df, aes(x=day
                                                    , y=value
                                                    , colour=treatment
                                                    , group = treatment)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se)
                , width=.1, position=pd) +
  geom_point() +
  geom_line() + 
  ggtitle("Chao1 diversity")

png("Figures/Chao1.diversity.png",height=4,width=9,units='in',res=300)
Chao1.diversity + 
  ggtitle("Chao1 diversity for incubated microcosm, rarefied to 6k")
dev.off()
#####
Observed <- plot_richness(physeq, measures = "Observed")
Observed.df <- Observed$data
# load summary function
source("Functions/summarySE.R")
summary.Observed.df <- summarySE(Observed.df, measurevar = "value", groupvars = c("treatment", "day"))
pd <- position_dodge(0.2) # move them .05 to the left and right
Observed.diversity <- ggplot(summary.Observed.df, aes(x=day
                                                , y=value
                                                , colour=treatment
                                                , group = treatment)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se)
                , width=.1, position=pd) +
  geom_point() +
  geom_line() + 
  ggtitle("Observed diversity")

png("Figures/Observed.diversity.png",height=4,width=9,units='in',res=300)
Observed.diversity + 
  ggtitle("Observed diversity for incubated microcosm, rarefied to 6k")
dev.off()