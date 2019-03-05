# This object has day 0 and amendement samples removed and is rarefied to 6000
physeq <- readRDS("data/RDS/IncPhyseqRareClusteredTree")
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