set.seed(42)

calculate_rarefaction_curves <- function(psdata, measures, depths, parallel=FALSE) {
  require('plyr') # ldply
  require('reshape2') # melt
  require('doParallel')
  
  # set parallel options if required
  if (parallel) {
    paropts  <- list(.packages=c("phyloseq", "reshape2"))
  } else {
    paropts  <- NULL
  }
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive() && ! parallel, 'text', 'none'), .parallel=parallel, .paropts=paropts)
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

# Loading required library and displaying core configuration
pkgTest <- function(x) {
    if (!require(x,character.only = TRUE))
    {
     install.packages(x,dep=TRUE)
       if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }

pkgTest('doParallel')
library(doParallel)


#####
#Read in your physeq object, mine is inc.physeq
library(phyloseq)
library(tidyverse)

inc.physeq <- readRDS("data/RDS/incubation_physeq_Aug18.RDS")

inc.physeq <- subset_samples(inc.physeq, day %in% c("7",
                                                 "14",
                                                 "21",
                                                 "35",
                                                 "49",
                                                 "97")) %>%
  filter_taxa(function(x) sum(x) >= 3, T)

inc.physeq
min(taxa_sums(inc.physeq))

# Use multiple cores
detectCores(all.tests=TRUE)
# 48

# Setting up and registering the cluster
cl <- makeCluster(detectCores(all.tests=TRUE))
registerDoParallel(cl)

rarefaction_curve_data <- calculate_rarefaction_curves(inc.physeq, c('Observed', 'Shannon'), rep(c(1, 10, 100, 1000, 1:100 * 1000), each = 10))
#####

# Shut down cluster
stopCluster(cl)
rm(cl)
registerDoSEQ()

summary(rarefaction_curve_data)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(inc.physeq)), by.x = 'Sample', by.y = 'row.names')
library('ggplot2')

plot <- ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = treatment,
    group = Sample
  )
) + geom_line(
) + geom_pointrange(
) + facet_wrap(
  facets = ~ Measure,
  scales = 'free_y'
) + geom_vline(xintercept = 5000) +
  geom_point(size = .1)
plot$layers[[2]] <- NULL
plot 

png("Figures/rare_curve_days7_97_mintax_3.png",height=5,width=8,units='in',res=600)
plot(plot)
dev.off()
