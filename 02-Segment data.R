
#### Segment armadillo acceleration ####

library(bayesmove)
library(tidyverse)
library(lubridate)
library(future)

### Import Data ###
dat<- read.csv("Binned Armadillo Acceleration Data.csv", as.is = T)
dat.list<- df_to_list(dat, "id")

# Only retain id and discretized step length (SL) and turning angle (TA) columns
dat.list2<- map(dat.list, subset, select = c(id, AC))



### Run Segmentation Model ###

set.seed(1)

alpha<- 1
ngibbs<- 15000
nbins<- 6


plan(multisession)  #run all MCMC chains in parallel
dat.res<- segment_behavior(data = dat.list2, ngibbs = ngibbs, nbins = nbins, alpha = alpha)
future:::ClusterRegistry("stop")  #close all threads and memory used
# takes 17.5 min to run 15000 iterations


# Trace-plots for the number of breakpoints and LML
traceplot(data = dat.res$nbrks, ngibbs = ngibbs, type = "nbrks")
traceplot(data = dat.res$LML, ngibbs = ngibbs, type = "LML")

## Based on traceplots of LML, all reached convergence



# Determine MAP for selecting breakpoints
MAP.est<- get_MAP(dat = dat.res$LML, nburn = 7500)
brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)

# How many breakpoints estimated per ID?
apply(brkpts[,-1], 1, function(x) length(purrr::discard(x, is.na)))


# Plot breakpoints over the data
plot_heatmap(data = dat.list2, nbins = nbins, brkpts = brkpts, title = TRUE, legend = TRUE)




### Export Results ###

# Assign track segments to all observations by ID
dat.seg<- assign_tseg(dat = dat.list, brkpts = brkpts)

# write.csv(dat.seg, "Segmented Armadillo Acceleration.csv", row.names = F)
