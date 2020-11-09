
#### Segment armadillo acceleration ####

library(bayesmove)
library(tidyverse)
library(lubridate)
library(future)

### Import Data ###
dat<- read.csv("Binned Armadillo Acceleration Data.csv", as.is = T)
dat$date<- as_datetime(dat$date)
dat.list<- df_to_list(dat, "id")

# #Remove days where SL & TA data is NA
# filter_by_day = function(data) {
#   cond<- unique(data[which(!is.na(data$SL)), "day"])
#   data<- data %>% 
#     filter(day %in% cond)
#   
#   data
# }
# 
# dat.list<- df_to_list(dat, "id") %>% 
#   map(., filter_by_day)

# Only retain id and discretized data streams
dat.list2<- map(dat.list, subset, select = c(id, AC, SL, TA))



### Run Segmentation Model ###

set.seed(1)

alpha<- 1
ngibbs<- 20000
nbins<- c(5,5,8)
breaks<- map(dat.list, ~find_breaks(., "day"))


plan(multisession)  #run all MCMC chains in parallel
dat.res<- segment_behavior(data = dat.list2, ngibbs = ngibbs, nbins = nbins, alpha = alpha,
                           breakpt = breaks)
future:::ClusterRegistry("stop")  #close all threads and memory used
# takes 37 min to run 20000 iterations
# takes 5 min for only SL and TA data
# takes 56 min for filtered AC, SL, TA data to run 20000 iterations


# Trace-plots for the number of breakpoints and LML
traceplot(data = dat.res$nbrks, ngibbs = ngibbs, type = "nbrks")
traceplot(data = dat.res$LML, ngibbs = ngibbs, type = "LML")

## Based on traceplots of LML, all reached convergence



# Determine MAP for selecting breakpoints
MAP.est<- get_MAP(dat = dat.res$LML, nburn = ngibbs/2)
brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)

# How many breakpoints estimated per ID?
apply(brkpts[,-1], 1, function(x) length(purrr::discard(x, is.na)))


# Plot breakpoints over the data
plot_breakpoints(data = dat.list, as_date = FALSE, var_names = c("AC","SL","TA"),
                 var_labels = c('AC Bins','SL Bins', 'TA Bins'), brkpts = brkpts)

plot_breakpoints(data = dat.list, as_date = TRUE, var_names = c("Activity.Count","step","angle"),
                 var_labels = c('Activity Counts','Step Length (m)', 'Turning Angle (rad)'),
                 brkpts = brkpts)



### Export Results ###

# Assign track segments to all observations by ID
dat.seg<- assign_tseg(dat = dat.list, brkpts = brkpts)

# write.csv(dat.seg, "Segmented Armadillo Acceleration.csv", row.names = F)
