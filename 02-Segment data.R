
#### Segment armadillo acceleration ####

library(bayesmove)
library(tidyverse)
library(lubridate)
library(future)

source('helper functions.R')


### Import Data ###
dat<- read.csv("Binned Armadillo Acceleration Data.csv", as.is = T)
dat$date<- as_datetime(dat$date)
# dat$day = yday(dat$date)

#Fix time1 var if not sequential after prior filtering steps
for (i in 1:dplyr::n_distinct(dat$id)) {
  ind<- which(dat$id == unique(dat$id)[i])
  
  if (dplyr::n_distinct(diff(dat[ind,]$time1)) > 1) 
    dat[ind,]$time1<- seq(1, length(ind))
}


#Create list
dat.list<- df_to_list(dat, "id")

#Determine active period
dat.list<- dat.list %>% 
  map(., nocturnal_period)


# #Remove days where SL & TA data is NA
# filter_by_day = function(data) {
#   
#   data.list<- df_to_list(data, "id")
#   dat.in<- list()  #data kept in for segmentation model
#   dat.out<- list()  #data held out for prediction
#   
#   for (i in 1:length(data.list)) {
#     cond<- unique(data.list[[i]][which(!is.na(data.list[[i]]$SL)), "day"])
#     
#     dat.in[[i]]<- data.list[[i]] %>%
#       filter(day %in% cond)
#     dat.out[[i]]<- data.list[[i]] %>%
#       filter(!day %in% cond)
#   }
# 
#   # dat.in<- bind_rows(dat.in)
#   # dat.out<- bind_rows(dat.out)
#   names(dat.in)<- names(dat.out)<- names(data.list)
#   
#   list(train = dat.in, test = dat.out)
# }
# 
# 
# dat.filt<- filter_by_day(data = dat)
# 
# #Compare AC distribs between train and test data
# tmp<- map(dat.filt, bind_rows) %>% 
#   map2(., c("train","test"), ~mutate(.x, type = .y)) %>% 
#   bind_rows()
# 
# ggplot(tmp, aes(Activity.Count)) +
#   geom_histogram() +
#   theme_bw() +
#   labs(x = "Activity Count", y = "Count") +
#   theme(axis.title = element_text(size = 16),
#         axis.text = element_text(size = 12),
#         strip.text = element_text(size = 14, face = "bold")) +
#   facet_wrap(~type, scales = "free_y")



# Only retain id and discretized data streams
# dat.list2<- map(dat.filt$train, subset, select = c(id, AC, SL, TA))
dat.list2<- map(dat.list, subset, select = c(id, AC, SL, TA))


### Run Segmentation Model ###

set.seed(1)

alpha<- 1
ngibbs<- 20000
nbins<- c(6,6,8)
breaks<- map(dat.list, ~find_breaks(., "period"))


plan(multisession)  #run all MCMC chains in parallel
dat.res<- segment_behavior(data = dat.list2, ngibbs = ngibbs, nbins = nbins, alpha = alpha,
                           breakpt = breaks)
future:::ClusterRegistry("stop")  #close all threads and memory used
# takes 55 min to run 20000 iterations
# takes 5 min for only SL and TA data
# takes 3.5 min on highly filtered data


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
