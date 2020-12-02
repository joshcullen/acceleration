
#### Pre-process armadillo data for bayesmove behavior model ####

library(bayesmove)
library(tidyverse)
library(lubridate)
library(furrr)
library(future)
library(tictoc)

source('helper functions.R')


### Import Armadillo Acceleration Data ###
files<- list.files(pattern = "*.csv")
files<- files[!grepl(pattern = "Armadillo", files)]
dat<- map_df(files, ~read.csv(.))


#inspect data
str(dat)
summary(dat)






### Prep Data for Model ###

dat2<- dat %>% 
  dplyr::select(-c(X,X.1)) %>%  #remove rowname cols
  mutate_at("Acquisition.Time",
            ~as.POSIXct(strptime(., format = "%Y.%m.%d %H:%M:%S", tz = "UTC"))) %>%  #make date POSIXct format
  mutate_at("Acquisition.Start.Time",
            ~as.POSIXct(strptime(., format = "%Y.%m.%d %H:%M:%S", tz = "UTC"))) %>%  #make date POSIXct format
  mutate_at("GPS.Fix.Time",
            ~as.POSIXct(strptime(., format = "%Y.%m.%d %H:%M:%S", tz = "UTC")))  #make date POSIXct format


#create data frame for activity counts and add dt column (units as min)
dat.AC<- dat2 %>% 
  drop_na(Activity.Count) %>%  #remove NA rows for AC
  group_split(id) %>% 
  map(., arrange, Acquisition.Time) %>%  #reorder rows by Acquisition.Time
  map(., distinct, Acquisition.Time, .keep_all = T) %>% #remove duplicate rows
  map(., ~mutate(., dt = difftime(Acquisition.Time, Acquisition.Start.Time, units = "min"),
                 .before = Activity.Count)) %>% 
  bind_rows()

table(dat.AC$dt)  #all but 4 of 140166 (0.004%) are at 5 min interval

dat.AC.list<- df_to_list(dat.AC, "id")
dat.AC.filt<- filter_time(dat.AC.list, int = 5)


#create data frame for activity counts and add dt column (units as min)
dat.GPS<- dat2 %>% 
  filter(is.na(Activity.Count)) %>%  #remove rows where AC is NA
  group_split(id) %>% 
  map(., arrange, Acquisition.Start.Time) %>%  #reorder rows by Acquisition.Time
  map(., distinct, Acquisition.Start.Time, .keep_all = T) %>% #remove duplicate rows
  map(., ~mutate(., dt = as.numeric(difftime(Acquisition.Start.Time, lag(Acquisition.Start.Time),
                                             units = "min")),
                 .before = Activity.Count)) %>%
  bind_rows() %>% 
  dplyr::select(-date) %>% 
  rename(date = Acquisition.Start.Time)


# Fill in large time gaps with NA and round times
dat.GPS<- dat.GPS %>% 
  fill_NA(data = ., int = 7) %>% 
  round_track_time(dat = ., id = "id", int = 7, tol = 1)  #round dt
table(dat.GPS$dt)  #all but 554 of 99282 (0.56%) are at 7 min interval; most others are at 12 min interval


# Filter times to time interval
dat.GPS.list<- df_to_list(dat.GPS, "id")
dat.GPS.filt<- filter_time(dat.GPS.list, int = 7)



### Merge AC and GPS data ###

plan(multisession)
tic()
dat.merged<- future_map2(dat.AC.filt, dat.GPS.filt,
                  ~merge_ac_gps(data.AC = .x, data.GPS = .y, tol = 1),
                  .progress = TRUE, .options = future_options(seed = TRUE)) %>% 
  bind_rows()
toc()
future:::ClusterRegistry("stop")  #close all threads and memory used
# takes 14 s to run for 99290 GPS observations



# Remove all rows where both AC and GPS coords are NA
dat.merged<- dat.merged %>% 
  filter(!is.na(Activity.Count) | !is.na(GPS.UTM.Easting))




### Inspect Data to Check for Aberrant Observations ###

# Viz distribution and time series of Activity.Count
ggplot(dat.merged, aes(Activity.Count)) +
  geom_density(aes(fill = id)) +
  theme_bw() +
  facet_wrap(~ id)

ggplot() +
  geom_path(data = dat.merged, aes(time1, Activity.Count, color = id)) +
  theme_bw() +
  facet_wrap(~id, scales = "free")

# All IDs show highly right-skewed distributions; use quantiles to define bin limits




# Calculate step lengths and turning angles
dat.merged2<- prep_data(dat = dat.merged, coord.names = c('GPS.UTM.Easting','GPS.UTM.Northing'),
                        id = "id")



# Viz distribution and time series of step lengths
ggplot(dat.merged2, aes(step)) +
  geom_density(aes(fill = id)) +
  theme_bw() +
  facet_wrap(~ id)

ggplot() +
  geom_path(data = dat.merged2, aes(time1, step, color = id)) +
  theme_bw() +
  facet_wrap(~id, scales = "free")

# Viz map of tracks to pickup on unusual movement patterns
ggplot() +
  geom_path(data = dat.merged2, aes(x, y, color = id)) + 
  theme_bw() +
  facet_wrap(~id, scales = "free")
# appears to be issues w/ some excessively long step lengths and outlier observations; these need to be removed



# Filter out steps > 2000 m
dat.merged3<- dat.merged2 %>% 
  filter(step <= 2000 | is.na(step))


# Viz map of tracks to pickup on unusual movement patterns
ggplot() +
  geom_path(data = dat.merged3, aes(x, y, color = id)) + 
  theme_bw() +
  facet_wrap(~id, scales = "free")
# still some outlier observations; these need to be indexed directly




## Index outliers for "tex" and "mazeboti"
ind<- c(which(dat.merged3$id == "tex" & dat.merged3$x < 600000),
        which(dat.merged3$id == "mazeboti" & dat.merged3$x < 600000 | dat.merged3$x > 650000))
dat.merged3<- dat.merged3[-ind,]



# Viz map of tracks to pickup on unusual movement patterns
ggplot(data = dat.merged3, aes(x, y)) +
  geom_point(aes(color = id), alpha = 0.6) +
  geom_path() + 
  theme_bw() +
  facet_wrap(~id, scales = "free")

ggplot(data = dat.merged3, aes(x, y)) +
  geom_path(aes(group = id, color = id)) + 
  geom_point(aes(color = id), alpha = 0.4) +
  theme_bw()
# looks much better, but may still want to filter out more observations later; stick w/ this for now



# Viz distribution and time series of step lengths
ggplot(dat.merged3, aes(step)) +
  geom_density(aes(fill = id)) +
  theme_bw() +
  facet_wrap(~ id)

ggplot() +
  geom_path(data = dat.merged3, aes(time1, step, color = id)) +
  theme_bw() +
  facet_wrap(~id, scales = "free_x")



#Viz distribution of turning angles
ggplot(dat.merged3, aes(angle)) +
  geom_density(aes(fill = id)) +
  theme_bw() +
  facet_wrap(~ id)

# Viz time series of turning angles
ggplot() +
  geom_path(data = dat.merged3, aes(time1, angle, color = id)) +
  theme_bw() +
  facet_wrap(~id, scales = "free")





## Need to remove observations where large consecutive number of observations have 0 for AC or that SL and TA are exactly equal to 0

## Index time periods to remove for "blanca", "gala", and "mazeboti"
cond1<- c(which(dat.merged3$id == "blanca" & 
                 dat.merged3$date >= as.POSIXct("2019-08-15 0:00:00", tz="UTC") & 
                 dat.merged3$date <= as.POSIXct("2019-11-09 0:00:00", tz="UTC")),
        which(dat.merged3$id == "gala" & 
                dat.merged3$date >= as.POSIXct("2019-08-24 0:00:00", tz="UTC") & 
                dat.merged3$date <= as.POSIXct("2019-10-01 0:00:00", tz="UTC")),
        which(dat.merged3$id == "mazeboti" & 
                dat.merged3$date >= as.POSIXct("2019-07-25 0:00:00", tz="UTC") & 
                dat.merged3$date <= as.POSIXct("2019-12-15 0:00:00", tz="UTC")))

cond2<- which(dat.merged3$Activity.Count == 0 & dat.merged3$step > 0 | dat.merged3$angle == 0)
cond3<- which(dat.merged3$step == 0 & dat.merged3$angle == 0)
# cond4<- which(!is.na(dat.merged3$Activity.Count) & is.na(dat.merged3$step) & 
#                 is.na(dat.merged3$angle))
# cond5<- which(is.na(dat.merged3$Activity.Count) & is.na(dat.merged3$step) & 
#                 is.na(dat.merged3$angle))

cond<- sort(unique(c(cond1, cond2, cond3)))
dat.merged4<- dat.merged3[-cond,]


# Viz map of tracks to pickup on unusual movement patterns
ggplot(data = dat.merged4, aes(x, y)) +
  geom_point(aes(color = id), alpha = 0.6) +
  geom_path() + 
  theme_bw() +
  facet_wrap(~id, scales = "free")

ggplot(data = dat.merged4, aes(x, y)) +
  geom_path(aes(group = id, color = id)) + 
  geom_point(aes(color = id), alpha = 0.4) +
  theme_bw()
# looks generally the same as before, but w/o spurious relocations


# Viz distribution and time series of step lengths
ggplot(dat.merged4, aes(Activity.Count)) +
  geom_density(aes(fill = id)) +
  theme_bw() +
  facet_wrap(~ id)

ggplot() +
  geom_path(data = dat.merged4, aes(time1, Activity.Count, color = id)) +
  theme_bw() +
  facet_wrap(~id, scales = "free_x")


# Viz distribution and time series of step lengths
ggplot(dat.merged4, aes(step)) +
  geom_density(aes(fill = id)) +
  theme_bw() +
  facet_wrap(~ id)

ggplot() +
  geom_path(data = dat.merged4, aes(time1, step, color = id)) +
  theme_bw() +
  facet_wrap(~id, scales = "free_x")


#Viz distribution of turning angles
ggplot(dat.merged4, aes(angle)) +
  geom_density(aes(fill = id)) +
  theme_bw() +
  facet_wrap(~ id)
#looks much better

# Viz time series of turning angles
ggplot() +
  geom_path(data = dat.merged4, aes(time1, angle, color = id)) +
  theme_bw() +
  facet_wrap(~id, scales = "free")



## Originally started w/ 140166 activity count obs and 18433 GPS obs, which was filtered to retain a total of 57643 obs




## Filter out days w/ no GPS locs
dat.merged4$day = yday(dat.merged4$date)

#Remove days where `step` and `angle` are NA
filter_by_day = function(data) {
  
  data.list<- df_to_list(data, "id")
  dat.in<- list()  #data kept in for segmentation model
  dat.out<- list()  #data held out for prediction
  
  for (i in 1:length(data.list)) {
    cond<- unique(data.list[[i]][which(!is.na(data.list[[i]]$step)), "day"])
    
    dat.in[[i]]<- data.list[[i]] %>%
      filter(day %in% cond)
    dat.out[[i]]<- data.list[[i]] %>%
      filter(!day %in% cond)
  }
  
  # dat.in<- bind_rows(dat.in)
  # dat.out<- bind_rows(dat.out)
  names(dat.in)<- names(dat.out)<- names(data.list)
  
  list(train = dat.in, test = dat.out)
}


dat.train_test<- filter_by_day(data = dat.merged4)  #filter for days w/ GPS data
dat.merged5<- dat.train_test$train %>%   #keep training data for analysis
  bind_rows()



## Find the distance between last and first recorded locations -- does the armadillo stay in same spot for ~ 20 hrs

#Find which locations to compare
ind.row<- which(!is.na(diff(dat.merged5$x)))
ind.diff<- which(diff(ind.row) > 10)

burrow.ind<- ind.row[sort(c(ind.diff, ind.diff+1))]
ind.ind<- seq(1, length(burrow.ind), by=2)  #need to index odd values of index
burrow.ind[ind.ind]<- burrow.ind[ind.ind] + 1  #need to increase end observation by 1 row 
burrow.locs<- dat.merged5[burrow.ind,] %>% 
  mutate(loc.type = rep(c("end","start"), nrow(.)/2))

# Viz map of locs
ggplot(data = burrow.locs, aes(x, y)) +
  geom_point(aes(color = loc.type), alpha = 0.6) +
  geom_path() + 
  theme_bw() +
  facet_wrap(~id, scales = "free")

#Calc distance between pairs of burrow locs
burrow.dist<- vector()
for (i in seq(2, nrow(burrow.locs), by = 2)) {
  tmp<- sqrt((burrow.locs$x[i] - burrow.locs$x[i-1])^2 + 
             (burrow.locs$y[i] - burrow.locs$y[i-1])^2)
  
  burrow.dist<- c(burrow.dist, tmp)
}

plot(density(burrow.dist))

## Apart from a few large outliers (likely due to the differences found between IDs when going through the DF), most differences in locs are relatively small. We could potentially assume that any observations on days where GPS measurements were recorded and Activity Counts are < 5 that armadillos are encamped in the burrow. We could then assign these observations with step length = 0 and TA ~ Unif(-pi, pi). This could potentially solve issue of Activity.Count data overwhelming these periods of time without tossing out all data w/o GPS measurements.



# ### Generate 0s for SL and random TA
# 
# #Find when there are >=10 consecutive NAs for step length
# is.na.rle <- rle(is.na(dat.merged5[,"step"]))
# is.na.rle$values <- is.na.rle$values & is.na.rle$lengths >= 10  #when >10 consecutive NAs
# ind<- which(inverse.rle(is.na.rle) == TRUE)
# 
# #Find which step length NAs are associated w/ activity count <= 5
# ind2<- which(dat.merged5[ind, "Activity.Count"] <= 5)
# 
# #Generate values for SL and TA
# dat.merged5[ind2, "step"]<- 0
# dat.merged5[ind2, "angle"]<- runif(length(ind2), min = -pi, max = pi)


### Define Bin Limits ###

# Activity Counts (start w/ 6)
act.bin.lims=quantile(dat.merged5[dat.merged5$Activity.Count != 0, "Activity.Count"],
                      c(0.25,0.50,0.65,0.75,0.90,1), na.rm=T)
act.bin.lims<- round(c(0, act.bin.lims))  #6 bins
# act.bin.lims=quantile(dat.merged4$Activity.Count, c(0,0.25,0.50,0.75,1), na.rm=T)  #4 bins

# Turning Angles (start w/ 8)
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

# Step Lengths (start w/ 6)
dist.bin.lims=quantile(dat.merged5[dat.merged5$step != 0, "step"],
                       c(0.25,0.50,0.75,0.95,0.99,1), na.rm=T) 
dist.bin.lims=c(0, dist.bin.lims) #6 bins



# Viz density distribs w/ bin limits
ggplot(dat.merged4, aes(Activity.Count)) +
  geom_density(fill = "cadetblue") +
  geom_vline(data = data.frame(lims = act.bin.lims), aes(xintercept = lims), 
             linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


ggplot(dat.merged4, aes(step)) +
  geom_density(fill = "lightblue") +
  geom_vline(data = data.frame(lims = dist.bin.lims), aes(xintercept = lims), 
             linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


ggplot(dat.merged4, aes(angle)) +
  geom_density(fill = "indianred") +
  geom_vline(data = data.frame(lims = angle.bin.lims), aes(xintercept = lims), 
             linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))






# Assign bins to observations
dat.list<- df_to_list(dat.merged5, "id")
dat_disc.list<- map(dat.list,
                       discrete_move_var,
                       lims = list(act.bin.lims, dist.bin.lims, angle.bin.lims),
                       varIn = c("Activity.Count", "step", "angle"),
                       varOut = c("AC","SL","TA"))


dat_disc.list2<- map(dat_disc.list, subset, select = c(id, AC, SL, TA))


#Viz discretized distrib
dat_disc.df<- dat_disc.list2 %>% 
  bind_rows() %>% 
  gather(key, value, -id)

param.prop<- dat_disc.df %>%
  group_by(key, value) %>%
  summarise(n=n()) %>%
  mutate(prop=n/nrow(dat_disc.df)) %>%
  ungroup() %>%   #if don't ungroup after grouping, ggforce won't work
  drop_na() %>% 
  data.frame()

param.prop[1:6, "value"]<- ((diff(act.bin.lims)/2) + act.bin.lims[1:6])
param.prop[7:12, "value"]<- ((diff(dist.bin.lims)/2) + dist.bin.lims[1:6])
param.prop[13:20, "value"]<- (diff(angle.bin.lims)/2) + angle.bin.lims[1:8]

ggplot(data = param.prop %>% filter(key == "AC"), aes(value, prop)) +
  geom_bar(stat = "identity", width = diff(act.bin.lims),
           fill = "cadetblue", color = "black") +
  labs(x = "Activity Count", y = "Proportion") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))

ggplot(data = param.prop %>% filter(key == "SL"), aes(value, prop)) +
  geom_bar(stat = "identity", width = diff(dist.bin.lims),
           fill = "lightblue", color = "black") +
  labs(x = "Step Length (m)", y = "Proportion") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))

ggplot(data = param.prop %>% filter(key == "TA"), aes(value, prop)) +
  geom_bar(stat = "identity", width = diff(angle.bin.lims),
           fill = "indianred", color = "black") +
  labs(x = "Turning Angle (rad)", y = "Proportion") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))





### Export Data ###
dat.out<- bind_rows(dat_disc.list)

# write.csv(dat.out, "Binned Armadillo Acceleration Data.csv", row.names = F)
