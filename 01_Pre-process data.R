
#### Pre-process armadillo data for bayesmove behavior model ####

library(bayesmove)
library(tidyverse)
library(lubridate)
library(furrr)
library(future)
library(tictoc)

source('helper functions.R')


### Import Armadillo Acceleration Data ###
dat<- list.files(pattern = "*.csv") %>% 
  map_df(~read.csv(.))


#inspect data
str(dat)
summary(dat)

#inspect what's happening with NA Activity Counts
ind.na<- which(is.na(dat$Activity.Count))
head(dat[ind.na,])





### Prep Data for Model ###

dat2<- dat %>% 
  dplyr::select(-c(X,X.1)) %>%  #remove rowname cols
  mutate_at("Acquisition.Time",
            ~as.POSIXct(strptime(., format = "%Y.%m.%d %H:%M:%S"))) %>%  #make date POSIXct format
  mutate_at("Acquisition.Start.Time",
            ~as.POSIXct(strptime(., format = "%Y.%m.%d %H:%M:%S"))) %>%  #make date POSIXct format
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

table(dat.AC$dt)  #all but 5 of 140166 (0.004%) are at 5 min interval

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
table(dat.GPS$dt)  #all but 554 of 99290 (0.56%) are at 7 min interval; most others are at 12 min interval


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



## Originally started w/ 140166 activity count obs and 18433 GPS obs, which was filtered to retain a total of 67872 obs




### Define Bin Limits ###

# Activity Counts (start w/ 5)
act.bin.lims=quantile(dat.merged3$Activity.Count,
                       c(0,0.75,0.80,0.85,0.95,1), na.rm=T)  #5 bins
#values from 0-18 assumed to represent little to no actual movement; reflective of remaining within a burrow

# Turning Angles (start w/ 8)
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

# Step Lengths (start w/ 5)
dist.bin.lims=quantile(dat.merged3$step,
                       c(0,0.50,0.75,0.95,0.99,1), na.rm=T) #5 bins



# Viz density distribs w/ bin limits
ggplot(dat.merged3, aes(Activity.Count)) +
  geom_density(fill = "cadetblue") +
  geom_vline(data = data.frame(lims = act.bin.lims), aes(xintercept = lims), 
             linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


ggplot(dat.merged3, aes(step)) +
  geom_density(fill = "lightblue") +
  geom_vline(data = data.frame(lims = dist.bin.lims), aes(xintercept = lims), 
             linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


ggplot(dat.merged3, aes(angle)) +
  geom_density(fill = "indianred") +
  geom_vline(data = data.frame(lims = angle.bin.lims), aes(xintercept = lims), 
             linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))



# Assign bins to observations
dat.list<- df_to_list(dat.merged3, "id")
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
  drop_na()

param.prop[1:5, "value"]<- ((diff(act.bin.lims)/2) + act.bin.lims[1:5])
param.prop[6:10, "value"]<- ((diff(dist.bin.lims)/2) + dist.bin.lims[1:5])
param.prop[11:18, "value"]<- (diff(angle.bin.lims)/2) + angle.bin.lims[1:8]

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
