
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





### Clean Data ###

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

dat.GPS<- round_track_time(dat = dat.GPS, id = "id", int = 7, tol = 1)  #round dt

table(dat.GPS$dt)  #all but 3078 of 1072 (0.02%) are at 7 min interval; most others are at 14, 21, or 28 min intervals due to failed signals or at 663 min during daytime period

dat.GPS.list<- df_to_list(dat.GPS, "id")
dat.GPS.filt<- filter_time(dat.GPS.list, int = 7)


plan(multisession)
tic()
dat.merged<- future_map2(dat.AC.filt, dat.GPS.filt,
                  ~merge_ac_gps(data.AC = .x, data.GPS = .y, tol = 1),
                  .progress = TRUE, .options = future_options(seed = TRUE)) %>% 
  bind_rows()
toc()
# takes 6 s to run for 49929 GPS observations



#### STEPS FOR MOVING FORWARD TO COMBINE ACTIVITY COUNTS AND GPS DATA ####
# 1) Try to include observations within long, daytime time gaps where AC is measured but GPS is not available; could provide further support for "Burrow" state






# Viz distribution of Activity.Count
ggplot(dat.merged, aes(Activity.Count)) +
  geom_density(aes(fill = id)) +
  theme_bw() +
  facet_wrap(~ id)

# All IDs show highly right-skewed distributions; use quantiles to define bin limits



#### NEED TO UPDATE FROM HERE; THIS INCLUDES CALCULATION OF STEP LENGTHS AND TURNING ANGLES TO BE BINNED ALONG WITH ACTIVITY COUNTS FOR ANALYSIS BY MODEL ####

### Define Bin Limits ###

# Define bin number and limits (start w/ 6)
act.bin.lims=quantile(filtered.df$Activity.Count,
                       c(0,0.50,0.75,0.80,0.85,0.95,1), na.rm=T)  #6 bins
act.bin.lims  #all 0s fall w/in bottom 50% of data

# Viz density distribs w/ bin limits
ggplot(filtered.df, aes(Activity.Count)) +
  geom_density(fill = "cadetblue") +
  geom_vline(data = data.frame(lims = act.bin.lims), aes(xintercept = lims), 
             linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))
  # facet_wrap(~ id)



# Assign bins to observations
dat_disc.list<- map(dat.list,
                       discrete_move_var,
                       lims = list(act.bin.lims),
                       varIn = c("Activity.Count"),
                       varOut = c("AC"))

#Viz discretized distrib
dat_disc.list2<- map(dat_disc.list, subset, select = c(id, AC))
dat_disc.df<- dat_disc.list2 %>% 
  bind_rows() %>% 
  gather(key, value, -id)

param.prop<- dat_disc.df %>%
  group_by(key, value) %>%
  summarise(n=n()) %>%
  mutate(prop=n/nrow(dat_disc.df)) %>%
  ungroup()  #if don't ungroup after grouping, ggforce won't work

param.prop[1:6, "value"]<- ((diff(act.bin.lims)/2) + act.bin.lims[1:6])


ggplot(data = param.prop , aes(value, prop)) +
  geom_bar(stat = "identity", width = (diff(act.bin.lims)),
           fill = "lightblue", color = "black") +
  # facet_zoom(xlim = c(0,3)) +
  labs(x = "Activity Counts", y = "Proportion") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))





### Export Data ###
dat.out<- bind_rows(dat_disc.list)

# write.csv(dat.out, "Binned Armadillo Acceleration Data.csv", row.names = F)
