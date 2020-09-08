
#### Pre-process armadillo data for bayesmove behavior model ####

library(bayesmove)
library(tidyverse)
library(lubridate)


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
  drop_na(Activity.Count) %>%  #remove NA rows for data of interest
  mutate_at("Acquisition.Time",
            ~as.POSIXct(strptime(., format = "%Y.%m.%d %H:%M:%S"))) %>%  #make date POSIXct format
  mutate_at("Acquisition.Start.Time",
            ~as.POSIXct(strptime(., format = "%Y.%m.%d %H:%M:%S")))  #make date POSIXct format

# Add dt col (units as min)
dat3<- dat2 %>% 
  group_split(id) %>% 
  map(., arrange, Acquisition.Time) %>%  #reorder rows by Acquisition.Time
  map(., distinct, Acquisition.Time, .keep_all = T) %>% #remove duplicate rows
  map(., ~mutate(., dt = difftime(Acquisition.Time, Acquisition.Start.Time, units = "min"))) %>% 
  bind_rows()
table(dat3$dt)  #all but 5 of 140166 (0.004%) are at 5 min interval


# Filter where dt == 5
dat.list<- df_to_list(dat3, "id")
dat.list<- filter_time(dat.list, int = 5)


# Viz distribution of Activity.Count
filtered.df<- bind_rows(dat.list)

ggplot(filtered.df, aes(Activity.Count)) +
  geom_density(aes(fill = id)) +
  theme_bw() +
  facet_wrap(~ id)

# All IDs show highly right-skewed distributions; use quantiles to define bin limits





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
