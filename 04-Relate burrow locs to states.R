### Relate VHF located burrow positions to estimated behavioral states

library(tidyverse)
library(lubridate)

source('helper functions.R')

#################
### Load data ###
#################

dat<- read.csv("Giant Armadillo state estimates.csv", as.is = T)
dat<-  dat %>% 
  rename(x = easting, y = northing) %>% 
  mutate(across(c('z.map','z.post.thresh','z.post.max'), factor,
                levels = c("Slow-Turn","Slow-Unif","Exploratory","Transit","Unclassified"))
  )
dat$date<- as_datetime(dat$date, tz = "UTC")
dat$month<- month.abb[month(dat$date)]
dat<- dat %>% 
  bayesmove::df_to_list("id") %>% 
  map(nocturnal_period) %>% 
  bind_rows()


burrow.locs<- read.csv("burrow_locations.csv", as.is = T)
burrow.locs$Nome<- gsub("Tex ", "Tex", burrow.locs$Nome)
burrow.locs<- burrow.locs %>% 
  filter(Nome %in% str_to_title(unique(dat$id))) %>% 
  rename(id = Nome, date = Data) %>% 
  mutate_at('id', tolower)
burrow.locs$date<- as_date(strptime(burrow.locs$date, "%m/%d/%Y", tz = "UTC"))

#fix some issues w data
ind<- c(26:38,52:60,85:91,96:102,111:112,123:126,141:151,157:168,190:197)
burrow.locs[ind, c('Y','X')]<- burrow.locs[ind, c('X','Y')]

#order rows by date (per ID)
burrow.locs<- burrow.locs %>% 
  group_by(id) %>% 
  arrange(date, .by_group = TRUE)

#add coords for time diff < 5 days when same burrow used
burrow.locs2<- burrow.locs %>%  #calc time diff per ID
  bayesmove::df_to_list("id") %>% 
  map(., ~mutate(., d.diff = as.numeric(difftime(.$date, lag(.$date), units = "days")))) %>% 
  bind_rows()

for (i in 2:nrow(burrow.locs2)) {
  if (burrow.locs2$d.diff[i] > 1 & 
      burrow.locs2$d.diff[i] <= 5 & 
      burrow.locs2$X[i] == burrow.locs2$X[i-1] &
      burrow.locs2$Y[i] == burrow.locs2$Y[i-1]) {
    
    burrow.locs2$seq.ind[i]<- 1
  } else {
    burrow.locs2$seq.ind[i]<- 0
  }
}

### only 4 of 208 observations match criteria for expansion of date windows; fine just sticking w/ strict date matches



#####################################################################
### Measure distance to recently used burrow in relation to state ###
#####################################################################

id1<- unique(dat$id)
dat.filt<- list()

for (i in 1:length(unique(dat$id))) {
  
  tmp.dat<- dat %>% 
    filter(id == id1[i])
  tmp.burrow<- burrow.locs %>% 
    filter(id == id1[i])
  for (j in 1:nrow(tmp.burrow)) {
    
    ind<- which(tmp.burrow$date[j] == date(tmp.dat$date))
    cond<- tmp.dat[ind, "period"] %>% 
      unique()
     
    boop<- tmp.dat %>% 
      filter(period %in% cond)
    
    ind2<- which(tmp.dat$period %in% cond)
    tmp.dat[ind2,"dist2burr"]<- sqrt((boop$x - tmp.burrow$X[j])^2 + (boop$y - tmp.burrow$Y[j])^2)
  }
  
  dat.filt[[i]]<- tmp.dat
}

dat.filt<- bind_rows(dat.filt)
dat.filt2<- dat.filt %>% 
  filter(!is.na(dist2burr))

# What is avg distance from most recent burrow per state?
dat.filt2 %>% 
  group_by(z.post.thresh) %>% 
  summarize(dist = mean(dist2burr))


# Plot movements compared to burrow
ggplot() +
  geom_point(data = dat.filt %>% 
               filter(id == 'emanuel') %>% 
               filter(date(date) == "2019-09-29"), aes(x, y, color = z.post.thresh)) +
  geom_point(data = burrow.locs %>% 
               filter(id == 'emanuel') %>% 
               filter(date(date) == "2019-09-29"), aes(X, Y), fill = "blue", shape = 25, size = 3) +
  theme_bw() +
  ggspatial::annotation_scale(location = 'tl') +
  facet_wrap(~z.post.thresh)


# Plot distribution of distances by state
dat.filt3<- dat.filt2 %>% 
  left_join(dat.filt2 %>% 
              group_by(id) %>% 
              summarise(N=n())) %>%
  mutate(Label=paste0(id,' (N = ',N,')'))

ggplot(data = dat.filt3) +
  geom_density(aes(dist2burr, color = z.post.thresh)) +
  theme_bw() +
  xlim(0,3000) +
  facet_wrap(~Label, scales = "free_y")

