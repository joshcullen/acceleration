### Compare camera trap observations in front of burrows to instances of Slow-Turn behaviors estimated by mixture model


library(tidyverse)
library(lubridate)

#################
### Load Data ###
#################

dat<- read.csv("Giant Armadillo state estimates.csv", as.is = T)
dat<-  dat %>% 
  rename(x = easting, y = northing) %>%
  mutate(across(c('z.map','z.post.thresh','z.post.max'), factor,
                levels = c("Slow-Turn","Slow-Unif","Exploratory","Transit","Unclassified"))
  )
dat$date<- as_datetime(dat$date, tz = "UTC")  #if not already in POSIX format


cam.trap<- read.csv("cameratrap_behav.csv", as.is = T)
cam.trap<- cam.trap %>% 
  rename(date1 = date) %>% 
  mutate(across('id', tolower))
cam.trap$id<- str_replace(cam.trap$id, "mazebotti", "mazeboti")
cam.trap$date.begin<- as_datetime(paste(cam.trap$date1, cam.trap$time.begin), tz = "UTC")
cam.trap$date.end<- as_datetime(paste(cam.trap$date1, cam.trap$time.end), tz = "UTC")
cam.trap$date.end[12]<- cam.trap$date.end[12] + days(1)  #need to progress day for this obs



################################################################
### Identify state estimates that fall on or near camera obs ###
################################################################

id1<- unique(dat$id)
dat.filt<- list()

for (i in 1:length(unique(dat$id))) {
  
  tmp.dat<- dat %>% 
    filter(id == id1[i])
  tmp.cam<- cam.trap %>% 
    filter(id == id1[i])
  
  ind<- numeric()
  for (j in 1:nrow(tmp.cam)) {
    ind<- c(ind,
            which(tmp.dat$date >= (tmp.cam$date.begin[j] - minutes(3)) & 
            tmp.dat$date <= (tmp.cam$date.end[j] + minutes(3)))
    )
  }
  
  dat.filt[[i]]<- tmp.dat[ind,]
}

dat.filt<- bind_rows(dat.filt)

#only 9 matches found; 8 are Unclassified and 1 is Exploratory