
#### Pre-process armadillo data for bayesmove behavior model ####
## Last revised 2022-09-29 by Josh Cullen

library(bayesmove)
library(tidyverse)
library(lubridate)
library(patchwork)


### Import Armadillo GPS Data ###

files<- list.files(path = "Data_raw/", pattern = "*.csv", full.names = TRUE)
files<- files[!grepl(pattern = "accel", files)]
nfiles <- length(files)

#get gps data
fim <- numeric()
for (i in 1:nfiles){
  #import the data
  dat <- read.csv(files[i], as.is=T)
  
  #eliminate NA's
  cond <- !is.na(dat$GPS.UTM.Northing)
  dat1 <- dat[cond,]
  
  #just keep high quality GPS data (satellite count>X)
  print(table(dat1[,c('GPS.Fix.Attempt','GPS.Satellite.Count')]))
  cond <- dat1$GPS.Satellite.Count > 4
  dat2 <- dat1[cond,]
  
  #things to keep
  nomes <- c('id','GPS.Fix.Time','GPS.UTM.Northing','GPS.UTM.Easting')
  dat4 <- dat2[,nomes]
  
  #store into a single file
  fim <- rbind(fim,dat4)
}



### Import Armadillo Acceleration Data ###

dat <- read.csv('Data_raw/accel_data_200720.csv', as.is=T)
table(dat$Predeployment.Data)
table(dat$bout)
table(dat$act)

boxplot(Activity.Count~act, data=dat)
boxplot(Activity.Count~bout, data=dat)
boxplot(Activity.Count~hora, data=dat)

#just predeployment=No
cond <- dat$Predeployment.Data == 'No'; mean(cond)
nomes <- c('id','Acquisition.Time','Activity.Count')
accel <- dat[cond,nomes]





### Further filter and merge GPS and acceleration data ###

dat<- fim

#for each individual, calculate SL, angle, and time interval
uni <- unique(dat$id)
nuni <- length(uni)

fim2 <- numeric()
for (i in 1:nuni){
  #get accel data for individual i
  cond <- accel$id == uni[i]
  accel1 <- accel[cond,]
  
  #get GPS data for individual i  
  cond <- dat$id == uni[i]
  dat1 <- dat[cond,]
  nobs <- nrow(dat1)
  
  #get time
  dat1$time1 <- as.POSIXct(strptime(dat1$GPS.Fix.Time, 
                                 format = "%Y.%m.%d %H:%M:%S", 
                                 tz = "UTC"))
  accel1$time1 <- as.POSIXct(strptime(accel1$Acquisition.Time, 
                                   format = "%Y.%m.%d %H:%M:%S", 
                                   tz = "UTC"))
  
  #order the data chronologically
  ind <- order(dat1$time1)
  unique(ind-c(1:length(ind)))
  dat2 <- dat1[ind,]
  
  fim <- matrix(NA,nobs,10)
  for (j in 3:nobs){
    #get activity count within less than 1 min of end time of SL
    dt1 <- dat2$time1[j] - accel1$time1
    attr1 <- attributes(dt1)$units
    ind <- numeric()
    if (attr1 == 'secs') ind <- which(abs(dt1) <= 60)
    if (attr1 == 'min') ind <- which(abs(dt1) <= 1)
    count2 <- NA
    if (length(ind) == 1) count2 <- accel1$Activity.Count[ind]
    if (length(ind) > 1) count2 <- accel1$Activity.Count[sample(ind,size=1)]
    
    #get time
    time.int1 <- dat2$time1[j-1] - dat2$time1[j-2]
    time1 <- as.numeric(time.int1)
    unit1 <- attributes(time.int1)$units
    time.int2 <- dat2$time1[j] - dat2$time1[j-1]
    time2 <- as.numeric(time.int2)
    unit2 <- attributes(time.int2)$units
    
    #get dx's and dy's
    dx1 <- (dat2$GPS.UTM.Easting[j-1] - dat1$GPS.UTM.Easting[j-2])
    dy1 <- (dat2$GPS.UTM.Northing[j-1] - dat2$GPS.UTM.Northing[j-2])
    dx2 <- (dat2$GPS.UTM.Easting[j] - dat1$GPS.UTM.Easting[j-1])
    dy2 <- (dat2$GPS.UTM.Northing[j] - dat2$GPS.UTM.Northing[j-1])
    
    #get step length
    sl2 <- sqrt((dx2^2) + (dy2^2))
    
    #get angles 1 and 2
    angle1 <- atan2(dy1,dx1)
    angle2 <- atan2(dy2,dx2)
    tmp <- angle2 - angle1
    
    #wrap this around if angle is bigger than pi or smaller than -pi
    tmp <- ifelse(tmp > pi, tmp - 2*pi, tmp)
    TA2 <- ifelse(tmp < -pi, tmp + 2*pi, tmp)
    
    #store results
    fim[j,] <- c(sl2,TA2,time1,unit1,time2,unit2,count2,
              dat2$GPS.UTM.Northing[j],dat2$GPS.UTM.Easting[j],
              dat2$time1[j])
  }
  colnames(fim) <- c('sl2','ta2','time1','unit1','time2','unit2','act.count','northing','easting','date_time')
  fim2 <- rbind(fim2, cbind(uni[i],fim))
}
colnames(fim2)[1] <- 'id'
mean(is.na(fim2[,'act.count']))


##-------------------------
dat<- data.frame(fim2)
dat[,c(2:4,6,8:11)] <- apply(dat[,c(2:4,6,8:11)], 2, as.numeric)

table(dat$unit1)
table(dat$unit2)
cond <- dat$unit1 == 'hours'
hist(dat$time1[cond])
cond <- dat$unit1 == 'secs'
hist(dat$time1[cond])

#just keep those observations in which the difference in fixes is measured in minutes
cond <- !is.na(dat$unit2) & dat$unit2 == 'mins'
dat1 <- dat[cond,]; mean(cond, na.rm=T)

hist(dat1$time2)
rtime2 <- round(dat1$time2)
table(rtime2)

#just keep measurements between 6 and 8 minutes
cond <- !is.na(dat1$time2) & dat1$time2 >= 6 & dat1$time2 <= 8
dat2 <- dat1[cond,]; mean(cond)

#if time1 is not between 6 and 8, make TA=NA
cond <- !is.na(dat2$time1) & dat2$time1 >= 6 & dat2$time1 <= 8
dat2$ta2[!cond] <- NA

mean(is.na(dat2$ta2)); mean(is.na(dat2$sl2))

#eliminate extreme SL's
hist(dat2$sl2)
qtt <- quantile(dat2$sl2,0.99) #cut-off of 800
hist(dat2$sl2[dat2$sl2 < qtt])
cond <- dat2$sl2 < 800
dat3 <- dat2[cond,]

#normalize SL by dividing by time
dat3$sl2 <- dat3$sl2 / dat3$time2


#if act.count==0, make it equal to NA (not sure if mal-functioning or not). For example, might have fallen to the ground
cond <- !is.na(dat3$act.count) & dat3$act.count == 0
dat3$act.count[cond] <- NA


dat3<- dat3[,c('id','sl2','ta2','act.count','northing','easting','date_time')]
dat3<- dat3 %>% 
  rename(date = date_time) %>% 
  mutate_at("date", as_datetime, tz = "UTC")




# Viz map of tracks
ggplot(data = dat3, aes(easting, northing)) +
  geom_point(aes(color = id), alpha = 0.6) +
  geom_path() + 
  theme_bw() +
  facet_wrap(~id, scales = "free")

ggplot(data = dat3, aes(easting, northing)) +
  geom_path(aes(group = id, color = id)) + 
  geom_point(aes(color = id), alpha = 0.4) +
  scale_color_brewer("", palette = "Dark2") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.grid = element_blank())




# Viz distribution and time series of activity counts
ggplot(dat3, aes(act.count)) +
  geom_density(aes(fill = id)) +
  theme_bw() +
  facet_wrap(~ id)

ggplot() +
  geom_path(data = dat3, aes(date, act.count, color = id)) +
  theme_bw() +
  facet_wrap(~id, scales = "free_x")


# Viz distribution and time series of step lengths
ggplot(dat3, aes(sl2)) +
  geom_density(aes(fill = id)) +
  theme_bw() +
  facet_wrap(~ id)

ggplot() +
  geom_path(data = dat3, aes(date, sl2, color = id)) +
  theme_bw() +
  facet_wrap(~id, scales = "free_x")


#Viz distribution of turning angles
ggplot(dat3, aes(ta2)) +
  geom_density(aes(fill = id)) +
  theme_bw() +
  facet_wrap(~ id)
#looks much better

# Viz time series of turning angles
ggplot() +
  geom_path(data = dat3, aes(date, ta2, color = id)) +
  theme_bw() +
  facet_wrap(~id, scales = "free")



## Originally started w/ 140166 activity count obs and 18433 GPS obs, which was filtered to retain a total of 9948 obs



##########################################################
### Discretize data streams for Bayesian mixture model ###
##########################################################

#discretize act.count
hist(dat3$act.count)
act.brk <- seq(from=0, to=300, by=50)
range(dat3$act.count, na.rm=T)
act.brk[length(act.brk)] <- 300.1  #to include observation if exactly equal to 300
tmp <- cut(dat3$act.count, breaks=act.brk)
dat3$act.cat <- as.numeric(tmp)
sum(!is.na(dat3$act.count) & is.na(dat3$act.cat))
boxplot(act.count~act.cat, data=dat3)

#discretize SL
hist(dat3$sl2)
mean(dat3$sl2 > 60)
sl.brk <- c(seq(from=0, to=60, by=10), max(dat3$sl2) + 1)  #to include observation if exactly equal to max speed
dat3$sl.cat <- NA
for (i in 2:length(sl.brk)){
  cond <- dat3$sl2 >= sl.brk[i-1] & dat3$sl2 < sl.brk[i]
  dat3$sl.cat[cond] <- i-1
}
table(dat3$sl.cat)

#discretize TA
hist(dat3$ta2)
range(dat3$ta2, na.rm=T)
ta.brk <- seq(from=-pi, to=pi, length.out=11)
ta.brk[length(ta.brk)] <- 3.142  #to include observation if exactly equal to pi
tmp <- cut(dat3$ta2, breaks=ta.brk)
dat3$ta.cat <- as.numeric(tmp)
table(dat3$ta.cat)
sum(!is.na(dat3$ta2) & is.na(dat3$ta.cat))
boxplot(ta2~ta.cat, data=dat3)









# Viz density distribs w/ bin limits
p.act <- ggplot(dat3, aes(act.count)) +
  geom_density(fill = "cadetblue") +
  geom_vline(data = data.frame(lims = act.brk), aes(xintercept = lims), 
             linetype = "dashed") +
  labs(x = "Activity count", y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(30, 5.5, 5.5, 5.5), "pt"))


p.speed <- ggplot(dat3, aes(sl2/60)) +  #plot in m/s instead of m/min
  geom_density(fill = "lightblue") +
  geom_vline(data = data.frame(lims = sl.brk), aes(xintercept = lims/60),  #plot in m/s instead of m/min
             linetype = "dashed") +
  labs(x = "Speed (m/s)", y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


p.ta <- ggplot(dat3, aes(ta2)) +
  geom_density(fill = "indianred") +
  geom_vline(data = data.frame(lims = ta.brk), aes(xintercept = lims), 
             linetype = "dashed") +
  labs(x = "Turning angle (radians)", y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


## Create composite plot
p.act / p.speed / p.ta + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag.position = c(0.08, 1), plot.tag = element_text(size = 18, hjust = 0, vjust = -0.4,
                                                                face = 'bold'))

# ggsave("Figures/Figure 3.tiff", width = 5, height = 7, units = "in", dpi = 400)



#Viz discretized distrib
dat_disc.df<- dat3 %>% 
  dplyr::select(act.cat, sl.cat, ta.cat) %>% 
  pivot_longer(cols = c(act.cat, sl.cat, ta.cat), names_to = "var", values_to = "value")

param.prop<- dat_disc.df %>%
  group_by(var, value) %>%
  summarise(n=n()) %>%
  drop_na() %>% 
  mutate(prop=n/sum(n)) %>%
  ungroup() %>%   #if don't ungroup after grouping, ggforce won't work
  data.frame()

param.prop[1:6, "value"]<- ((diff(act.brk)/2) + act.brk[1:6])
param.prop[7:13, "value"]<- ((diff(sl.brk)/2) + sl.brk[1:7])
param.prop[14:23, "value"]<- (diff(ta.brk)/2) + ta.brk[1:10]

ggplot(data = param.prop %>% filter(var == "act.cat"), aes(value, prop)) +
  geom_bar(stat = "identity", width = diff(act.brk),
           fill = "cadetblue", color = "black") +
  labs(x = "Activity Count", y = "Proportion") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))

ggplot(data = param.prop %>% filter(var == "sl.cat"), aes(value/60, prop)) +
  geom_bar(stat = "identity", width = diff(sl.brk/60),
           fill = "lightblue", color = "black") +
  labs(x = "Step Length (m/min)", y = "Proportion") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))

ggplot(data = param.prop %>% filter(var == "ta.cat"), aes(value, prop)) +
  geom_bar(stat = "identity", width = diff(ta.brk),
           fill = "indianred", color = "black") +
  labs(x = "Turning Angle (rad)", y = "Proportion") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))





### Export Data ###

# write.csv(dat3, "Data_processed/Giant Armadillo Binned Data.csv", row.names = F)
