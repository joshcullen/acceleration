### Identify 'ecological seasons' using mixture model for clustering

library(tidyverse)
library(lubridate)
library(raster)
library(viridis)
library(bayesmove)


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
dat$month<- factor(dat$month, levels = month.abb[c(5:12,1)])



#NDWI
setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")

ndwi<- brick('GiantArm_ndwi_monthly.grd')
ndwi.df<- as.data.frame(ndwi, xy = T)
ndwi.df2<- pivot_longer(ndwi.df, cols = -c(x,y), names_to = "month", values_to = "ndwi")
ndwi.df2$month<- factor(ndwi.df2$month, levels = names(ndwi))



#NDVI

ndvi<- brick('GiantArm_ndvi_monthly.grd')
compareRaster(ndwi, ndvi)
ndvi.df<- as.data.frame(ndvi, xy = T)
ndvi.df2<- pivot_longer(ndvi.df, cols = -c(x,y), names_to = "month", values_to = "ndvi")
ndvi.df2$month<- factor(ndvi.df2$month, levels = names(ndvi))

setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration")




##############################
### Range Standardize Data ###
##############################

rng.std<- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

sl.ta<- dat %>% 
  dplyr::select(sl2, ta2, month) %>% 
  mutate(across(c('sl2','ta2'), rng.std, na.rm = T)) %>% 
  group_by(month) %>% 
  summarize(mean.sl = mean(sl2, na.rm = T), var.sl = var(sl2, na.rm = T),
            min.sl = min(sl2, na.rm = T), max.sl = max(sl2, na.rm = T),
            mean.ta = mean(ta2, na.rm = T), var.ta = var(ta2, na.rm = T),
            min.ta = min(ta2, na.rm = T), max.ta = max(ta2, na.rm = T)) %>% 
  ungroup()

env.covs<- cbind(ndvi = ndvi.df2$ndvi, ndwi.df2[,c('ndwi','month')]) %>% 
  mutate(across(c('ndvi','ndwi'), rng.std, na.rm = T)) %>% 
  group_by(month) %>% 
  summarize(mean.ndvi = mean(ndvi, na.rm = T), var.ndvi = var(ndvi, na.rm = T),
            min.ndvi = min(ndvi, na.rm = T), max.ndvi = max(ndvi, na.rm = T),
            mean.ndwi = mean(ndwi, na.rm = T), var.ndwi = var(ndwi, na.rm = T),
            min.ndwi = min(ndwi, na.rm = T), max.ndwi = max(ndwi, na.rm = T)) %>% 
  ungroup()


## Combine summarized datasets and discretize into bins with a binwidth of 0.1
dat.strms<- left_join(sl.ta, env.covs, by = 'month')

bin.lims<- seq(0, 1, 0.1)
dat.strms2<- discrete_move_var(dat.strms[,-1],
                               lims = rep(list(bin.lims), ncol(dat.strms[,-1])),
                               varIn = names(dat.strms[,-1]),
                               varOut = paste0(names(dat.strms[,-1]), 1))
dat.strms2<- dat.strms2[,grepl(1, names(dat.strms2))] %>%  #remove original continuous data
  data.frame()



#################
### Run model ###
#################

set.seed(123)

# Define model params
alpha=0.1  #prior
ngibbs=20000  #number of Gibbs sampler iterations
nburn=ngibbs/2  #number of burn-in iterations
nmaxclust=4  #number of maximum possible states (clusters) present

# Run model
dat.res<- cluster_obs(dat=dat.strms2[,c('mean.ndvi1','min.ndvi1')],
                      alpha=alpha, ngibbs=ngibbs, nmaxclust=nmaxclust,
                      nburn=nburn)
# takes 33 s to run 20,000 iterations

# Inspect traceplot of log-likelihood
post.seq=(nburn + 1):ngibbs
plot(dat.res$loglikel, type = "l")
plot(dat.res$loglikel[post.seq], type = "l")

# Inspect likely number of groups (seasons)
theta<- dat.res$theta[post.seq,]
colnames(theta)<- 1:ncol(theta)
theta1<- colMeans(theta)
theta1<- sort(theta1, decreasing = TRUE)
cumsum(theta1)

# Inspect state assignments by each method
z.post<- as.matrix(dat.res$z.posterior)
z.post2<- t(apply(z.post, 1, function(x) x/sum(x)))
thresh<- 0.75
z.post3<- apply(z.post2, 1, function(x) ifelse(max(x) > thresh, which(x > thresh), NA))
z.post4<- apply(z.post2, 1, function(x) which.max(x))
dat.res$z.MAP

## Jan - July appears to represent flooding; Aug - Dec represents dry




##################################
### Aggregate into new seasons ###
##################################

season.ind<- ifelse(names(ndvi) %in% month.abb[1:7], "Flood", "Dry")
ndvi.season<- stackApply(ndvi, season.ind, fun = mean)
names(ndvi.season)<- c("Flood","Dry")

breaks<- seq(0, 1, by=0.001)
cols<- colorRampPalette(RColorBrewer::brewer.pal(9, 'RdYlGn'))(length(breaks)-1)
rasterVis::levelplot(ndvi.season, at=breaks, col.regions=cols, main="NDVI")



ndwi.season<- stackApply(ndwi, season.ind, fun = mean)
names(ndwi.season)<- c("Flood","Dry")

breaks<- seq(-1, 1, by=0.0001)
cols<- colorRampPalette(RColorBrewer::brewer.pal(9, 'RdBu'))(length(breaks)-1)
rasterVis::levelplot(ndwi.season, at=breaks, col.regions=cols, main="NDWI")





############################################################
### Extract Values of Continuous Variables by LULC Class ###
############################################################

#LULC
flood.lulc<- raster('cheiann_UTM1.tif')
names(flood.lulc)<- 'lulc'

#crop raster closer to armadillo relocs
flood.lulc<- crop(flood.lulc, ndvi.season)
flood.lulc.df<- as.data.frame(flood.lulc, xy = TRUE)


dry.lulc<- raster('secann_UTM1.tif')
names(dry.lulc)<- 'lulc'
dry.lulc<- crop(dry.lulc, flood.lulc)
dry.lulc.df<- as.data.frame(dry.lulc, xy = TRUE)

# Plot seasonal LULC
lulc<- rbind(dry.lulc.df, flood.lulc.df)
lulc$season<- rep(c("Dry","Wet"), each = nrow(dry.lulc.df))


ggplot() +
  geom_raster(data = lulc, aes(x, y, fill = factor(lulc))) +
  scale_fill_manual("", values = c("darkgreen","burlywood4","darkolivegreen3","lightskyblue1"),
                    na.value = "transparent",
                    labels = c("Forest", "Closed Savanna", "Open Savanna", "Floodable","")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold")) +
  facet_wrap(~ season)



flood.lulc.list<- list()
for (i in 1:max(getValues(flood.lulc), na.rm = T)) {
  
  ind<- which(getValues(flood.lulc) == i)
  
  flood.lulc.list[[i]]<- data.frame(extract(ndvi.season$Flood, ind), class = i)
}

flood.lulc.contvals<- bind_rows(flood.lulc.list)
flood.lulc.contvals$class<- factor(flood.lulc.contvals$class)
levels(flood.lulc.contvals$class)<- c("Forest", "Closed Savanna", "Open Savanna", "Floodable")
flood.lulc.contvals<- flood.lulc.contvals %>% 
  pivot_longer(cols = -class, names_to = "covar", values_to = "value")

ggplot(data=flood.lulc.contvals, aes(x=value, fill = class)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d("") +
  labs(x = "Value", y = "Density") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))





dry.lulc.list<- list()
for (i in 1:max(getValues(dry.lulc), na.rm = T)) {
  
  ind<- which(getValues(dry.lulc) == i)
  
  dry.lulc.list[[i]]<- data.frame(extract(ndvi.season$Dry, ind), class = i)
}

dry.lulc.contvals<- bind_rows(dry.lulc.list)
dry.lulc.contvals$class<- factor(dry.lulc.contvals$class)
levels(dry.lulc.contvals$class)<- c("Forest", "Closed Savanna", "Open Savanna", "Floodable")
dry.lulc.contvals<- dry.lulc.contvals %>% 
  pivot_longer(cols = -class, names_to = "covar", values_to = "value")

ggplot(data=dry.lulc.contvals, aes(x=value, fill = class)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d("") +
  labs(x = "Value", y = "Density") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))
