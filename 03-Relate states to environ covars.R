### Relate state estimates to environmental covariates

library(tidyverse)
library(lubridate)
library(raster)
library(viridis)
library(ks)
library(cowplot)

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
dat$month<- factor(dat$month, levels = month.abb[c(5:12,1)])
dat$season<- ifelse(dat$month %in% c(month.abb[3:5]), "Fall",
                    ifelse(dat$month %in% c(month.abb[6:8]), "Winter",
                           ifelse(dat$month %in% c(month.abb[9:11]), "Spring", "Summer")))
dat$season<- factor(dat$season, levels = c("Fall","Winter","Spring","Summer"))
dat$season2<- ifelse(dat$month %in% month.abb[c(10:12,1:3)], "Rainy", "Dry")



### Extract all available environ covars and check for correlations before deciding which to include

#LULC
flood.lulc<- raster('cheiann_UTM1.tif')
names(flood.lulc)<- 'lulc'

#crop raster closer to armadillo relocs
flood.lulc<- crop(flood.lulc, extent(dat %>% 
                           summarize(xmin = min(x) - 3000,
                                     xmax = max(x) + 3000,
                                     ymin = min(y) - 3000,
                                     ymax = max(y) + 3000) %>% 
                           unlist()))
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


#NDWI
setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")

ndwi<- brick('GiantArm_ndwi_season.grd')
ndwi<- resample(ndwi, flood.lulc, method = "bilinear")
compareRaster(flood.lulc, ndwi)

ndwi.df<- as.data.frame(ndwi, xy = T)
ndwi.df2<- pivot_longer(ndwi.df, cols = -c(x,y), names_to = "season", values_to = "ndwi")
ndwi.df2$season<- factor(ndwi.df2$season, levels = names(ndwi))



#NDVI

ndvi<- brick('GiantArm_ndvi_season.grd')
ndvi<- resample(ndvi, flood.lulc, method = "bilinear")
compareRaster(flood.lulc, ndvi)

ndvi.df<- as.data.frame(ndvi, xy = T)
ndvi.df2<- pivot_longer(ndvi.df, cols = -c(x,y), names_to = "season", values_to = "ndvi")
ndvi.df2$season<- factor(ndvi.df2$season, levels = names(ndvi))

setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration")



#Tasseled Cap Brightness

bright<- brick('GiantArm_tcbright_season.grd')
bright<- resample(bright, flood.lulc, method = "bilinear")
compareRaster(flood.lulc, bright)

bright.df<- as.data.frame(bright, xy = T)
bright.df2<- pivot_longer(bright.df, cols = -c(x,y), names_to = "season", values_to = "bright")
bright.df2$season<- factor(bright.df2$season, levels = names(bright))



#Tasseled Cap Greenness

green<- brick('GiantArm_tcgreen_season.grd')
green<- resample(green, flood.lulc, method = "bilinear")
compareRaster(flood.lulc, green)

green.df<- as.data.frame(green, xy = T)
green.df2<- pivot_longer(green.df, cols = -c(x,y), names_to = "season", values_to = "green")
green.df2$season<- factor(green.df2$season, levels = names(green))



#Tasseled Cap Wetness

wet<- brick('GiantArm_tcwet_season.grd')
wet<- resample(wet, flood.lulc, method = "bilinear")
compareRaster(flood.lulc, wet)

wet.df<- as.data.frame(wet, xy = T)
wet.df2<- pivot_longer(wet.df, cols = -c(x,y), names_to = "season", values_to = "wet")
wet.df2$season<- factor(wet.df2$season, levels = names(wet))



#Elevation
dem<- raster('giantarm_dem.tif')
names(dem)<- 'elev'
dem<- resample(dem, flood.lulc, method = "bilinear")
compareRaster(flood.lulc, dem)
dem.df<- as.data.frame(dem, xy = TRUE)



## Viz behavioral states over (Dry) LULC map

#by z.post.thresh
ggplot(data = dat, aes(x, y)) +
  geom_raster(data = dry.lulc.df, aes(x, y, fill = factor(lulc))) +
  scale_fill_manual("", values = c("darkgreen","burlywood4","darkolivegreen3","lightskyblue1"),
                    na.value = "transparent",
                    labels = c("Forest", "Closed Savanna", "Open Savanna", "Floodable","")) +
  geom_path(aes(group = id), color = "grey20", size = 0.5) +
  geom_point(aes(color = z.post.thresh), size = 1, na.rm = T) +
  scale_color_manual("", values = c(viridis(n=4, option = 'inferno'), "grey")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))


ggplot(data = dat, aes(x, y)) +
  geom_raster(data = ndwi.df2, aes(x, y, fill = ndwi)) +
  scale_fill_viridis_c("", na.value = "transparent", limits = c(-1,1)) +
  geom_path(aes(group = id), color = "grey20", size = 0.5) +
  geom_point(aes(color = z.post.thresh), size = 1, na.rm = T) +
  scale_color_manual("", values = c(viridis(n=4, option = 'inferno'), "grey")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  facet_wrap(~ season)



######################################
### Fit and viz KUD for each state ###
######################################

## extract results
computeKDE = function(xy, prob=c("50%", "95%")) {
  
  kd <- kde(xy, compute.cont=TRUE)
  
  #define another function to be mapped
  get_contour = function(kd, prob) {
    contour <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                     z=estimate, levels=cont[prob])[[1]])
    as_tibble(contour) %>% 
      mutate(prob = prob)
  }

  dat_out <- map_dfr(prob, ~get_contour(kd, .)) %>% 
    group_by(prob) %>% 
    mutate(n_val = 1:n()) %>% 
    ungroup()

return(dat_out)
}

xy.slow_turn <- subset(dat, z.post.thresh == "Slow-Turn")[, c("x", "y")]
xy.slow_unif <- subset(dat, z.post.thresh == "Slow-Unif")[, c("x", "y")]
xy.expl <- subset(dat, z.post.thresh == "Exploratory")[, c("x", "y")]
xy.transit<- subset(dat, z.post.thresh == "Transit")[, c("x", "y")]

# subsample, for less autocorrelation
# xy.summer <- xy.summer[sample(1:nrow(xy.summer), 400), ]
# xy.winter <- xy.winter[sample(1:nrow(xy.winter), 400), ]

# compute kde
kde.slow_turn <- computeKDE(xy.slow_turn, prob = c(50,95))
kde.slow_turn$state<- "Slow-Turn"
kde.slow_unif <- computeKDE(xy.slow_unif, prob = c(50,95))
kde.slow_unif$state<- "Slow-Unif"
kde.expl <- computeKDE(xy.expl, prob = c(50,95))
kde.expl$state<- "Exploratory"
kde.transit <- computeKDE(xy.transit, prob = c(50,95))
kde.transit$state<- "Transit"

kde<- rbind(kde.slow_turn, kde.slow_unif, kde.expl, kde.transit)
kde$state<- factor(kde$state, levels = c("Slow-Turn", "Slow-Unif", "Exploratory",
                                         "Transit"))

#facet of all mapped KDE by state
ggplot() +
  geom_raster(data = flood.lulc.df, aes(x, y, fill = factor(lulc)), alpha = 0.5) +
  scale_fill_manual("", values = c("darkgreen","burlywood4","darkolivegreen3","lightskyblue1"),
                    na.value = "transparent",
                    labels = c("Forest", "Closed Savanna", "Open Savanna", "Floodable","")) +
  # geom_point(data = dat %>% 
  #              rename(state = z.post.thresh) %>% 
  #              filter(state != "Unclassified"),
  #            aes(x, y, color = id), size = 0.25, alpha = 0.75) +
  # scale_color_brewer("", palette = "Dark2") +
  geom_path(data = kde, aes(x, y, group = prob), colour = "black", size = 1) +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  coord_equal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold")) +
  facet_wrap(~ state, ncol = 2)




##########################################
### Extract environ covars for all obs ###
##########################################

dry.ind<- which(dat$season2 == "Dry")

# Proportions LULC w/in 30 m buffer (by season)

#Flood
flood.extr.lulc<- extract(flood.lulc, dat[,c('x','y')], buffer = 30)  #extract lulc w/in 30 m buffer
dry.extr.lulc<- extract(dry.lulc, dat[,c('x','y')], buffer = 30)

extr.lulc2<- map(flood.extr.lulc, ~{prop.table(table(.))}) %>%  #convert to proportions lulc
  map(., ~{
    mat<- matrix(0, 1, 4)
    mat[1, as.numeric(names(.x))]<- .x
    
    df<- data.frame(mat)
    df
    }) %>% 
  bind_rows()
dry.lulc2<- map(dry.extr.lulc, ~{prop.table(table(.))}) %>%  #convert to proportions lulc
  map(., ~{
    mat<- matrix(0, 1, 4)
    mat[1, as.numeric(names(.x))]<- .x
    
    df<- data.frame(mat)
    df
  }) %>% 
  bind_rows()

extr.lulc2[dry.ind,]<- dry.lulc2[dry.ind,]
names(extr.lulc2)<- c("Forest", "Closed_Savanna", "Open_Savanna", "Floodable")




# Mean elevation w/in 30 m buffer
extr.dem<- extract(dem, dat[,c('x','y')], buffer = 30, fun = mean)


# Mean NDWI, NDVI, Brightness, Greenness, and Wetness w/in 30 m buffer
seasons<- names(ndwi)
dyn.covars<- list(ndwi, ndvi, bright, green, wet)
names(dyn.covars)<- c("ndwi", "ndvi", "bright", "green", "wet")
extr.dyn.covars<- matrix(NA, nrow(dat), length(dyn.covars))
for (i in 1:length(seasons)) {
  ind<- which(dat$season == seasons[i])
  tmp<- map(dyn.covars, function(x) x[[i]]) %>% 
    stack()
  
  extr.dyn.covars[ind,]<- extract(tmp, dat[ind, c('x','y')], buffer = 30, fun = mean)
}
extr.dyn.covars<- data.frame(extr.dyn.covars)
names(extr.dyn.covars)<- names(dyn.covars)


# Merge extracted covariates to dat
dat2<- cbind(dat, extr.lulc2, elev = extr.dem, extr.dyn.covars)


# Check correlations among covariates (remove vars if |corr| > 0.7)
PerformanceAnalytics::chart.Correlation(dat2[,17:26])
## based on high corrs, removing NDWI, NDVI, and Brightness
## due to low level of added information, also removing LULC and elevation

######################################################
### Exploratory Data Analysis of Behavioral States ###
######################################################

## Proportion of behavior by ID

df.ID<- dat2 %>% 
  group_by(id, z.post.thresh) %>% 
  summarize(n = n()) %>% 
  mutate(freq = n / sum(n))

ggplot(df.ID, aes(z.post.thresh, freq, fill = z.post.thresh)) +
  geom_col(color = "black") +
  scale_fill_manual("", values = c(viridis(n=4, option = 'inferno'), "grey"),
                    guide = guide_legend(reverse = TRUE)) +
  theme_bw() +
  labs(x = "", y = "Proportion of Observations") +
  facet_wrap(~ id) +
  coord_flip() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))



## Proportion of behavior by season

df.season<- dat2 %>% 
  group_by(season, z.post.thresh) %>% 
  summarize(n = n()) %>% 
  mutate(freq = n / sum(n))

ggplot(df.season, aes(z.post.thresh, freq, fill = z.post.thresh)) +
  geom_col(color = "black") +
  scale_fill_manual("", values = c(viridis(n=4, option = 'inferno'), "grey"),
                    guide = guide_legend(reverse = TRUE)) +
  theme_bw() +
  labs(x = "", y = "Proportion of Observations") +
  facet_wrap(~ season) +
  coord_flip() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))



## Proportion of behavior by time of night

df.TON<- dat2 %>% 
  filter(z.post.thresh != "Unclassified") %>% 
  # mutate(date = date - hours(4)) %>% 
  mutate(hour1 = hour(date)) %>% 
  mutate_at('hour1', factor, levels = c(12:23,0:11)) %>% 
  mutate_at('z.post.thresh', droplevels) %>% 
  group_by(hour1, z.post.thresh, .drop = F) %>% 
  tally() %>% 
  mutate(freq = n / sum(n))

ggplot(data = df.TON, aes(hour1, n, fill = z.post.thresh)) +
  annotation_raster(alpha("grey30", .5),
                    xmin = 7, xmax = 19,
                    ymin = -Inf, ymax = Inf) +
  geom_bar(stat = "identity",
           color = "black") +
  scale_fill_manual("", values = c(viridis(n=4, option = 'inferno'), "grey")) +
  theme_bw() +
  labs(x = "Hour", y = "Number of Observations") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10))

ggplot(data = df.TON, aes(hour1, freq, fill = z.post.thresh)) +
  annotation_raster(alpha("grey30", .5),
                    xmin = 7, xmax = 19,
                    ymin = -Inf, ymax = Inf) +
  geom_bar(stat = "identity",
           color = "black") +
  scale_fill_manual("", values = c(viridis(n=4, option = 'inferno'), "grey")) +
  theme_bw() +
  labs(x = "Hour", y = "Proportion of Observations") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold")) +
  facet_wrap(~ z.post.thresh)



## Proportion of LULC by behavioral state

#filter for only "pure" land class (100% w/in 30 m buffer)
df.lulc<- dat2 %>% 
  filter(Forest == 1 | Closed_Savanna == 1 | Open_Savanna == 1 | Floodable == 1) %>% 
  dplyr::select(z.post.thresh, Forest, Closed_Savanna, Open_Savanna, Floodable) %>% 
  pivot_longer(cols = -z.post.thresh, names_to = "lulc", values_to = "value") %>% 
  filter(value == 1) %>% 
  dplyr::select(-value) %>% 
  group_by(lulc, z.post.thresh) %>% 
  tally() %>% 
  filter(z.post.thresh != "Unclassified") %>% 
  mutate(prop = n/sum(n))

df.lulc$lulc<- factor(df.lulc$lulc, levels = c("Forest", "Closed_Savanna", "Open_Savanna",
                                               "Floodable"))

ggplot(df.lulc, aes(lulc, n, fill = z.post.thresh)) +
  geom_bar(stat = "identity",
           color = "black") +
  scale_fill_manual("", values = viridis(n=4, option = 'inferno')) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_x_discrete(labels = c("Forest", "Closed Savanna", "Open Savanna", "Floodable")) +
  labs(x = "", y = "Number of Observations") +
  coord_flip() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))


ggplot(df.lulc, aes(lulc, prop, fill = z.post.thresh)) +
  geom_bar(stat = "identity",
           color = "black") +
  scale_fill_manual("", values = viridis(n=4, option = 'inferno')) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_x_discrete(labels = c("Forest", "Closed Savanna", "Open Savanna", "Floodable")) +
  labs(x = "", y = "Proportion of Observations") +
  coord_flip() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) +
  facet_wrap(~ z.post.thresh)


## Elevation by behavioral state
ggplot(dat2, aes(x = elev, color = z.post.thresh)) +
  geom_density(alpha = 0.25, size = 1) +
  scale_color_manual("", values = c(viridis(n=4, option = 'inferno'), "grey")) +
  labs(x = "Mean Elevation (w/in 30 m buffer)", y = "Density") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))


## Greenness by behavioral state
ggplot(dat2, aes(x = green, color = z.post.thresh)) +
  geom_density(alpha = 0.25, size = 1) +
  scale_color_manual("", values = c(viridis(n=4, option = 'inferno'), "grey")) +
  labs(x = "Mean Greenness (w/in 30 m buffer)", y = "Density") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))


## Wetness by behavioral state
ggplot(dat2, aes(x = wet, color = z.post.thresh)) +
  geom_density(alpha = 0.25, size = 1) +
  scale_color_manual("", values = c(viridis(n=4, option = 'inferno'), "grey")) +
  labs(x = "Mean Wetness (w/in 30 m buffer)", y = "Density") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))




#######################################################
### Analyze relationships between states and covars ###
#######################################################

plot(dat2$z.post.thresh, dat2$Forest)
plot(dat2$z.post.thresh, dat2$Closed_Savanna)
plot(dat2$z.post.thresh, dat2$Open_Savanna)
plot(dat2$z.post.thresh, dat2$Floodable)
plot(dat2$z.post.thresh, dat2$elev)
plot(dat2$z.post.thresh, dat2$green)
plot(dat2$z.post.thresh, dat2$wet)

#remove observations w/ "Unclassified" state (based on z.post.thresh)
dat3<- dat2 %>% 
  filter(z.post.thresh != "Unclassified") %>% 
  mutate(across(.cols = c(Forest:wet), scale))




## Compare Slow-Turn vs Slow-Unif (Closed Savanna and blanca are reference classes)
dat.st.su<- dat3 %>% 
  filter(z.post.thresh == "Slow-Turn" | z.post.thresh == "Slow-Unif")
dat.st.su$state<- ifelse(dat.st.su$z.post.thresh == "Slow-Turn", 0, 1)
table(dat.st.su$id, dat.st.su$state)  #check n per combo
st.su.mod<- glm(state ~ green + wet + id,
                data = dat.st.su, family = binomial)
summary(st.su.mod)

# Slow-Unif as ref
## odds ratios and 95% CI
coeffs.st.su<- data.frame(exp(cbind(fit = coef(st.su.mod), confint(st.su.mod))))
coeffs.st.su<- coeffs.st.su[2:3,]
names(coeffs.st.su)[2:3]<- c("Lower","Upper")
coeffs.st.su$coef.names<- rownames(coeffs.st.su)
coeffs.st.su$coef.names<- factor(coeffs.st.su$coef.names,
                                 levels = unique(coeffs.st.su$coef.names))

p.st.su<- ggplot(data=coeffs.st.su, aes(x=coef.names, y=fit, ymin=Lower, ymax=Upper)) +
  geom_hline(yintercept = 1) +
  geom_errorbar(position = position_dodge(0.55), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.55), size=2) +
  ylim(c(0,2)) +
  annotate(geom = "text", x = 1.5, y = 1.5, label = "bold(Slow-Uniform)",
           parse = T, size = 4) +
  annotate(geom = "text", x = 1.5, y = 0.5, label = "bold(Slow-Turn)",
           parse = T, size = 4) +
  theme_bw() +
  scale_x_discrete(labels = c("Greenness","Wetness")) +
  coord_flip() +
  labs(x="", y="Odds Ratio") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        panel.grid = element_blank())




## Compare Slow-Turn vs Exploratory
dat.st.exp<- dat3 %>% 
  filter(z.post.thresh == "Slow-Turn" | z.post.thresh == "Exploratory")
dat.st.exp$state<- ifelse(dat.st.exp$z.post.thresh == "Slow-Turn", 0, 1)
table(dat.st.exp$id, dat.st.exp$state)  #check n per combo
st.exp.mod<- glm(state ~ green + wet + id,
                data = dat.st.exp, family = binomial)
summary(st.exp.mod)

# Exploratory as ref
## odds ratios and 95% CI
coeffs.st.exp<- data.frame(exp(cbind(fit = coef(st.exp.mod), confint(st.exp.mod))))
coeffs.st.exp<- coeffs.st.exp[2:3,]
names(coeffs.st.exp)[2:3]<- c("Lower","Upper")
coeffs.st.exp$coef.names<- rownames(coeffs.st.exp)
coeffs.st.exp$coef.names<- factor(coeffs.st.exp$coef.names,
                                 levels = unique(coeffs.st.exp$coef.names))

p.st.exp<- ggplot(data=coeffs.st.exp, aes(x=coef.names, y=fit, ymin=Lower, ymax=Upper)) +
  geom_hline(yintercept = 1) +
  geom_errorbar(position = position_dodge(0.55), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.55), size=2) +
  ylim(c(0.25,1.75)) +
  annotate(geom = "text", x = 1.5, y = 1.3, label = "bold(Exploratory)",
           parse = T, size = 4) +
  annotate(geom = "text", x = 1.5, y = 0.65, label = "bold(Slow-Turn)",
           parse = T, size = 4) +
  theme_bw() +
  scale_x_discrete(labels = c("Greenness","Wetness")) +
  coord_flip() +
  labs(x="", y="Odds Ratio") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        panel.grid = element_blank())



## Compare Slow-Turn vs Transit
dat.st.t<- dat3 %>% 
  filter(z.post.thresh == "Slow-Turn" | z.post.thresh == "Transit")
dat.st.t$state<- ifelse(dat.st.t$z.post.thresh == "Slow-Turn", 0, 1)
table(dat.st.t$id, dat.st.t$state)  #check n per combo
st.t.mod<- glm(state ~ green + wet + id,
                 data = dat.st.t, family = binomial)
summary(st.t.mod)

# Transit as ref
## odds ratios and 95% CI
coeffs.st.t<- data.frame(exp(cbind(fit = coef(st.t.mod), confint(st.t.mod))))
coeffs.st.t<- coeffs.st.t[2:3,]
names(coeffs.st.t)[2:3]<- c("Lower","Upper")
coeffs.st.t$coef.names<- rownames(coeffs.st.t)
coeffs.st.t$coef.names<- factor(coeffs.st.t$coef.names,
                                  levels = unique(coeffs.st.t$coef.names))

p.st.t<- ggplot(data=coeffs.st.t, aes(x=coef.names, y=fit, ymin=Lower, ymax=Upper)) +
  geom_hline(yintercept = 1) +
  geom_errorbar(position = position_dodge(0.55), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.55), size=2) +
  ylim(c(0,2)) +
  annotate(geom = "text", x = 1.5, y = 1.6, label = "bold(Transit)",
           parse = T, size = 4) +
  annotate(geom = "text", x = 1.5, y = 0.4, label = "bold(Slow-Turn)",
           parse = T, size = 4) +
  theme_bw() +
  scale_x_discrete(labels = c("Greenness","Wetness")) +
  coord_flip() +
  labs(x="", y="Odds Ratio") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        panel.grid = element_blank())



## Compare Slow-Unif vs Exploratory
dat.su.exp<- dat3 %>% 
  filter(z.post.thresh == "Slow-Unif" | z.post.thresh == "Exploratory")
dat.su.exp$state<- ifelse(dat.su.exp$z.post.thresh == "Slow-Unif", 0, 1)
table(dat.su.exp$id, dat.su.exp$state)  #check n per combo
su.exp.mod<- glm(state ~ green + wet + id,
                 data = dat.su.exp, family = binomial)
summary(su.exp.mod)

# Exploratory as ref
## odds ratios and 95% CI
coeffs.su.exp<- data.frame(exp(cbind(fit = coef(su.exp.mod), confint(su.exp.mod))))
coeffs.su.exp<- coeffs.su.exp[2:3,]
names(coeffs.su.exp)[2:3]<- c("Lower","Upper")
coeffs.su.exp$coef.names<- rownames(coeffs.su.exp)
coeffs.su.exp$coef.names<- factor(coeffs.su.exp$coef.names,
                                  levels = unique(coeffs.su.exp$coef.names))

p.su.exp<- ggplot(data=coeffs.su.exp, aes(x=coef.names, y=fit, ymin=Lower, ymax=Upper)) +
  geom_hline(yintercept = 1) +
  geom_errorbar(position = position_dodge(0.55), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.55), size=2) +
  ylim(c(0.5,1.5)) +
  annotate(geom = "text", x = 1.5, y = 1.3, label = "bold(Exploratory)",
           parse = T, size = 4) +
  annotate(geom = "text", x = 1.5, y = 0.7, label = "bold(Slow-Uniform)",
           parse = T, size = 4) +
  theme_bw() +
  scale_x_discrete(labels = c("Greenness","Wetness")) +
  coord_flip() +
  labs(x="", y="Odds Ratio") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        panel.grid = element_blank())




## Compare Slow-Unif vs Transit
dat.su.t<- dat3 %>% 
  filter(z.post.thresh == "Slow-Unif" | z.post.thresh == "Transit")
dat.su.t$state<- ifelse(dat.su.t$z.post.thresh == "Slow-Unif", 0, 1)
table(dat.su.t$id, dat.su.t$state)  #check n per combo
su.t.mod<- glm(state ~ green + wet + id,
                 data = dat.su.t, family = binomial)
summary(su.t.mod)

# Transit as ref
## odds ratios and 95% CI
coeffs.su.t<- data.frame(exp(cbind(fit = coef(su.t.mod), confint(su.t.mod))))
coeffs.su.t<- coeffs.su.t[2:3,]
names(coeffs.su.t)[2:3]<- c("Lower","Upper")
coeffs.su.t$coef.names<- rownames(coeffs.su.t)
coeffs.su.t$coef.names<- factor(coeffs.su.t$coef.names,
                                  levels = unique(coeffs.su.t$coef.names))

p.su.t<- ggplot(data=coeffs.su.t, aes(x=coef.names, y=fit, ymin=Lower, ymax=Upper)) +
  geom_hline(yintercept = 1) +
  geom_errorbar(position = position_dodge(0.55), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.55), size=2) +
  ylim(c(0,2)) +
  annotate(geom = "text", x = 1.5, y = 1.4, label = "bold(Transit)",
           parse = T, size = 4) +
  annotate(geom = "text", x = 1.5, y = 0.5, label = "bold(Slow-Uniform)",
           parse = T, size = 4) +
  theme_bw() +
  scale_x_discrete(labels = c("Greenness","Wetness")) +
  coord_flip() +
  labs(x="", y="Odds Ratio") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        panel.grid = element_blank())




## Compare Exploratory vs Transit
dat.exp.t<- dat3 %>% 
  filter(z.post.thresh == "Exploratory" | z.post.thresh == "Transit")
dat.exp.t$state<- ifelse(dat.exp.t$z.post.thresh == "Exploratory", 0, 1)
table(dat.exp.t$id, dat.exp.t$state)  #check n per combo
exp.t.mod<- glm(state ~ green + wet + id,
                 data = dat.exp.t, family = binomial)
summary(exp.t.mod)

# Transit as ref
## odds ratios and 95% CI
coeffs.exp.t<- data.frame(exp(cbind(fit = coef(exp.t.mod), confint(exp.t.mod))))
coeffs.exp.t<- coeffs.exp.t[2:3,]
names(coeffs.exp.t)[2:3]<- c("Lower","Upper")
coeffs.exp.t$coef.names<- rownames(coeffs.exp.t)
coeffs.exp.t$coef.names<- factor(coeffs.exp.t$coef.names,
                                  levels = unique(coeffs.exp.t$coef.names))

p.exp.t<- ggplot(data=coeffs.exp.t, aes(x=coef.names, y=fit, ymin=Lower, ymax=Upper)) +
  geom_hline(yintercept = 1) +
  geom_errorbar(position = position_dodge(0.55), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.55), size=2) +
  ylim(c(0.25,1.75)) +
  annotate(geom = "text", x = 1.5, y = 1.4, label = "bold(Transit)",
           parse = T, size = 4) +
  annotate(geom = "text", x = 1.5, y = 0.5, label = "bold(Exploratory)",
           parse = T, size = 4) +
  theme_bw() +
  scale_x_discrete(labels = c("Greenness","Wetness")) +
  coord_flip() +
  labs(x="", y="Odds Ratio") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        panel.grid = element_blank())



### Make composite plot of all relationships

plot_grid(p.st.su, p.st.exp, p.st.t,
          p.su.exp, p.su.t, p.exp.t,
          nrow = 2)




##########################################################
### Explore Relationships between Vegetation and Water ###
##########################################################

#Create 2D KDE of NDVI-NDWI and Greenness-Wetness to investigate relationships and how this could impact interpretation of results

ggplot(data = dat2 %>% filter(z.post.thresh != "Unclassified"), aes(ndvi,ndwi)) +
  geom_density2d_filled(contour_var = "ndensity") +
  theme_bw() +
  labs(x = "NDVI", y = "NDWI") +
  facet_grid(z.post.thresh ~ season) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"))

ggplot(data = dat2 %>% filter(z.post.thresh != "Unclassified"), aes(green,wet)) +
  geom_density2d_filled(contour_var = "ndensity") +
  theme_bw() +
  labs(x = "Greenness", y = "Wetness") +
  facet_grid(z.post.thresh ~ season) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"))



###############################################
### Calculate daily displacements by season ###
###############################################

# Define nocturnal activity periods; if not already converted, substract 4 hours from UTC time
dat2$date<- dat2$date - hours(4)
dat2<- nocturnal_period(dat2)

#NEED TO FIX NOCTURNAL_PERIOD()


############################################################
### Extract Values of Continuous Variables by LULC Class ###
############################################################

covs<- stack(green[[1]], wet[[1]], ndvi[[1]], ndwi[[1]], bright[[1]], dem)
names(covs)<- c("Greenness", "Wetness", "NDVI", "NDWI", "Brightness", "Elevation")
flood.lulc.list<- list()

for (i in 1:max(getValues(flood.lulc), na.rm = T)) {
  
  ind<- which(getValues(flood.lulc) == i)
  
  flood.lulc.list[[i]]<- data.frame(extract(covs, ind), class = i)
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
        legend.text = element_text(size = 10)) +
  facet_wrap(~ covar, scales = "free")

