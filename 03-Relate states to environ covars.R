### Relate state estimates to environmental covariates

library(tidyverse)
library(lubridate)
library(raster)
library(viridis)
library(mgcv)
library(wesanderson)
library(rnaturalearth)
library(rnaturalearthhires)
library(cowplot)
library(sf)
library(ggspatial)
library(brms)
library(ggdist)

source('helper functions.R')


#################
### Load data ###
#################

dat<- read.csv("Giant Armadillo state estimates.csv", as.is = T)
dat<-  dat %>% 
  rename(x = easting, y = northing) %>% 
  mutate(across(c('z.map','z.post.thresh','z.post.max'), factor,
                levels = c("VE","Local Search","Exploratory","Transit","Unclassified"))
  )
dat$date<- as_datetime(dat$date, tz = "UTC")
dat$month<- month.abb[month(dat$date)]
dat$month<- factor(dat$month, levels = month.abb[c(5:12,1)])
dat$season<- ifelse(dat$month %in% month.abb[c(10:12,1:3)], "Wet", "Dry")


#Plot of all tracks together
pal2<- c(wes_palettes$Darjeeling1, wes_palettes$Darjeeling2[2:1])  #palette
dat$id<- str_to_title(dat$id)

ind<- dat %>% 
  bayesmove::prep_data(., coord.names = c("x", "y"), id = "id") %>% 
  summarize(large.step = which(step > 800)) %>%   #find large gaps in tracks (based on filtering step where SL > 800 m were removed; 99th quantile)
  unlist()

# Inset NA rows at these large gaps (SL > 800)
dat.sf<- dat
for (i in 1:length(ind)) {
  dat.sf<- add_row(dat.sf, id = dat.sf$id[ind[i]], .after = ind[i])
  ind<- ind+1
}
# dat.sf<- dat.sf %>% 
#   st_as_sf(., coords = c("x", "y"), crs = 32721) %>% 
#   group_by(id) %>% 
#   dplyr::summarize(do_union = F) %>%
#   sf::st_cast("LINESTRING")


#ESRI REST API for world satellite imagery
esri_world <- paste0('https://services.arcgisonline.com/arcgis/rest/services/',
                     'World_Imagery/MapServer/tile/${z}/${y}/${x}.jpeg')

# This version (using geom_spatial_path) gives a warning, but still provides the correct map (for now)
main.map<- ggplot() +
  annotation_map_tile(type = esri_world, zoom = 13, progress = "none") +
  geom_spatial_path(data = dat.sf, aes(x, y, group = id, color = id), size = 0.5, crs = 32721) +
  scale_color_manual("",
                     values = pal2,
                     labels = c("Blanca (n=2074)",
                                "Emanuel (n=1032)",
                                "Gala (n=2008)",
                                "Mafalda (n=1439)",
                                "Mazeboti (n=1883)",
                                "Sara (n=596)",
                                "Tex (n=916)")) +
  labs(x = "Easting", y = "Northing") +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_continuous(expand = c(0.2, 0.2)) +
  annotation_scale(location = "br", width_hint = 0.5, style = "ticks",
                   line_col = "white", text_col = "white") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 8),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(fill = NA, color = viridis(3, option = "mako")[2],
                                    size = 2)) +
  guides(color = guide_legend(override.aes = list(size = 1)))

#access global map layer for inset
SA<- ne_countries(scale = 'medium', continent = 'south america', returnclass = "sf")
SA.tr<- st_transform(SA, 32721)
# SA.tr<- spTransform(SA, CRS("EPSG:32721"))
# # SA.df<- fortify(SA.tr)
# SA.sf<- st_as_sf(SA.tr)

inset1<- ggplot() +
  geom_sf(data = SA.tr, fill = "grey75", color = "white", size = 0.25) +
  geom_rect(aes(xmin = min(dat$x) - 3e5,
                xmax = max(dat$x) + 3e5,
                ymin = min(dat$y) - 3e5,
                ymax = max(dat$y) + 3e5),
            color = viridis(3, option = "mako")[1], fill = NA, size = 0.5) +
  scale_x_continuous(expand = c(0,0), limits = c(-2.5e6, st_bbox(SA.tr)$xmax)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_void()

inset2<- ggplot() +
  geom_sf(data = SA.tr, fill = "grey75", color = "white", size = 0.25) +
  geom_rect(aes(xmin = min(dat$x),
                xmax = max(dat$x),
                ymin = min(dat$y),
                ymax = max(dat$y)),
            color = viridis(3, option = "mako")[2],
            fill = viridis(3, option = "mako")[2],
            size = 1) +
  scale_x_continuous(expand = c(0,0), limits = c(min(dat$x) - 3e5, max(dat$x) + 3e5)) +
  scale_y_continuous(expand = c(0,0), limits = c(min(dat$y) - 3e5, max(dat$y) + 3e5)) +
  theme_void() +
  theme(panel.border = element_rect(color = viridis(3, option = "mako")[1], fill = NA,
                                    size = 1))

ggdraw(main.map) + 
  draw_plot(inset1, x = 0.09, y = 0.6, width = 0.15, height = 0.3) +
  draw_plot(inset2, x = 0.24, y = 0.68, width = 0.17, height = 0.17)

# ggsave("Figure 1_new.png", width = 9, height = 5, units = "in", dpi = 330)





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
  theme(legend.position = "top",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_blank()) +
  facet_wrap(~ season)

# ggsave("Figure 2.png", width = 7, height = 4, units = "in", dpi = 330)



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

# ggsave("Giant Arm behav_lulc.png", width = 7, height = 5, units = "in", dpi = 330)


# Example plot of region subset
ggplot() +
  geom_raster(data = dry.lulc.df, aes(x, y, fill = factor(lulc))) +
  scale_fill_manual("", values = c("darkgreen","burlywood4","darkolivegreen3","lightskyblue1"),
                    na.value = "transparent",
                    labels = c("Forest", "Closed Savanna", "Open Savanna", "Floodable","")) +
  scale_x_continuous(expand = c(0,0), limits = c(626250,627200)) +
  scale_y_continuous(expand = c(0,0), limits = c(7875000,7876000)) +
  labs(x="Easting", y="Northing") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# ggsave("Giant Arm lulc_subset.png", width = 6, height = 4, units = "in", dpi = 330)




##########################################
### Extract environ covars for all obs ###
##########################################

dry.ind<- which(dat$season == "Dry")

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


# Merge extracted covariates to dat
dat2<- cbind(dat, extr.lulc2)





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
  mutate(date = date - hours(4)) %>%
  mutate(hour1 = hour(date)) %>% 
  mutate_at('hour1', factor, levels = c(12:23,0:11)) %>% 
  mutate_at('z.post.thresh', droplevels) %>% 
  group_by(hour1, z.post.thresh, .drop = F) %>% 
  tally() %>% 
  mutate(freq = n / sum(n))

p.no<- ggplot(data = df.TON, aes(hour1, n, fill = z.post.thresh)) +
  annotation_raster(alpha("grey30", .5),
                    xmin = 6, xmax = 18,
                    ymin = -Inf, ymax = Inf) +
  geom_bar(stat = "identity",
           color = "black") +
  scale_fill_manual("", values = c(viridis(n=4, option = 'inferno'), "grey")) +
  theme_bw() +
  labs(x = "Time of day (h)", y = "Number of Observations") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = "top",
        panel.grid.minor = element_blank())


p.prop<- ggplot(data = df.TON, aes(hour1, freq, fill = z.post.thresh)) +
  annotation_raster(alpha("grey30", .5),
                    xmin = 6, xmax = 18,
                    ymin = -Inf, ymax = Inf) +
  geom_bar(stat = "identity",
           color = "black") +
  scale_fill_manual("", values = c(viridis(n=4, option = 'inferno'), "grey"), guide = "none") +
  theme_bw() +
  labs(x = "Time of day (h)", y = "Proportion of Observations") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  facet_wrap(~ z.post.thresh)


plot_grid(p.no, p.prop, nrow = 2)
# ggsave("Figure 4.png", width = 12, height = 8, dpi = 330)


## Reclassify 'Time of Night' per activity period
# 1 = 1st hour
# 2 = last hour
# 3 = everything in between

dat2$hour<- hour(dat2$date)

foo<- dat2 %>% 
  bayesmove::df_to_list("id") %>% 
  map(nocturnal_period)


for (i in 1:length(foo)) {
  tmp<- foo[[i]]
  TON<- rep(NA, nrow(tmp))
  
  for (j in 1:n_distinct(tmp$period)) {
    ind<- which(tmp$period == j)
    tmp1<- tmp[ind,]
    
    hrs<- unique(tmp1$hour)
    TON[ind]<- ifelse(tmp1$hour == hrs[1], 1,
                      ifelse(tmp1$hour == hrs[length(hrs)], 2, 3))
    
  }
  
  foo[[i]]$TON<- TON
}

foo2<- bind_rows(foo)


## Calculate avg activity duration per night
test<- foo %>% 
  bind_rows %>% 
  group_by(id, period) %>% 
  slice(c(1,n())) %>% 
  summarize(dur = difftime(date[2], date[1], units = "hours"))

test %>% 
  filter(dur < 13) %>%  #need to filter out a few ~35 hr durations (issue w/ period assignment)
  ungroup() %>% 
  summarize(mean = mean(dur))

df.TON2<- foo2 %>% 
  filter(z.post.thresh != "Unclassified") %>% 
  mutate_at('TON', factor, levels = c(1,3,2)) %>% 
  mutate_at('z.post.thresh', droplevels) %>% 
  group_by(TON, z.post.thresh, .drop = F) %>% 
  tally() %>% 
  mutate(freq = n / sum(n))

ggplot(data = df.TON2, aes(TON, n, fill = z.post.thresh)) +
  geom_bar(stat = "identity",
           color = "black") +
  scale_fill_manual("", values = c(viridis(n=4, option = 'inferno'), "grey")) +
  theme_bw() +
  labs(x = "Relative Time of Night", y = "Number of Observations") +
  scale_x_discrete(labels = c("Start", "Mid", "Late")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10))


ggplot(data = df.TON2, aes(TON, freq, fill = z.post.thresh)) +
  geom_bar(stat = "identity",
           color = "black") +
  scale_fill_manual("", values = c(viridis(n=4, option = 'inferno'), "grey")) +
  theme_bw() +
  labs(x = "Relative Time of Night", y = "Proportion of Observations") +
  scale_x_discrete(labels = c("Start", "Mid", "Late")) +
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


# Proportion behavior vs proportion LULC

df.state.forest<- dat2 %>% 
  filter(z.post.thresh != 'Unclassified') %>% 
  mutate(across(c(Forest:Floodable), factor)) %>% 
  group_by(Forest, z.post.thresh) %>% 
  summarize(n = n()) %>% 
  group_by(Forest) %>% 
  mutate(prop = n/sum(n))

ggplot(df.state.forest, aes(as.numeric(as.character(Forest)), prop, fill = z.post.thresh)) +
  geom_area() +
  scale_fill_viridis_d(option = 'inferno') +
  labs(x = 'Proportion of Forest', y = 'Proportion of States') +
  theme_bw()


df.state.closed_savanna<- dat2 %>% 
  filter(z.post.thresh != 'Unclassified') %>% 
  mutate(across(c(Forest:Floodable), factor)) %>% 
  group_by(Closed_Savanna, z.post.thresh) %>% 
  summarize(n = n()) %>% 
  group_by(Closed_Savanna) %>% 
  mutate(prop = n/sum(n))

ggplot(df.state.closed_savanna, aes(as.numeric(as.character(Closed_Savanna)), prop,
                                    fill = z.post.thresh)) +
  geom_area() +
  scale_fill_viridis_d(option = 'inferno') +
  labs(x = 'Proportion of Closed Savanna', y = 'Proportion of States') +
  theme_bw()


df.state.open_savanna<- dat2 %>% 
  filter(z.post.thresh != 'Unclassified') %>% 
  mutate(across(c(Forest:Floodable), factor)) %>% 
  group_by(Open_Savanna, z.post.thresh) %>% 
  summarize(n = n()) %>% 
  group_by(Open_Savanna) %>% 
  mutate(prop = n/sum(n))

ggplot(df.state.open_savanna, aes(as.numeric(as.character(Open_Savanna)), prop,
                                  fill = z.post.thresh)) +
  geom_area() +
  scale_fill_viridis_d(option = 'inferno') +
  labs(x = 'Proportion of Open Savanna', y = 'Proportion of States') +
  theme_bw()


df.state.floodable<- dat2 %>% 
  filter(z.post.thresh != 'Unclassified') %>% 
  mutate(across(c(Forest:Floodable), factor)) %>% 
  group_by(Floodable, z.post.thresh) %>% 
  summarize(n = n()) %>% 
  group_by(Floodable) %>% 
  mutate(prop = n/sum(n))

ggplot(df.state.floodable, aes(as.numeric(as.character(Floodable)), prop,
                               fill = z.post.thresh)) +
  geom_area() +
  scale_fill_viridis_d(option = 'inferno') +
  labs(x = 'Proportion of Floodable', y = 'Proportion of States') +
  theme_bw()



##########################################################################################
### Analyze relationships between states and covars as multinomial logistic regression ###
##########################################################################################


#remove observations w/ "Unclassified" state (based on z.post.thresh)
dat3<- dat2 %>% 
  filter(z.post.thresh != "Unclassified") %>%  #remove Unclassified obs
  mutate(across(id, factor))  #convert ID to factor


## prior predictive check
fit0 <-
  brm(data = dat3, 
      family = categorical(link = logit, refcat = 'VE'),
      z.post.thresh ~ Forest + I(Forest^2) + Open_Savanna + I(Open_Savanna^2) + 
        Floodable + I(Floodable^2) + (1|id),
      prior = c(prior(normal(1, 2), class = Intercept, dpar = muLocalSearch),
                prior(normal(1, 2), class = Intercept, dpar = muExploratory),
                prior(normal(1, 2), class = Intercept, dpar = muTransit),
                prior(normal(0, 5), class = b, dpar = muLocalSearch),
                prior(normal(0, 5), class = b, dpar = muExploratory),
                prior(normal(0, 5), class = b, dpar = muTransit),
                prior(exponential(1), class = sd, dpar = muLocalSearch),
                prior(exponential(1), class = sd, dpar = muExploratory),
                prior(exponential(1), class = sd, dpar = muTransit)),
      iter = 2000, warmup = 1000, cores = 4, chains = 4,
      seed = 123, sample_prior = "only")

pp_check(fit0, 
         type = "bars", 
         ndraws = 200, 
         size = 1/2, 
         fatten = 2)

## run model
fit1 <-
  brm(data = dat3, 
      family = categorical(link = logit, refcat = 'VE'),
      z.post.thresh ~ Forest + I(Forest^2) + Open_Savanna + I(Open_Savanna^2) + 
        Floodable + I(Floodable^2) + (1|id),
      prior = c(prior(normal(1, 2), class = Intercept, dpar = muLocalSearch),
                prior(normal(1, 2), class = Intercept, dpar = muExploratory),
                prior(normal(1, 2), class = Intercept, dpar = muTransit),
                prior(normal(0, 5), class = b, dpar = muLocalSearch),
                prior(normal(0, 5), class = b, dpar = muExploratory),
                prior(normal(0, 5), class = b, dpar = muTransit),
                prior(exponential(1), class = sd, dpar = muLocalSearch),
                prior(exponential(1), class = sd, dpar = muExploratory),
                prior(exponential(1), class = sd, dpar = muTransit)),
      iter = 2000, warmup = 1000, cores = 4, chains = 4,
      seed = 123)

print(fit1)
plot(fit1)  #traceplots and density plots for each param

ce.95 <- conditional_effects(
  fit1, 
  categorical = TRUE,
  effects = c("Forest", "Open_Savanna", "Floodable"),
  prob = 0.95,  #95% quantile
  method = "fitted",
  resolution = 1000)
ce.50 <- conditional_effects(
  fit1, 
  categorical = TRUE,
  effects = c("Forest", "Open_Savanna", "Floodable"),
  prob = 0.5,  #50% quantile
  method = "fitted",
  resolution = 1000)


ce.95.df <- ce.95 %>% 
  map({. %>% dplyr::select(effect1__:upper__)}) %>% 
  bind_rows(.id = "lulc") %>% 
  mutate(lulc = gsub(":.*", "", lulc)) %>% 
  mutate(lulc = gsub("Open_Savanna", "Open Savanna", lulc)) %>% 
  mutate(across(lulc, factor, levels = c('Forest', 'Open Savanna', 'Floodable'))) %>% 
  mutate(prob = 0.95)
ce.50.df <- ce.50 %>% 
  map({. %>% dplyr::select(effect1__:upper__)}) %>% 
  bind_rows(.id = "lulc") %>% 
  mutate(lulc = gsub(":.*", "", lulc)) %>% 
  mutate(lulc = gsub("Open_Savanna", "Open Savanna", lulc)) %>% 
  mutate(across(lulc, factor, levels = c('Forest', 'Open Savanna', 'Floodable'))) %>% 
  mutate(prob = 0.5)

ggplot(ce.95.df, aes(effect1__, estimate__, group = effect2__)) +
  geom_line(aes(color = effect2__), size = 2) +
  scale_color_viridis_d("", option = "inferno", begin = 0, end = 0.8) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = effect2__), alpha = 0.6) +
  scale_fill_viridis_d("", option = "inferno") +
  geom_ribbon(data = ce.50.df, aes(ymin = lower__, ymax = upper__, fill = effect2__),
              alpha = 0.5) +
  theme_bw() +
  ylab("Probability") +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1)) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 20),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.text = element_text(size = 16),
        legend.position = "top",
        panel.spacing.x = unit(1, "lines")) +
  facet_wrap(~ lulc, nrow = 1, strip.position = "bottom")

# ggsave('Figure 5_new.png', width = 12, height = 9, units = "in", dpi = 330)


## posterior predictive check
pp_check(fit1, 
         type = "bars", 
         ndraws = 200, 
         size = 1/2, 
         fatten = 2)



## Create density plot of parameter estimates from posterior
fit1.post <- fit1 %>% 
  as_draws_df(variable = "^b_", regex = TRUE) %>% 
  pivot_longer(cols = 1:21, names_to = "coeff", values_to = "value") %>% 
  mutate(state = case_when(str_detect(string = coeff, pattern = "LocalSearch") ~ "Local Search",
                           str_detect(string = coeff, pattern = "Exploratory") ~ "Exploratory",
                           str_detect(string = coeff, pattern = "Transit") ~ "Transit"),
         coeff.name = case_when(str_detect(string = coeff, pattern = "Intercept") ~ "Intercept",
                                str_detect(string = coeff, pattern = "Forest$") ~ "Forest",
                                str_detect(string = coeff, pattern = "ForestE2") ~ "Forest^2",
                                str_detect(string = coeff, pattern = "Open_Savanna$") ~ "Open Savanna",
                                str_detect(string = coeff, pattern = "Open_SavannaE2") ~ "Open Savanna^2",
                                str_detect(string = coeff, pattern = "Floodable$") ~ "Floodable",
                                str_detect(string = coeff, pattern = "FloodableE2") ~ "Floodable^2")
  )

fit1.post$state <- factor(fit1.post$state, levels = c("Local Search", "Exploratory", "Transit"))
fit1.post$coeff.name <- factor(fit1.post$coeff.name,
                               levels = c("Intercept", "Forest", "Forest^2", "Open Savanna",
                                          "Open Savanna^2", "Floodable", "Floodable^2"))


ggplot(fit1.post, aes(x = value, y = coeff.name, fill = state, color = state)) +
  geom_vline(xintercept = 0, size = 1) +
  stat_slab(alpha = 0.5) +
  labs(x = "Value", y = "") +
  scale_fill_manual("State",
                    values = viridis(n=4, begin = 0, end = 0.95, option = "inferno")[2:4]) +
  scale_color_manual("State",
                     values = viridis(n=4, begin = 0, end = 0.90, option = "inferno")[2:4]) +
  scale_x_continuous(breaks = seq(-8, 8, by = 2), limits = c(-8,8)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "top")

# ggsave('Figure S2.png', width = 8, height = 6, units = "in", dpi = 330)






########################
### Multinomial GAMM ###
########################

## If needed for a future purpose, multinomial regressions can be run in {mgcv}; still would need to perform additional post-hoc data summarization to assess whether particular states can be attributed to specific land classes

dat3$z.post.thresh<- fct_drop(dat3$z.post.thresh)
dat3$state<- as.numeric(dat3$z.post.thresh) - 1


table(dat3$id, dat3$state)  #check n per combo
mod1<- gam(list(state ~ s(Forest, k = 3, bs = "cs") + s(Open_Savanna, k = 3, bs = "cs") + 
                  s(Floodable, k = 3, bs = "cs") + s(id, bs = "re"),
                ~ s(Forest, k = 3, bs = "cs") + s(Open_Savanna, k = 3, bs = "cs") + 
                  s(Floodable, k = 3, bs = "cs") + s(id, bs = "re"),
                ~ s(Forest, k = 3, bs = "cs") + s(Open_Savanna, k = 3, bs = "cs") + 
                  s(Floodable, k = 3, bs = "cs") + s(id, bs = "re")),
           data = dat3, family = multinom(K=3), method = "REML")

summary(mod1)
plot(mod1, pages=1, trans = plogis)



## Manually create partial effects plots
pred.multinom.mod<- predict(mod1, newdata = new.dat, type = "response", se.fit = TRUE,
                     exclude = "s(id)")


pred.multinom.df<- pred.multinom.mod$fit %>% 
  data.frame() %>% 
  rename(VE = X1, LocalSearch = X2, Exploratory = X3, Transit = X4) %>% 
  mutate(#lower = plogis(fit - 1.96*se),
         #upper = plogis(fit + 1.96*se),
         x = rep(x.seq, 3),
         LULC = factor(rep(c("Forest", "Open Savanna", "Floodable"), each = 100),
                       levels = c("Forest", "Open Savanna", "Floodable"))
  ) %>% 
  pivot_longer(cols = VE:Transit, names_to = "state", values_to = "prob") #%>% 
  # mutate(fit = plogis(fit),
  #        mod = "Pr(Transit)")


ggplot(pred.multinom.df, aes(x, prob, color = state)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ LULC, ncol = 1)

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
dat2<- dat2 %>% 
  bayesmove::df_to_list("id") %>% 
  map(nocturnal_period) %>% 
  bind_rows()

# Calculate daily distance traveled by ID
daily.dist<- dat2 %>% 
  group_by(id, period) %>% 
  summarize(dist = sum(sl2), month = unique(month))


# Calculate daily max displacement from initial location
daily.displ<- dat2 %>% 
  bayesmove::df_to_list("id") %>% 
  map(bayesmove::df_to_list, "period")

for (i in 1:length(daily.displ)) {
  
  tmp<- daily.displ[[i]]
  n.period<- length(tmp)
  
  for (j in 1:n.period) {
    
    boop<- tmp[[j]]
    boop$displ<- sqrt((boop$x - boop$x[1])^2 + (boop$y - boop$y[1])^2)
    
    tmp[[j]]<- boop
  }
  
  daily.displ[[i]]<- tmp
}
  
daily.displ2<- flatten(daily.displ) %>% 
  bind_rows() %>% 
  mutate_at("period", as.numeric)

daily.displ3<- daily.displ2 %>% 
  group_by(id, period) %>% 
  summarize(displ = max(displ, na.rm = T), month = unique(month))



## Plot monthly trends by ID

#distance
ggplot(daily.dist, aes(month, dist, fill = id)) +
  # geom_point(size = 2, alpha = 0.5) +
  geom_boxplot() +
  labs(y="Total distance traveled per day (m)", x = "Month") +
  theme_bw() +
  facet_wrap(~id)

#max displacement
ggplot(daily.displ3[-332,], aes(month, displ, fill = id)) +
  # geom_point(size = 2, alpha = 0.5) +
  geom_boxplot() +
  labs(y="Max displacement per day (m)", x = "Month") +
  theme_bw() +
  facet_wrap(~id)


## Summarize daily max displacement (mean and range)
daily.displ3<- daily.displ3 %>% 
  filter(displ < 10000)  #remove weird high value from mazeboti
daily.displ3 %>% 
  group_by(id) %>% 
  summarize(mean = mean(displ), lo = min(displ), hi = max(displ))

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

