### Relate state estimates to environmental covariates

library(tidyverse)
library(lubridate)
library(raster)
library(viridis)
library(mgcv)
library(mgcViz)
library(wesanderson)
library(rnaturalearth)
library(rnaturalearthhires)
library(cowplot)
library(sf)
library(ggspatial)

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
dat$season<- ifelse(dat$month %in% c(month.abb[3:5]), "Fall",
                    ifelse(dat$month %in% c(month.abb[6:8]), "Winter",
                           ifelse(dat$month %in% c(month.abb[9:11]), "Spring", "Summer")))
dat$season<- factor(dat$season, levels = c("Fall","Winter","Spring","Summer"))
dat$season2<- ifelse(dat$month %in% month.abb[c(10:12,1:3)], "Wet", "Dry")


#Plot of all tracks together
pal2<- c(wes_palettes$Darjeeling1, wes_palettes$Darjeeling2[2:1])  #palette
dat$id<- str_to_title(dat$id)

dat.sf<- dat %>% 
  st_as_sf(., coords = c("x", "y"), crs = 32721) %>% 
  group_by(id) %>% 
  dplyr::summarize(do_union = F) %>%
  sf::st_cast("LINESTRING")


#ESRI REST API for world satellite imagery
esri_world <- paste0('https://services.arcgisonline.com/arcgis/rest/services/',
                     'World_Imagery/MapServer/tile/${z}/${y}/${x}.jpeg')

main.map<- ggplot() +
  annotation_map_tile(type = esri_world, zoom = 13, progress = "none") +
  layer_spatial(dat.sf, size = 0.5, aes(color = id)) +
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
SA<- ne_countries(scale = 'medium', continent = 'south america')
SA.tr<- spTransform(SA, CRS("+init=epsg:32721"))
# SA.df<- fortify(SA.tr)
SA.sf<- st_as_sf(SA.tr)

inset1<- ggplot() +
  geom_sf(data = SA.sf, fill = "grey75", color = "white", size = 0.25) +
  geom_rect(aes(xmin = min(dat$x) - 3e5,
                xmax = max(dat$x) + 3e5,
                ymin = min(dat$y) - 3e5,
                ymax = max(dat$y) + 3e5),
            color = viridis(3, option = "mako")[1], fill = NA, size = 0.5) +
  scale_x_continuous(expand = c(0,0), limits = c(-2.5e6, st_bbox(SA.sf)$xmax)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_void()

inset2<- ggplot() +
  geom_sf(data = SA.sf, fill = "grey75", color = "white", size = 0.25) +
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
  group_by(season2, z.post.thresh) %>% 
  summarize(n = n()) %>% 
  mutate(freq = n / sum(n))

ggplot(df.season, aes(z.post.thresh, freq, fill = z.post.thresh)) +
  geom_col(color = "black") +
  scale_fill_manual("", values = c(viridis(n=4, option = 'inferno'), "grey"),
                    guide = guide_legend(reverse = TRUE)) +
  theme_bw() +
  labs(x = "", y = "Proportion of Observations") +
  facet_wrap(~ season2) +
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
  scale_fill_manual("", values = c(viridis(n=4, option = 'inferno'), "grey"), guide = F) +
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



#######################################################
### Analyze relationships between states and covars ###
#######################################################


#remove observations w/ "Unclassified" state (based on z.post.thresh)
dat3<- dat2 %>% 
  filter(z.post.thresh != "Unclassified") %>% 
  mutate(across(.cols = c(Forest:Floodable), scale)) %>%
  mutate(across(c(id,season2), factor))

#set new scales and labels for plots
att.forest <- attributes(dat3$Forest)
att.cs <- attributes(dat3$Closed_Savanna)
att.os <- attributes(dat3$Open_Savanna)
att.flood <- attributes(dat3$Floodable)

gam.labels <- seq(0, 1, 0.2)
forest.breaks <- scale(gam.labels, att.forest[2], att.forest[3])[,1]
cs.breaks <- scale(gam.labels, att.cs[2], att.cs[3])[,1]
os.breaks <- scale(gam.labels, att.os[2], att.os[3])[,1]
flood.breaks <- scale(gam.labels, att.flood[2], att.flood[3])[,1]

state.pal<- viridis(n=4, option = 'inferno')



## Compare VE vs all others
dat.ve<- dat3
dat.ve$state<- ifelse(dat.ve$z.post.thresh == "VE", 1, 0)
table(dat.ve$id, dat.ve$state)  #check n per combo
ve.mod<- gam(state ~ s(Forest, k = 5, bs = "cs") + s(Closed_Savanna, k = 5, bs = "cs") +
                  s(Open_Savanna, k = 5, bs = "cs") + s(Floodable, k = 5, bs = "cs") +
                  s(id, bs = "re"),
                data = dat.ve, family = binomial, method = "REML")
summary(ve.mod)

# prep data for plotting w/ mgcViz
ve.viz<- getViz(ve.mod)

## define theme for GAMM ggplots
theme_gam <- function(){ 
  font <- "Arial"   #assign font family up front
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 18),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 14),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10))
      
      
    )
}

p.ve1<- plot( sm(ve.viz, 1) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.5, fill = state.pal[1]) + 
  l_fitLine() + 
  labs(y = expression(paste(italic(Pr), ' VE')), x = '') +
  ylim(0.35, 0.72) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = forest.breaks) +
  theme_gam()
p.ve2<- plot( sm(ve.viz, 2) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.5, fill = state.pal[1]) + 
  l_fitLine() + 
  labs(y = '', x = '') +
  ylim(0.35, 0.72) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = cs.breaks) +
  theme_gam()
p.ve3<- plot( sm(ve.viz, 3) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.5, fill = state.pal[1]) + 
  l_fitLine() + 
  labs(y = '', x = '') +
  ylim(0.35, 0.72) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = os.breaks) +
  theme_gam()
p.ve4<- plot( sm(ve.viz, 4) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.5, fill = state.pal[1]) + 
  l_fitLine() + 
  labs(y = '', x = '') +
  ylim(0.35, 0.72) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = flood.breaks) +
  theme_gam()

#define covariates in vector for axis labels
# lulc<- c("Forest", "Closed Savanna", "Open Savanna", "Floodable")
# par(mfrow = c(1,4))
# for (i in 1:4) {
#   y.title<- ifelse(i == 1, "Probability of Local Search", "")
#   plot(st.su.mod, seWithMean = T, rug = T, shade = T, select = i,
#        xlab = lulc[i], ylab = y.title, cex.axis = 1, cex.lab = 1)
#   abline(h = 0.5, lty = 2)
# }





## Compare Local Search vs all others
dat.ls<- dat3
dat.ls$state<- ifelse(dat.ls$z.post.thresh == "Local Search", 1, 0)
table(dat.ls$id, dat.ls$state)  #check n per combo
ls.mod<- gam(state ~ s(Forest, k = 5, bs = "cs") + s(Closed_Savanna, k = 5, bs = "cs") +
                  s(Open_Savanna, k = 5, bs = "cs") + s(Floodable, k = 5, bs = "cs") +
                  s(id, bs = "re"),
                data = dat.ls, family = binomial, method = "REML")
summary(ls.mod)

# prep data for plotting w/ mgcViz
ls.viz<- getViz(ls.mod)

p.ls1<- plot( sm(ls.viz, 1) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.5, fill = state.pal[2]) + 
  l_fitLine() + 
  labs(y = expression(paste(italic(Pr), ' Local Search')), x = '') +
  ylim(0.4, 0.6) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = forest.breaks) +
  theme_gam()
p.ls2<- plot( sm(ls.viz, 2) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.5, fill = state.pal[2]) + 
  l_fitLine() + 
  labs(y = '', x = '') +
  ylim(0.4, 0.6) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = cs.breaks) +
  theme_gam()
p.ls3<- plot( sm(ls.viz, 3) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.5, fill = state.pal[2]) + 
  l_fitLine() + 
  labs(y = '', x = '') +
  ylim(0.4, 0.6) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = os.breaks) +
  theme_gam()
p.ls4<- plot( sm(ls.viz, 4) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.5, fill = state.pal[2]) + 
  l_fitLine() + 
  labs(y = '', x = '') +
  ylim(0.44, 0.56) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = flood.breaks) +
  theme_gam()

# for (i in 1:4) {
#   plot(st.exp.mod, trans = plogis, seWithMean = T, rug = T, shade = T, select = i,
#        xlab = lulc[i], ylab = "Probability of Exploratory", cex.axis = 1, cex.lab = 1)
#   abline(h = 0.5, lty = 2)
# }




## Compare Exploratory vs all others
dat.exp<- dat3
dat.exp$state<- ifelse(dat.exp$z.post.thresh == "Exploratory", 1, 0)
table(dat.exp$id, dat.exp$state)  #check n per combo
exp.mod<- gam(state ~ s(Forest, k = 5, bs = "cs") + s(Closed_Savanna, k = 5, bs = "cs") +
                   s(Open_Savanna, k = 5, bs = "cs") + s(Floodable, k = 5, bs = "cs") +
                   s(id, bs = "re"),
                 data = dat.exp, family = binomial, method = "REML")
summary(exp.mod)

# prep data for plotting w/ mgcViz
exp.viz<- getViz(exp.mod)

p.exp1<- plot( sm(exp.viz, 1) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.5, fill = state.pal[3]) + 
  l_fitLine() + 
  labs(y = expression(paste(italic(Pr), ' Exploratory')), x = '') +
  ylim(0.3, 0.7) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = forest.breaks) +
  theme_gam()
p.exp2<- plot( sm(exp.viz, 2) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.5, fill = state.pal[3]) + 
  l_fitLine() + 
  labs(y = '', x = '') +
  ylim(0.3, 0.7) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = cs.breaks) +
  theme_gam()
p.exp3<- plot( sm(exp.viz, 3) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.5, fill = state.pal[3]) + 
  l_fitLine() + 
  labs(y = '', x = '') +
  ylim(0.3, 0.7) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = os.breaks) +
  theme_gam()
p.exp4<- plot( sm(exp.viz, 4) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.5, fill = state.pal[3]) + 
  l_fitLine() + 
  labs(y = '', x = '') +
  ylim(0.3, 0.7) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = flood.breaks) +
  theme_gam()


# for (i in 1:4) {
#   plot(st.t.mod, trans = plogis, seWithMean = T, rug = T, shade = T, select = i,
#        xlab = lulc[i], ylab = "Probability of Transit", cex.axis = 1, cex.lab = 1)
#   abline(h = 0.5, lty = 2)
# }



## Compare Transit vs all others
dat.t<- dat3
dat.t$state<- ifelse(dat.t$z.post.thresh == "Transit", 1, 0)
table(dat.t$id, dat.t$state)  #check n per combo
t.mod<- gam(state ~ s(Forest, k = 5, bs = "cs") + s(Closed_Savanna, k = 5, bs = "cs") +
                   s(Open_Savanna, k = 5, bs = "cs") + s(Floodable, k = 5, bs = "cs") +
                   s(id, bs = "re"),
                 data = dat.t, family = binomial, method = "REML")
summary(t.mod)

# prep data for plotting w/ mgcViz
t.viz<- getViz(t.mod)

p.t1<- plot( sm(t.viz, 1) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.75, fill = state.pal[4]) + 
  l_fitLine() + 
  labs(y = expression(paste(italic(Pr), ' Transit')), x = 'Forest') +
  ylim(0.4, 0.7) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = forest.breaks) +
  theme_gam()
p.t2<- plot( sm(t.viz, 2) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.75, fill = state.pal[4]) + 
  l_fitLine() + 
  labs(y = '', x = 'Closed Savanna') +
  ylim(0.4, 0.7) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = cs.breaks) +
  theme_gam()
p.t3<- plot( sm(t.viz, 3) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.75, fill = state.pal[4]) + 
  l_fitLine() + 
  labs(y = '', x = 'Open Savanna') +
  ylim(0.4, 0.7) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = os.breaks) +
  theme_gam()
p.t4<- plot( sm(t.viz, 4) , trans = plogis) + 
  # l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 0.75, fill = state.pal[4]) + 
  l_fitLine() + 
  labs(y = '', x = 'Floodable') +
  ylim(0.4, 0.7) +
  # l_points(shape = 19, size = 1, alpha = 0.1) + 
  scale_x_continuous(labels = gam.labels, breaks = flood.breaks) +
  theme_gam()

# for (i in 1:4) {
#   plot(su.exp.mod, trans = plogis, seWithMean = T, rug = T, shade = T, select = i,
#        xlab = lulc[i], ylab = "Probability of Exploratory", cex.axis = 1, cex.lab = 1)
#   abline(h = 0.5, lty = 2)
# }




## Create panel plot
# png('Figure 5.png', width = 12, height = 8, units = "in", res = 330)
gridPrint(p.ve1, p.ve2, p.ve3, p.ve4,
          p.ls1, p.ls2, p.ls3, p.ls4,
          p.exp1, p.exp2, p.exp3, p.exp4,
          p.t1, p.t2, p.t3, p.t4,
          nrow = 4, ncol = 4)
# dev.off()










### Behavioral valuation of landscape
dry.lulc.df2<- dry.lulc.df
dry.lulc.df2$lulc<- ifelse(dry.lulc.df2$lulc == 3, 2, dry.lulc.df2$lulc)
  
pal1<- viridis(4, option = "inferno")

ggplot() +
  geom_raster(data = dry.lulc.df2, aes(x, y, fill = factor(lulc))) +
  scale_fill_manual("", values = c(pal1[1], "grey85", pal1[4]),
                    na.value = "transparent",
                    labels = c("Slow-Turn", "Mixed", "Transit","")) +
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



########################
### Multinomial GAMM ###
########################

## If needed for a future purpose, multinomial regressions can be run in {mgcv}; still would need to perform additional post-hoc data summarization to assess whether particular states can be attributed to specific land classes

dat3$z.post.thresh<- fct_drop(dat3$z.post.thresh)
dat3$state<- as.numeric(dat3$z.post.thresh) - 1


table(dat3$id, dat3$state)  #check n per combo
mod1<- gam(list(state ~ s(Forest, k = 5, bs = "cs") + s(Closed_Savanna, k = 5, bs = "cs") +
                  s(Open_Savanna, k = 5, bs = "cs") + s(Floodable, k = 5, bs = "cs") +
                  s(id, bs = "re"),
                ~ s(Forest, k = 5, bs = "cs") + s(Closed_Savanna, k = 5, bs = "cs") +
                  s(Open_Savanna, k = 5, bs = "cs") + s(Floodable, k = 5, bs = "cs") +
                  s(id, bs = "re"),
                ~ s(Forest, k = 5, bs = "cs") + s(Closed_Savanna, k = 5, bs = "cs") +
                  s(Open_Savanna, k = 5, bs = "cs") + s(Floodable, k = 5, bs = "cs") +
                  s(id, bs = "re")),
           data = dat3, family = multinom(K=3), method = "REML")

summary(mod1)
plot(mod1, pages=1)

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

