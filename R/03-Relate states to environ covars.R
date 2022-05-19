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

source('R/helper functions.R')


#################
### Load data ###
#################

dat<- read.csv("Data_processed/Giant Armadillo state estimates.csv", as.is = T)
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

# ggsave("Figures/Figure 1.png", width = 9, height = 5, units = "in", dpi = 400)





### Load LU/LC layers for flood and dry seasons ###

#LULC
flood.lulc<- raster('Environ_vars/cheiann_UTM1.tif')
names(flood.lulc)<- 'lulc'

#crop raster closer to armadillo relocs
flood.lulc<- crop(flood.lulc, extent(dat %>% 
                           summarize(xmin = min(x) - 3000,
                                     xmax = max(x) + 3000,
                                     ymin = min(y) - 3000,
                                     ymax = max(y) + 3000) %>% 
                           unlist()))
flood.lulc.df<- as.data.frame(flood.lulc, xy = TRUE)


dry.lulc<- raster('Environ_vars/secann_UTM1.tif')
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

# ggsave("Figures/Figure 2.tiff", width = 7, height = 4, units = "in", dpi = 400)



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
# ggsave("Figures/Figure 4.png", width = 12, height = 8, dpi = 400)








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
      seed = 123, backend = "cmdstanr")
#took 11.5 min to run

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

# ggsave('Figures/Figure 5.png', width = 12, height = 9, units = "in", dpi = 400)


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

# ggsave('Figures/Figure S1.png', width = 8, height = 6, units = "in", dpi = 330)



## Create plots of expected (mean) behavioral state probability when in 100% of each LU/LC

newdata <- data.frame(Forest = c(1, 0, 0, 0),
                       `I(Forest^2)` = c(1, 0, 0, 0),
                       Open_Savanna = c(0, 1, 0, 0),
                       `I(Open_Savanna^2)` = c(0, 1, 0, 0),
                       Floodable = c(0, 0, 1, 0),
                       `I(Floodable^2)` = c(0, 0, 1, 0),
                      id = NA
                      )
names(newdata) <- c('Forest', "I(Forest^2)", 'Open_Savanna', "I(Open_Savanna^2)",
                    'Floodable', "I(Floodable^2)", "id")

pred.lulc.1 <- fit1 %>% 
  posterior_epred(newdata, re_formula = NA)


# For plotting full distribution of state proportions
dims <- dim(pred.lulc.1)
pred.distr <- expand.grid(state = c('VE', 'Local Search', 'Exploratory', 'Transit'),
                          class = c('Forest', 'Open Savanna', 'Floodable', 'Closed Savanna'),
                          row1 = 1:4000) %>% 
  arrange(state, class) %>% 
  mutate(prob = array(apply(pred.lulc.1, 3, rbind), dims[1] * dims[2] * dims[3]))


ggplot(pred.distr, aes(class, prob, fill = state)) +
  geom_violin(position = position_dodge(0.3)) +
  scale_fill_viridis_d("", option = "inferno") +
  labs(x = "", y = "Probability") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "top")

# ggsave('Figures/Figure S2.png', width = 8, height = 6, units = "in", dpi = 330)


