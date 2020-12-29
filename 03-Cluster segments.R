#### Segment armadillo acceleration ####

library(bayesmove)
library(tidyverse)
library(lubridate)

### Import Data ###

dat<- read.csv("Segmented Armadillo Acceleration.csv", as.is = T)
dat$date<- as_datetime(dat$date)


# Select only id, tseg, AC, SL, and TA
dat2<- dat[,c("id","tseg","AC","SL","TA")]

# Summarize observations by track segment
nbins<- c(6,6,8)
obs<- summarize_tsegs(dat = dat2, nbins = nbins)
head(obs)




### Run LDA model ###

set.seed(1)

# Prepare for Gibbs sampler
ngibbs<- 1000  #number of MCMC iterations for Gibbs sampler
nburn<- ngibbs/2  #number of iterations for burn-in
nmaxclust<- max(nbins) - 1  #one fewer than max number of bins used for data streams
ndata.types<- length(nbins)  #number of data types

# Priors
gamma1<- 0.1
alpha<- 0.1

# Run LDA model
res<- cluster_segments(dat=obs, gamma1=gamma1, alpha=alpha,
                       ngibbs=ngibbs, nmaxclust=nmaxclust,
                       nburn=nburn, ndata.types=ndata.types)
# takes 1 min to run on 1000 iterations
# takes 12 s to run using only SL and TA

# Check traceplot of log likelihood
plot(res$loglikel, type='l', xlab = "Iteration", ylab = "Log Likelihood")




### Visualize results from theta and phi ###

# Extract proportions of behaviors per track segment
theta.estim<- extract_prop(res = res, ngibbs = ngibbs, nburn = nburn, nmaxclust = nmaxclust)

# Convert to data frame for ggplot2
theta.estim_df<- theta.estim %>%
  as.data.frame() %>%
  pivot_longer(., cols = 1:all_of(nmaxclust), names_to = "behavior", values_to = "prop") %>%
  modify_at("behavior", factor)
levels(theta.estim_df$behavior)<- 1:nmaxclust

# Plot results
ggplot(theta.estim_df, aes(behavior, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1),
              alpha = 0.3) +
  labs(x="\nBehavior", y="Proportion of Total Behavior\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))


# Calculate mean proportions per behavior
(theta.means<- round(colMeans(theta.estim), digits = 3))

# Calculate cumulative sum
cumsum(theta.means)  #first 2 clusters account for 94% of obs; stick with these




# Extract bin estimates from phi matrix
behav.res<- get_behav_hist(dat = res, nburn = nburn, ngibbs = ngibbs, nmaxclust = nmaxclust,
                           var.names = c("Activity Count","Step Length", "Turning Angle"))

# Plot histograms of proportion data
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(2), "grey35",
                               "grey35","grey35","grey35","grey35"), guide = FALSE) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")




# Reformat proportion estimates for all track segments
theta.estim.long<- expand_behavior(dat = dat, theta.estim = theta.estim, obs = obs,
                                   nbehav = 2, behav.names = c("High", "Low"),
                                   behav.order = 2:1)


tol<- 60*24*7  #1 week in mins
units<- "min"
# Add gaps when dt > 1 week (in minutes)
theta.estim.long2<- insert_NAs(data = theta.estim.long, int = 7, units = "mins")



# Plot results
ggplot(theta.estim.long) +
  geom_area(aes(x=date, y=prop, fill = behavior), color = "black", size = 0.25,
              position = "fill") +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id, nrow = 7, scales = "free_x")




### Merge results, visualize, and export ###

# Load 'original' cleaned data
dat.orig<- read.csv("Binned Armadillo Acceleration Data.csv", as.is = T)
# dat.orig$obs<- dat.orig$time1 

# Convert segmented dataset into list
# dat$obs<- dat$time1
dat.list<- df_to_list(dat = dat, ind = "id")

# Merge results with original data
dat.out<- assign_behavior(dat.orig = dat.orig,
                             dat.seg.list = dat.list,
                             theta.estim.long = theta.estim.long,
                             behav.names = c("Low","High"))
# dat.out$behav<- factor(dat.out$behav, levels = c("Low", "High"))




ggplot() +
  geom_path(data = dat.out, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = dat.out, aes(x, y, fill=behav), size=2.5, pch=21, alpha=0.7) +
  geom_point(data = dat.out %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3,
             stroke = 1.25) +
  geom_point(data = dat.out %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  scale_fill_viridis_d("") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top") +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(~id, scales = "free")



## Explore in Shiny

dat.out$date<- as_datetime(dat.out$date)
dat.out2<- filter(dat.out, !is.na(x))
dat.out2$behav<- as.numeric(dat.out2$behav)

shiny_tracks(dat.out2, 32721)
