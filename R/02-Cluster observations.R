
#### Estimate behavioral states from armadillo data using bayesmove mixture model ####
## Last revised 2022-09-29 by Josh Cullen

library(bayesmove)
library(tidyverse)
library(lubridate)


#############################
### Import armadillo data ###
#############################

dat<- read.csv('Data_processed/Giant Armadillo Binned Data.csv', as.is = T)
dat$date<- as_datetime(dat$date, tz = "Etc/GMT+4")  #change timezone to subtract 4 hours


dat2<- subset(dat, select = c(sl.cat, ta.cat, act.cat))

## Map available points for each ID
ggplot(dat, aes(easting, northing, color = id)) +
  geom_path(aes(group = id)) +
  geom_point(size = 1, alpha = 0.75) +
  theme_bw() +
  coord_fixed()






#################################################
### Run non-parametric Bayesian mixture model ###
#################################################

set.seed(3)

# Define model params
alpha = 0.1  #prior
ngibbs = 20000  #number of Gibbs sampler iterations
nburn = ngibbs/2  #number of burn-in iterations
nmaxclust = 10  #number of maximum possible states (clusters) present

# Run model
dat.res<- cluster_obs(dat=dat2, alpha=alpha, ngibbs=ngibbs, nmaxclust=nmaxclust,
                      nburn=nburn)
# takes 18 min to run 20,000 iterations

# Inspect traceplot of log-likelihood
post.seq = (nburn + 1):ngibbs
plot(dat.res$loglikel, type = "l")
plot(dat.res$loglikel[post.seq], type = "l")

plot(dat.res$gamma1, type = "l")
plot(dat.res$gamma1[post.seq], type = "l")

par(mfrow=c(2,2), ask = T)
for (i in 1:nmaxclust) {
  plot(dat.res$theta[,i], type = "l", main = paste("State", i))
}
par(mfrow = c(1,1), ask=F)


## Inspect and Plot results

theta<- dat.res$theta[post.seq,]
colnames(theta)<- 1:ncol(theta)
theta1<- colMeans(theta)
theta1<- sort(theta1, decreasing = TRUE)
cumsum(theta1)  #3 or 4 states seem possible

max.gr=4
mean(rowSums(theta[,1:max.gr])) #98%
median(rowSums(theta[,1:max.gr])) #98%
range(rowSums(theta[,1:max.gr])) #93%-100%


### Check that Phi's match up with same state across iterations (i.e., no state flipping)
par(mfrow = c(4,3))
for (i in sample(post.seq, 4)) {
  plot(matrix(dat.res$phi[[1]][i,], nrow = nmaxclust)[4,], type = "h")
  plot(matrix(dat.res$phi[[2]][i,], nrow = nmaxclust)[4,], type = "h")
  plot(matrix(dat.res$phi[[3]][i,], nrow = nmaxclust)[4,], type = "h")
}
par(mfrow=c(1,1))
## Checks out, no state flipping when plotting random samples from posterior



## Viz state-dependent distributions
# Extract bin estimates for each possible state from the `phi` matrix of the model results
behav.res<- get_behav_hist(dat = dat.res, nburn = nburn, ngibbs = ngibbs, nmaxclust = nmaxclust,
                           var.names = c("Speed","Turning Angle","Activity Counts"))

ggplot(behav.res %>% filter(behav %in% 1:4), aes(x = bin, y = prop, fill = behav)) +
  geom_col(color = "black") +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  facet_grid(behav ~ var, scales = "free_x")

behav.res2<- behav.res %>% 
  filter(behav %in% 1:max.gr) %>% 
  mutate_at("behav", factor, levels = c(3,1,2,4))
levels(behav.res2$behav)<- c("VE", "Local Search", "Exploratory", "Transit")


# Rename bin labels on continuous scale
act.lims<- data.frame(brk = seq(from = 0, to = 300, by = 50),
                      var = 'Activity Counts')
sl.lims<- data.frame(brk = round(c(seq(from = 0, to = 60, by = 10), max(dat$sl2)) / 60, 2),
                     var = 'Speed')  #convert to m/s
ta.lims<- data.frame(brk = round(seq(from = -pi, to = pi, length.out = 11), 2),
                    var = 'Turning Angle')
lims<- rbind(act.lims, sl.lims, ta.lims)

new.lab<- data.frame()
for (i in 1:length(unique(behav.res2$var))) {
  ind<- unique(behav.res2$var)[i]
  tmp<- lims %>% 
    filter(var == ind)
  
  for (j in 2:nrow(tmp)) {
    new.lab<- rbind(new.lab, 
                    c(paste0(tmp[j-1, 'brk'], ', ', tmp[j, 'brk']),
                      unique(tmp$var),
                      j-1))
  }
}
names(new.lab)<- c('labs', 'var', 'bin')
new.lab$bin<- as.integer(new.lab$bin)
behav.res3<- left_join(behav.res2, new.lab, by = c('var', 'bin'))
behav.res3$labs<- factor(behav.res3$labs, levels = unique(behav.res3$labs))

# Plot state-dependent distributions 
ggplot(behav.res3, aes(x = labs, y = prop, fill = behav)) +
  geom_col(color = "black") +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold"),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_viridis_d(option = 'inferno', guide = "none") +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav ~ var, scales = "free_x")

# ggsave("Figures/Figure 4.png", width = 7, height = 7, units = "in", dpi = 400)



## Attribute behaviors to states and extract each of the different estimates
# Using MAP estimate, threshold of 75% assignments from posterior, and most common state
z.post<- as.matrix(dat.res$z.posterior)
z.post2<- t(apply(z.post, 1, function(x) x/sum(x)))
thresh<- 0.75
z.post3<- apply(z.post2, 1, function(x) ifelse(max(x) > thresh, which(x > thresh), NA))
z.post4<- apply(z.post2, 1, function(x) which.max(x))

# Compare posterior estimates against MAP
table(z.post3, dat.res$z.MAP)  #estimates are close, but some discrepancies
table(z.post4, dat.res$z.MAP)  #estimates are close, but some discrepancies


## Map estimated states
dat<- dat %>% 
  mutate(z.map = dat.res$z.MAP,
         z.post.thresh = z.post3,
         z.post.max = z.post4)
dat$z.map<- ifelse(dat$z.map > max.gr, NA, dat$z.map)
dat$z.post.max<- ifelse(dat$z.post.max > max.gr, NA, dat$z.post.max)

# Attribute behaviors to states
dat2<- dat %>% 
  mutate(across(c('z.map','z.post.thresh','z.post.max'),
                ~case_when(. == 1 ~ "Local Search",
                           . == 2 ~ "Exploratory",
                           . == 3 ~ "VE",
                           . == 4 ~ "Transit",
                           is.na(.) ~ "Unclassified")
                )) %>% 
  mutate(across(c('z.map','z.post.thresh','z.post.max'),
                factor, levels = c('VE','Local Search','Exploratory','Transit',
                                   'Unclassified')
                )) %>% 
  mutate_at("id", str_to_title)


#by ID
ggplot(dat2, aes(easting, northing, color = id)) +
  geom_path(aes(group = id), color = "grey80", size = 0.5) +
  geom_point(size = 1, alpha = 0.3) +
  theme_bw() +
  facet_grid(~id, scales = "free", space = "free")

#by z.map
ggplot(dat2, aes(easting, northing, color = z.map)) +
  scale_color_manual("", values = c(viridis::viridis(4), "grey")) +
  geom_path(aes(group = id), color = "grey80", size = 0.5) +
  geom_point(size = 1) +
  theme_bw()

#by z.post.thresh
ggplot(dat2, aes(easting, northing, fill = z.post.thresh)) +
  scale_fill_manual("State", values = c(viridis::viridis(4, option = "inferno"), "grey")) +
  geom_path(aes(group = id), color = "grey80", size = 0.5) +
  geom_point(size = 1, shape = 21) +
  theme_bw() +
  labs(x = "Easting", y = "Northing") +
  facet_wrap(~id, scales = "free") +
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5)))



ggplot(dat2 %>% filter(id == "Mafalda"), aes(easting, northing, fill = z.post.thresh)) +
  scale_fill_manual("State", values = c(viridis::viridis(4, option = "inferno"), "grey")) +
  geom_path(aes(group = id), color = "grey80", size = 0.5) +
  geom_point(size = 2, shape = 21) +
  theme_bw() +
  labs(x = "Easting", y = "Northing") +
  facet_wrap(~id, scales = "free") +
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5)))



#by z.post.max
ggplot(dat2, aes(easting, northing, color = z.post.max)) +
  scale_color_manual("", values = c(viridis::viridis(4), "grey")) +
  geom_path(aes(group = id), color = "grey80", size = 0.5) +
  geom_point(size = 1) +
  theme_bw()


dat2 %>%   # for estimates based on threshold
  group_by(z.post.thresh) %>% 
  tally() %>% 
  mutate(prop = n/sum(n))

dat2 %>%   # for MAP estimates
  group_by(z.map) %>% 
  tally() %>% 
  mutate(prop = n/sum(n))

### Over 1/3 of all observations have an "Unclassified" state based on posterior estimates
### Only 0.7% of all observations are "Unclassified" based on MAP estimate




### Export results
behav.res<- behav.res %>% 
  filter(behav %in% 1:max.gr) %>% 
  mutate(across(c('behav'),
                ~case_when(. == 1 ~ "Local Search",
                           . == 2 ~ "Exploratory",
                           . == 3 ~ "VE",
                           . == 4 ~ "Transit",
                           is.na(.) ~ "Unclassified")
  ))

# write.csv(dat2, "Data_processed/Giant Armadillo state estimates.csv", row.names=F)
# write.csv(behav.res, "Data_processed/Giant Armadillo behavior distribs.csv", row.names=F)







######################################################################
### Run non-parametric Bayesian mixture model using only SL and TA ###
######################################################################

## Run model w/ only SL and TA for comparison

dat3<- subset(dat, select = c(sl.cat, ta.cat))  #check states for only SL and TA

set.seed(3)

# Define model params
alpha = 0.1  #prior
ngibbs = 20000  #number of Gibbs sampler iterations
nburn = ngibbs/2  #number of burn-in iterations
nmaxclust = 10  #number of maximum possible states (clusters) present

# Run model
dat.res2<- cluster_obs(dat=dat3, alpha=alpha, ngibbs=ngibbs, nmaxclust=nmaxclust,
                      nburn=nburn)
# takes 15 min to run 20,000 iterations

# Inspect traceplot of log-likelihood
post.seq <- (nburn + 1):ngibbs
plot(dat.res2$loglikel, type = "l")
plot(dat.res2$loglikel[post.seq], type = "l")


## Inspect and Plot results
theta<- dat.res2$theta[post.seq,]
colnames(theta)<- 1:ncol(theta)
theta1<- colMeans(theta)
theta1<- sort(theta1, decreasing = TRUE)
cumsum(theta1)  #3 states optimal

max.gr=3

## Viz state-dependent distributions
# Extract bin estimates for each possible state from the `phi` matrix of the model results
behav.res<- get_behav_hist(dat = dat.res2, nburn = nburn, ngibbs = ngibbs, nmaxclust = nmaxclust,
                           var.names = c("Speed","Turning Angle"))

ggplot(behav.res %>% filter(behav %in% 1:3), aes(x = bin, y = prop, fill = behav)) +
  geom_col(color = "black") +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  facet_grid(behav ~ var, scales = "free_x")

behav.res2<- behav.res %>% 
  filter(behav %in% 1:max.gr) %>% 
  mutate_at("behav", factor, levels = c(2,1,3))
levels(behav.res2$behav)<- c("VE", "Local Search", "Exploratory")

# Rename bin labels on continuous scale
sl.lims<- data.frame(brk = round(c(seq(from = 0, to = 60, by = 10), max(dat$sl2)) / 60, 2),
                     var = 'Speed')  #convert to m/s
ta.lims<- data.frame(brk = round(seq(from = -pi, to = pi, length.out = 11), 2),
                     var = 'Turning Angle')
lims<- rbind(sl.lims, ta.lims)

new.lab<- data.frame()
for (i in 1:length(unique(behav.res2$var))) {
  ind<- unique(behav.res2$var)[i]
  tmp<- lims %>% 
    filter(var == ind)
  
  for (j in 2:nrow(tmp)) {
    new.lab<- rbind(new.lab, 
                    c(paste0(tmp[j-1, 'brk'], ', ', tmp[j, 'brk']),
                      unique(tmp$var),
                      j-1))
  }
}
names(new.lab)<- c('labs', 'var', 'bin')
new.lab$bin<- as.integer(new.lab$bin)
behav.res3<- left_join(behav.res2, new.lab, by = c('var', 'bin'))
behav.res3$labs<- factor(behav.res3$labs, levels = unique(behav.res3$labs))

# Plot state-dependent distributions 
ggplot(behav.res3, aes(x = labs, y = prop, fill = behav)) +
  geom_col(color = "black") +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold"),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_viridis_d(option = 'inferno', guide = "none") +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav ~ var, scales = "free_x")

# ggsave("Figures/Figure S1.png", width = 7, height = 7, units = "in", dpi = 330)
