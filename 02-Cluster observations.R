
library(bayesmove)
library(tidyverse)
library(lubridate)


#############################
### Import armadillo data ###
#############################

dat<- read.csv('Giant Armadillo Binned Data.csv', as.is = T)
dat$date<- as_datetime(dat$date, tz = "Etc/GMT+4")  #change timezone to subtract 4 hours


dat2<- subset(dat, select = c(sl.cat, ta.cat, act.cat))
# dat3<- subset(dat, select = c(sl.cat, ta.cat))  #check states for only SL and TA

## Map available points for each ID
ggplot(dat, aes(easting, northing, color = id)) +
  geom_path(aes(group = id)) +
  geom_point(size = 1, alpha = 0.75) +
  theme_bw() +
  coord_fixed()



# Run model
set.seed(3)

# Define model params
alpha=0.1  #prior
ngibbs=20000  #number of Gibbs sampler iterations
nburn=ngibbs/2  #number of burn-in iterations
nmaxclust=10  #number of maximum possible states (clusters) present

# Run model
dat.res<- cluster_obs(dat=dat2, alpha=alpha, ngibbs=ngibbs, nmaxclust=nmaxclust,
                      nburn=nburn)
# takes 18 min to run 20,000 iterations

# Inspect traceplot of log-likelihood
post.seq=(nburn + 1):ngibbs
plot(dat.res$loglikel, type = "l")
plot(dat.res$loglikel[post.seq], type = "l")



## Inspect and Plot results

theta<- dat.res$theta[post.seq,]
colnames(theta)<- 1:ncol(theta)
theta1<- colMeans(theta)
theta1<- sort(theta1, decreasing = TRUE)
cumsum(theta1)  #3 states seem optimal; possibly 4

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
## Checks out, no state flipping when plotting random samples from posterior



## Viz state-dependent distributions
# Extract bin estimates for each possible state from the `phi` matrix of the model results
behav.res<- get_behav_hist(dat = dat.res, nburn = nburn, ngibbs = ngibbs, nmaxclust = nmaxclust,
                           var.names = c("Speed","Turning Angle","Activity Counts"))
# behav.res<- get_behav_hist(dat = dat.res, nburn = nburn, ngibbs = ngibbs, nmaxclust = nmaxclust,
#                            var.names = c("Speed","Turning Angle"))

behav.res$behav<- factor(behav.res$behav, levels = 1:nmaxclust)

# Plot state-dependent distributions 
ggplot(behav.res %>% filter(behav %in% 1:max.gr), aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity', color = "black") +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(4, option = 'inferno'), rep("grey35", 6)),
                    guide = FALSE) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  scale_x_continuous(breaks = 1:10) +
  facet_grid(behav ~ var, scales = "free_x")




## Attribute behaviors to states and extract each of the different estimates
# Using MAP estimate, threshold of 90% assignments from posterior, and most common state
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
                ~case_when(. == 1 ~ "Slow-Unif",
                           . == 2 ~ "Exploratory",
                           . == 3 ~ "Slow-Turn",
                           . == 4 ~ "Transit",
                           is.na(.) ~ "Unclassified")
                )) %>% 
  mutate(across(c('z.map','z.post.thresh','z.post.max'),
                factor, levels = c('Slow-Turn','Slow-Unif','Exploratory','Transit',
                                   'Unclassified')
                ))


#by ID
ggplot(dat2, aes(easting, northing, color = id)) +
  geom_path(aes(group = id), color = "grey80", size = 0.5) +
  geom_point(size = 1) +
  theme_bw()

#by z.map
ggplot(dat2, aes(easting, northing, color = z.map)) +
  scale_color_manual("", values = c(viridis::viridis(4), "grey")) +
  geom_path(aes(group = id), color = "grey80", size = 0.5) +
  geom_point(size = 1) +
  theme_bw()

#by z.post.thresh
ggplot(dat2, aes(easting, northing, color = z.post.thresh)) +
  scale_color_manual("", values = c(viridis::viridis(4), "grey")) +
  geom_path(aes(group = id), color = "grey80", size = 0.5) +
  geom_point(size = 1) +
  theme_bw()

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
                ~case_when(. == 1 ~ "Slow-Unif",
                           . == 2 ~ "Exploratory",
                           . == 3 ~ "Slow-Turn",
                           . == 4 ~ "Transit",
                           is.na(.) ~ "Unclassified")
  ))

# write.csv(dat2, "Giant Armadillo state estimates.csv",row.names=F)
# write.csv(behav.res, "Giant Armadillo behavior distribs.csv", row.names=F)