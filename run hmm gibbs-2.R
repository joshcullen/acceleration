library('MCMCpack')
library('Rcpp')
library(tidyverse)
set.seed(3)

source('mixmov_function.R')
source('mixmov_gibbs.R')
sourceCpp('aux1.cpp')

#import the data
dat0=read.csv('GA edited data3.csv',as.is=T)
dat=data.matrix(dat0[,c('sl.cat','ta.cat','act.cat')])

#prior
alpha=0.1

#initialize parameters
nmaxclust=10

#MCMC stuff
ngibbs=20000
nburn=ngibbs/2

#run gibbs
mod=mixture_movement(dat=dat,alpha=alpha,ngibbs=ngibbs,nmaxclust=nmaxclust,nburn=nburn)

#look at convergence
post.seq=(nburn + 1):ngibbs
plot(mod$loglikel[post.seq], type='l')


## Inspect and Plot results
# MAP.iter<- nburn + which.max(mod$loglikel[seq1])

theta<- mod$theta[post.seq,]
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
  plot(matrix(mod$phi[[1]][i,], nrow = nmaxclust)[4,], type = "h")
  plot(matrix(mod$phi[[2]][i,], nrow = nmaxclust)[4,], type = "h")
  plot(matrix(mod$phi[[3]][i,], nrow = nmaxclust)[4,], type = "h")
}
## Checks out, no state flipping when plotting random samples from posterior


sl<- mod$phi[[1]][post.seq,]
sl<- matrix(colMeans(sl), nrow = nmaxclust)
rowSums(sl)  #check that data stored properly
sl<- sl[1:4,]

ta<- mod$phi[[2]][post.seq,]
ta<- matrix(colMeans(ta), nrow = nmaxclust)
rowSums(ta)  #check that data stored properly
ta<- ta[1:4,]

act<- mod$phi[[3]][post.seq,]
act<- matrix(colMeans(act), nrow = nmaxclust)
rowSums(act)  #check that data stored properly
act<- act[1:4,]

par(mfrow = c(4,3))
for (i in 1:nrow(sl)) {
  plot(sl[i,], type = "h")
  plot(ta[i,], type = "h")
  plot(act[i,], type = "h")
}

# Add to single data frame
sl.df<- data.frame(sl, state = 1:nrow(sl))
names(sl.df)<- gsub('X', '', names(sl.df))
sl.long<- pivot_longer(sl.df, cols = -state, names_to = "bin", values_to = "value")
sl.long$var<- "Speed"

ta.df<- data.frame(ta, state = 1:nrow(ta))
names(ta.df)<- gsub('X', '', names(ta.df))
ta.long<- pivot_longer(ta.df, cols = -state, names_to = "bin", values_to = "value")
ta.long$var<- "Turning Angle"

act.df<- data.frame(act, state = 1:nrow(act))
names(act.df)<- gsub('X', '', names(act.df))
act.long<- pivot_longer(act.df, cols = -state, names_to = "bin", values_to = "value")
act.long$var<- "Activity Count"

phi.df<- rbind(act.long, sl.long, ta.long)



## Determine likely state for each observation
# Using threshold of 0.75% assignments from posterior
z.post<- as.matrix(mod$z.posterior)
z.post2<- t(apply(z.post, 1, function(x) x/sum(x)))
thresh<- 0.75
z.post3<- apply(z.post2, 1, function(x) ifelse(max(x) > thresh, which(x > thresh), NA))
z.post4<- apply(z.post2, 1, function(x) which.max(x))

# Compare posterior estimates against MAP
table(z.post3, mod$z.MAP)  #estimates are close, but some discrepancies
table(z.post4, mod$z.MAP)  #estimates are close, but some discrepancies


## Map estimated states
dat0<- dat0 %>% 
  mutate(z.map = mod$z.MAP,
         z.post.thresh = z.post3,
         z.post.max = z.post4)
dat0$z.map<- ifelse(dat0$z.map > max.gr, NA, dat0$z.map)
dat0$z.post.max<- ifelse(dat0$z.post.max > max.gr, NA, dat0$z.post.max)

#by ID
ggplot(dat0, aes(easting, northing, color = id)) +
  geom_path(aes(group = id), color = "grey80", size = 0.5) +
  geom_point(size = 1) +
  theme_bw()

#by z.map
ggplot(dat0, aes(easting, northing, color = factor(z.map))) +
  geom_path(aes(group = id), color = "grey80", size = 0.5) +
  geom_point(size = 1) +
  theme_bw()

#by z.post.thresh
ggplot(dat0, aes(easting, northing, color = factor(z.post.thresh))) +
  geom_path(aes(group = id), color = "grey80", size = 0.5) +
  geom_point(size = 1) +
  theme_bw()

#by z.post.max
ggplot(dat0, aes(easting, northing, color = factor(z.post.max))) +
  geom_path(aes(group = id), color = "grey80", size = 0.5) +
  geom_point(size = 1) +
  theme_bw()

### Over 1/3 of all observations have an "Unclassified" state based on posterior estimates
### Only 0.7% of all observations are "Unclassified" based on MAP estimate
### Could either adjust threshold for state classification or deal with messy results


# write.csv(dat0, "Giant Armadillo state estimates.csv",row.names=F)
# write.csv(phi.df,"Giant Armadillo behavior distribs.csv",row.names=F)
