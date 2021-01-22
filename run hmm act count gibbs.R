rm(list=ls(all=TRUE))
library('MCMCpack')
library('Rcpp')
set.seed(3)

setwd('U:\\GIT_models\\mixture_movem')
source('mixmov_function.R')
source('mixmov_gibbs.R')
sourceCpp('aux1.cpp')

#import the data
setwd('U:\\josh cullen\\giant armadillo\\derived data')
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
seq1=1:ngibbs
seq1=nburn:ngibbs
plot(mod$loglikel[seq1], type='l')

setwd('U:\\josh cullen\\giant armadillo\\hmm results act count')
write.csv(mod$phi[[1]][seq1,],'phi1.csv',row.names=F)
write.csv(mod$phi[[2]][seq1,],'phi2.csv',row.names=F)
write.csv(mod$phi[[3]][seq1,],'phi3.csv',row.names=F)
write.csv(mod$theta[seq1,],'theta.csv',row.names=F)
write.csv(mod$loglikel,'llk.csv',row.names=F)
write.csv(mod$z,'z.csv',row.names=F)
write.csv(mod$gamma1[seq1],'gamma1.csv',row.names=F)