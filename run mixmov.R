# rm(list=ls(all=TRUE))
library('MCMCpack')
library('Rcpp')
set.seed(2)

setwd('U:\\GIT_models\\mixture_movem')
source('mixmov_function.R')
source('mixmov_gibbs.R')
sourceCpp('aux1.cpp')
dat=read.csv('fake data.csv',as.is=T)

#prior
alpha=0.1

#prepare for gibbs
ngibbs=10000
nburn=ngibbs/2
nmaxclust=15
model1=mixture_movement(dat=dat,alpha=alpha,ngibbs=ngibbs,nmaxclust=nmaxclust,nburn=nburn)


