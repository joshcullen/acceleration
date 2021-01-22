# rm(list=ls(all=TRUE))
# setwd('U:\\josh cullen\\giant armadillo\\data')
dat=read.csv('accel_data_200720.csv',as.is=T)
table(dat$Predeployment.Data)
table(dat$bout)
table(dat$act)

boxplot(Activity.Count~act,data=dat)
boxplot(Activity.Count~bout,data=dat)
boxplot(Activity.Count~hora,data=dat)

#just predeployment=No
cond=dat$Predeployment.Data=='No'; mean(cond)
nomes=c('id','Acquisition.Time','Activity.Count')
dat1=dat[cond,nomes]

#export results
# setwd('U:\\josh cullen\\giant armadillo\\derived data')
write.csv(dat1,'accel_data.csv',row.names=F)
