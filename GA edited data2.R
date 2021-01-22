# rm(list=ls())
# setwd('U:\\josh cullen\\giant armadillo\\derived data')
dat=read.csv('GA edited data1.csv',as.is=T)

table(dat$unit1)
table(dat$unit2)
cond=dat$unit1=='hours'
hist(dat$time1[cond])
cond=dat$unit1=='secs'
hist(dat$time1[cond])

#just keep those observations in which the difference in fixes is measured in minutes
cond=!is.na(dat$unit2) & dat$unit2=='mins'
dat1=dat[cond,]; mean(cond,na.rm=T)

hist(dat1$time2)
rtime2=round(dat1$time2)
table(rtime2)

#just keep measurements between 6 and 8 minutes
cond=!is.na(dat1$time2) & dat1$time2>= 6 & dat1$time2 <= 8
dat2=dat1[cond,]; mean(cond)

#if time1 is not between 6 and 8, make TA=NA
cond=!is.na(dat2$time1) & dat2$time1>= 6 & dat2$time1 <= 8
dat2$ta2[!cond]=NA

mean(is.na(dat2$ta2)); mean(is.na(dat2$sl2))

#eliminate extreme SL's
hist(dat2$sl2)
qtt=quantile(dat2$sl2,0.99) #cut-off of 800
hist(dat2$sl2[dat2$sl2<qtt])
cond=dat2$sl2< 800
dat3=dat2[cond,]

#normalize SL by dividing by time
dat3$sl2=dat3$sl2/dat3$time2

#relationship between sl2 and act.count
teste=dat3
plot(act.count~sl2,data=teste) 
teste$denis=NA
cond=!is.na(teste$act.count) & teste$act.count==0
teste$denis[cond]=0
cond=!is.na(teste$act.count) & teste$act.count>0
teste$denis[cond]=1
boxplot(sl2~denis,data=teste)

#if act.count==0, make it equal to NA (not sure if mal-functioning or not). For example, might have fallen to the ground
cond=!is.na(dat3$act.count) & dat3$act.count==0
dat3$act.count[cond]=NA

#save results
# setwd('U:\\josh cullen\\giant armadillo\\derived data')
write.csv(dat3[,c('id','sl2','ta2','act.count','northing','easting','date_time')],
          'GA edited data2.csv', row.names=F)
