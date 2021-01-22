# rm(list=ls())
# setwd('U:\\josh cullen\\giant armadillo\\derived data')
dat3=read.csv('GA edited data2.csv',as.is=T)

#discretize act.count
hist(dat3$act.count)
act.brk=seq(from=0,to=300,by=50)
range(dat3$act.count,na.rm=T)
act.brk[length(act.brk)]=300.1
tmp=cut(dat3$act.count,breaks=act.brk)
dat3$act.cat=as.numeric(tmp)
sum(!is.na(dat3$act.count) & is.na(dat3$act.cat))
boxplot(act.count~act.cat,data=dat3)

#discretize SL
hist(dat3$sl2)
mean(dat3$sl2>60)
sl.brk=c(seq(from=0,to=60,by=10),max(dat3$sl2)+1)
dat3$sl.cat=NA
for (i in 2:length(sl.brk)){
  cond=dat3$sl2 >= sl.brk[i-1] & dat3$sl2 < sl.brk[i]
  dat3$sl.cat[cond]=i-1
}
table(dat3$sl.cat)

#discretize TA
hist(dat3$ta2)
range(dat3$ta2,na.rm=T)
ta.brk=seq(from=-pi,to=pi,length.out=11)
ta.brk[length(ta.brk)]=3.142 #to include observation if exactly equal to pi
tmp=cut(dat3$ta2,breaks=ta.brk)
dat3$ta.cat=as.numeric(tmp)
table(dat3$ta.cat)
sum(!is.na(dat3$ta2) & is.na(dat3$ta.cat))
boxplot(ta2~ta.cat,data=dat3)

#export results
# setwd('U:\\josh cullen\\giant armadillo\\derived data')
write.csv(dat3[,c('id','sl.cat','ta.cat','act.cat','northing','easting','date_time')],
          'GA edited data3.csv', row.names=F)
