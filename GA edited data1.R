# rm(list=ls())
# setwd('U:\\josh cullen\\giant armadillo\\derived data')
dat=read.csv('GA edited data.csv',as.is=T)
accel=read.csv('accel_data.csv',as.is=T)

#for each individual, calculate SL, angle, and time interval
uni=unique(dat$id)
nuni=length(uni)

fim2=numeric()
for (i in 1:nuni){
  #get accel data for individual i
  cond=accel$id==uni[i]
  accel1=accel[cond,]
  
  #get GPS data for individual i  
  cond=dat$id==uni[i]
  dat1=dat[cond,]
  nobs=nrow(dat1)
  
  #get time
  dat1$time1=as.POSIXct(strptime(dat1$GPS.Fix.Time, 
                                 format = "%Y.%m.%d %H:%M:%S", 
                                 tz = "UTC"))
  accel1$time1=as.POSIXct(strptime(accel1$Acquisition.Time, 
                                   format = "%Y.%m.%d %H:%M:%S", 
                                   tz = "UTC"))
  
  #order the data chronologically
  ind=order(dat1$time1)
  unique(ind-c(1:length(ind)))
  dat2=dat1[ind,]
  
  fim=matrix(NA,nobs,10)
  for (j in 3:nobs){
    #get activity count within less than 1 min of end time of SL
    dt1=dat2$time1[j]-accel1$time1
    attr1=attributes(dt1)$units
    ind=numeric()
    if (attr1=='secs') ind=which(abs(dt1)<=60)
    if (attr1=='min') ind=which(abs(dt1)<=1)
    count2=NA
    if (length(ind)==1) count2=accel1$Activity.Count[ind]
    if (length(ind) >1) count2=accel1$Activity.Count[sample(ind,size=1)]
    
    #get time
    time.int1=dat2$time1[j-1]-dat2$time1[j-2]
    time1=as.numeric(time.int1)
    unit1=attributes(time.int1)$units
    time.int2=dat2$time1[j]-dat2$time1[j-1]
    time2=as.numeric(time.int2)
    unit2=attributes(time.int2)$units

    #get dx's and dy's
    dx1=(dat2$GPS.UTM.Easting[j-1]-dat1$GPS.UTM.Easting[j-2])
    dy1=(dat2$GPS.UTM.Northing[j-1]-dat2$GPS.UTM.Northing[j-2])
    dx2=(dat2$GPS.UTM.Easting[j]-dat1$GPS.UTM.Easting[j-1])
    dy2=(dat2$GPS.UTM.Northing[j]-dat2$GPS.UTM.Northing[j-1])
    
    #get step length
    sl2=sqrt((dx2^2)+(dy2^2))
    
    #get angles 1 and 2
    angle1=atan2(dy1,dx1)
    angle2=atan2(dy2,dx2)
    tmp=angle2-angle1
    
    #wrap this around if angle is bigger than pi or smaller than -pi
    tmp=ifelse(tmp>pi,tmp-2*pi,tmp)
    TA2=ifelse(tmp< -pi, tmp+2*pi,tmp)

    #store results
    fim[j,]=c(sl2,TA2,time1,unit1,time2,unit2,count2,
              dat2$GPS.UTM.Northing[j],dat2$GPS.UTM.Easting[j],
              dat2$time1[j])
  }
  colnames(fim)=c('sl2','ta2','time1','unit1','time2','unit2','act.count','northing','easting','date_time')
  fim2=rbind(fim2,cbind(uni[i],fim))
}
colnames(fim2)[1]='id'
mean(is.na(fim2[,'act.count']))

# setwd('U:\\josh cullen\\giant armadillo\\derived data')
write.csv(fim2,'GA edited data1.csv',row.names=F)
