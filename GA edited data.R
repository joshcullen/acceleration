# rm(list=ls())
# setwd('U:\\josh cullen\\giant armadillo\\data')
files<- list.files(pattern = "*.csv")
files<- files[!grepl(pattern = "Armadillo", files)]
nfiles=length(files)

#get gps data
fim=numeric()
for (i in 2:nfiles){
  #import the data
  dat=read.csv(files[i],as.is=T)
  
  #eliminate NA's
  cond=!is.na(dat$GPS.UTM.Northing)
  dat1=dat[cond,]
  
  #just keep high quality GPS data (satellite count>X)
  print(table(dat1[,c('GPS.Fix.Attempt','GPS.Satellite.Count')]))
  cond=dat1$GPS.Satellite.Count>4
  dat2=dat1[cond,]
  
  #things to keep
  nomes=c('id','GPS.Fix.Time','GPS.UTM.Northing','GPS.UTM.Easting')
  dat4=dat2[,nomes]
  
  #store into a single file
  fim=rbind(fim,dat4)
}

# setwd('U:\\josh cullen\\giant armadillo\\derived data')
write.csv(fim,'GA edited data.csv',row.names=F)
