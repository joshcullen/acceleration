rm(list=ls(all=TRUE))
nmaxclust=10

#get categories
setwd('U:\\josh cullen\\giant armadillo\\derived data')
sl.cat=read.csv('SL categories.csv',as.is=T)$x
sl.cat1=sl.cat[-length(sl.cat)]+(diff(sl.cat)[1]/2)
sl.cat1[length(sl.cat1)]='>60'

ta.cat=read.csv('TA categories.csv',as.is=T)$x
ta.cat1=ta.cat[-length(ta.cat)]+(diff(ta.cat)[1]/2)
ta.cat1=round(ta.cat1,2)

act.cat=read.csv('ACT categories.csv',as.is=T)$x
act.cat1=act.cat[-length(act.cat)]+(diff(act.cat)[1]/2)
act.cat1=round(act.cat1,2)

#how many groups?
setwd('U:\\josh cullen\\giant armadillo\\hmm results act count')
theta=read.csv('theta.csv',as.is=T)
theta1=colMeans(theta)
plot(theta1,type='h')

max.gr=4
mean(rowSums(theta[,1:max.gr])) #97%
median(rowSums(theta[,1:max.gr])) #97%
range(rowSums(theta[,1:max.gr])) #93%-100%

#get phi1
phi1=read.csv('phi1.csv',as.is=T)
ncat1=ncol(phi1)/nmaxclust; ncat1
phi1a=matrix(colMeans(phi1),nmaxclust,ncat1)
phi1b=phi1a[1:max.gr,]

#get phi2
phi2=read.csv('phi2.csv',as.is=T)
ncat2=ncol(phi2)/nmaxclust; ncat2
phi2a=matrix(colMeans(phi2),nmaxclust,ncat2)
phi2b=phi2a[1:max.gr,]

#get phi3
phi3=read.csv('phi3.csv',as.is=T)
ncat3=ncol(phi3)/nmaxclust; ncat3
phi3a=matrix(colMeans(phi3),nmaxclust,ncat3)
phi3b=phi3a[1:max.gr,]

#plot results
setwd('U:\\josh cullen\\giant armadillo\\hmm results act count')
png('summary results.png',width=700,height=1000)
par(mfcol=c(max.gr,3),mar=c(4,3,1,1),oma=c(1,4,4,1))
ylim1=range(phi1b[1:max.gr,])
for (i in 1:max.gr){
  plot(phi1b[i,],type='h',ylim=ylim1,xaxt='n',cex.axis=2,xlab='')
  text(x=7,y=ylim1[2]-0.1,i,cex=3)
  axis(side=1,at=1:7,sl.cat1,las=2,cex.axis=2,xlab='')
}

ylim1=range(phi2b[1:max.gr,])
for (i in 1:max.gr){
  plot(phi2b[i,],type='h',ylim=ylim1,xaxt='n',cex.axis=2,xlab='')
  text(x=9,y=ylim1[2]-0.05,i,cex=3)
  axis(side=1,at=1:10,ta.cat1,las=2,cex.axis=2,xlab='')
}

ylim1=range(phi3b[1:max.gr,])
for (i in 1:max.gr){
  plot(phi3b[i,],type='h',ylim=ylim1,xaxt='n',cex.axis=2,xlab='')
  text(x=5,y=ylim1[2]-0.05,i,cex=3)
  axis(side=1,at=1:length(act.cat1),act.cat1,las=2,cex.axis=2,xlab='')
}
mtext(side=2,at=0.5,outer=T,line=1,'Probability',cex=3)
mtext(side=3,at=c(0.18,0.5,0.85),outer=T,line=1,
      c('Step length','Turning Angle','Act. count'),cex=2)
dev.off()

#it is easier to compare SL's side by side
png('summary results SL.png',width=700,height=700)
par(mfcol=c(1,1),mar=c(4,4,1,1))
ylim1=range(phi1b[1:max.gr,])
ncat=ncol(phi1b)
plot(NA,xlim=c(1,ncat+0.4),ylim=ylim1,xlab='',xaxt='n',cex.axis=2,ylab='')
for (i in 1:max.gr){
  for (j in 1:ncat){
    lines(rep(j,2)+((i-1)/10),c(0,phi1b[i,j]),col=i)    
  }
}
legend(x=4,y=0.8,lty=1,col=1:4,paste('State',1:4),cex=2)
axis(side=1,at=1:7,sl.cat1,las=2,cex.axis=2,xlab='')
dev.off()