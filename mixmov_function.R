sample.z=function(nobs,nmaxclust, dat, ltheta, lphi, ndata.types){
  lprob=matrix(ltheta,nobs,nmaxclust,byrow=T)
  for (i in 1:nmaxclust){
    for (j in 1:ndata.types){
      #account for NA
      tmp=lphi[[j]][i,dat[,j]]
      cond=is.na(dat[,j])
      tmp[cond]=0
      
      #finish calculation
      lprob[,i]=lprob[,i]+tmp
    }
  }
  max1=apply(lprob,1,max)
  lprob=lprob-max1
  tmp=exp(lprob)
  prob=tmp/rowSums(tmp)
  
  z=rmultinom1(prob=prob, randu=runif(nobs)) 
  z+1
}
#-----------------------------------
sample.v=function(z,gamma1,nmaxclust){
  tmp=table(z)
  tmp1=rep(0,nmaxclust)
  tmp1[as.numeric(names(tmp))]=tmp
  theta=v=rep(NA,nmaxclust)
  aux=1
  for (i in 1:(nmaxclust-1)){
    seq1=(i+1):nmaxclust
    v[i]=rbeta(1,tmp1[i]+1,sum(tmp1[seq1])+gamma1)
    theta[i]=v[i]*aux
    aux=aux*(1-v[i])
  }
  v[nmaxclust]=1
  theta[nmaxclust]=aux
  list(theta=theta,v=v)
}
#-----------------------------------
sample.phi=function(alpha,nmaxclust,ncat.dat,ndata.types,nmat){
  phi=list()
  for (j in 1:ndata.types){
    tmp2=matrix(NA,nmaxclust,ncat.dat[j])
    for (i in 1:nmaxclust){
      tmp2[i,]=rdirichlet(1,nmat[[j]][i,]+alpha)
    }
    phi[[j]]=tmp2
  }
  phi
}
#------------------------------------
get.llk=function(phi,theta,ndata.types,dat,nobs,nmaxclust){
  prob=matrix(theta,nobs,nmaxclust,byrow=T)
  for (i in 1:nmaxclust){
    for (j in 1:ndata.types){
      #account for NA      
      tmp=phi[[j]][i,dat[,j]]
      cond=is.na(dat[,j])
      tmp[cond]=1
      
      #finish calculation
      prob[,i]=prob[,i]*tmp
    }
  }
  sum(log(rowSums(prob)))
}
#------------------------------------
sample.gamma=function(v,ngroup,gamma.possib){
  #calculate the log probability associated with each possible value of gamma
  cond=v>0.9999999
  v[cond]=0.9999999
  
  ngamma=length(gamma.possib)
  soma=sum(log(1-v[-ngroup]))
  k=(ngroup-1)*(lgamma(1+gamma.possib)-lgamma(gamma.possib))
  res=k+(gamma.possib-1)*soma
  #check this code: sum(dbeta(v[-ngroup],1,gamma.possib[5],log=T))
  
  #exponentiate and normalize probabilities
  res=res-max(res)
  res1=exp(res)
  res2=res1/sum(res1)
  
  #sample from a categorical distribution
  tmp=rmultinom(1,size=1,prob=res2)
  ind=which(tmp==1)
  gamma.possib[ind]
}