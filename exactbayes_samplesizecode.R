
#### fully bayesian power ####

# internal function#

power_fbayes<-function(r.alloc,n,parE,theta)
{
  set.seed(123)
  
  Est_Prob<-rep(0,n_star)
  
  nP=r.alloc[1]*n
  nR=r.alloc[2]*n
  nE=r.alloc[3]*n
  
  tE<-tR<-tP<-1
  
  count<-0
  
  ########## Test Procedure ##############
  for(i in 1:n_star)
  {
    count1=0
    count2=0
    
    xP=rbinom(1, nP, parP)
    xR=rbinom(1, nR, parR)
    xE=rbinom(1, nE, parE)
    
    while (count1<T)
    {
      aE<-1; bE<-1; aR<-1; bR<-1; aP<-1; bP<-1
      parEgdata=rbeta(1,aE+xE,bE+nE-xE)
      parRgdata=rbeta(1,aR+xR,bR+nR-xR)
      parPgdata=rbeta(1,aP+xP,bP+nP-xP)
      
      dEPn=parEgdata-parPgdata
      dRPn=parRgdata-parPgdata
      
      if(dRPn>0)
      { 
        count1=count1+1
        
        if (dEPn/dRPn>theta)
          count2=count2+1
      }
    }
    
    Est_Prob[i]=count2/T
    
    if(Est_Prob[i]>p_star)
      count=count+1
  }
  
  power<-count/n_star
  return(power)
} 
# sample size calculation #
samplesize_fn_fbayes<-function(r.alloc,parE,theta,a_max)
{
  a_min<-max(1,(a_max-400))
  a_max<- a_max+100
  a<-a_min:a_max
  
  power<-NULL
  
  for(i in a)
  {
    power<-c(power,power_fbayes(r.alloc=r.alloc,n=i,parE=parE,theta=theta))
  }
  n<-try(a[min(which(power>=0.8))],silent=T)
  if(class(n)=="try-error") n<-NA
  
  return(n)
}


############## full Bayes ###################
n_star<-1000 ### No. of simulations ###
T<-1000 ### samples ###

p_star<-0.975

parR<-.7
parP<-.1

r.alloc.list<-matrix(c(1,1,1,1,2,2,1,2,3),3,3)
load("n.table.freq.RData")
fullbayes.table<-NULL
counter=1
for (i.r.alloc in c(1:1))
{
  r.alloc<-r.alloc.list[i.r.alloc,]
  for (theta in c(0.8))
  {
    for (parE in c(.9, .85, .8, .75, .7))
    {
      n<-samplesize_fn_fbayes(r.alloc=r.alloc,parE=parE,theta=theta,a_max=n.table.freq[counter,6])
      counter=counter+1 
      fullbayes.table<-rbind(fullbayes.table,c(r.alloc=r.alloc,parE=parE,theta=theta,n=n))
    }
  }
}

fullbayes.table