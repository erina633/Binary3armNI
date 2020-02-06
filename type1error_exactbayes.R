#### fully bayesian power ####
power_fbayes<-function(r.alloc, n, parE, theta)
{
  set.seed(123)
  n_star<- 10000
  p_star<-0.975
  T_val = 1000
  parR<-.7
  parP<-.1
  
  Est_Prob = rep(0, n_star)
  
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
    
    
    while (count1<T_val)
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
        
        if (dEPn/dRPn<=((parE - parP)/(parR - parP)))
          count2=count2+1
      }
    }
    
    Est_Prob[i]=count2/T_val
    
    if(Est_Prob[i]>p_star)
      count=count+1
  }
  
  power<-count/n_star
  return(power)
} 

n1 = c(22, 32, 51, 94, 188)

r.allocmat = matrix(c(1,1,1,2,2,1,3, 2,1),3,3)

parE = c(.9, .85, .8, .75, .7)
theta = .8

fb_type1 = rep(NA, length(parE))
for(i in 1: length(parE)){
  fb_type1[i] = power_fbayes(r.allocmat[,1], n1, parE[i], theta)
}
fb_type1

