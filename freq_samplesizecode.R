power_freq<-function(r.alloc,n,parE,theta)
{
  nP=r.alloc[1]*n # sample size in the placebo arm
  nR=r.alloc[2]*n # sample size in the reference arm
  nE=r.alloc[3]*n # sample size in the experimental arm
  
  parE_null<-parP + theta*(parR-parP) # value of parE under null.
  
  mu_null=parE_null-parP
  
  mu=parE-parP
  mv=parR-parP
  
  sigma2u_null<-parE_null*(1-parE_null)/nE + parP*(1-parP)/nP
  sigma2u<-parE*(1-parE)/nE + parP*(1-parP)/nP
  
  sigma2v<-parR*(1-parR)/nR + parP*(1-parP)/nP
  
  rho_null<-(parP*(1-parP)/nP)/(sqrt(sigma2u_null*sigma2v))
  rho<-(parP*(1-parP)/nP)/(sqrt(sigma2u*sigma2v))
  
  alpha<-(-mv/sqrt(sigma2v))
  
  phi=dnorm(alpha,0,1)
  
  c<-1-pnorm(alpha,0,1)
  
  E1_null<-mu_null+sqrt(sigma2u_null)*rho_null*phi/c
  
  E1<-mu+sqrt(sigma2u)*rho*phi/c
  
  E2<-mv+sqrt(sigma2v)*phi/c
  
  V1_null<-sigma2u_null*(1+((rho_null^2)*alpha*phi/c)-(rho_null*phi/c)^2)
  
  V1<-sigma2u*(1+((rho^2)*alpha*phi/c)-(rho*phi/c)^2)
  
  V2<-sigma2v*(1-(phi/c)*((phi/c)-alpha))
  
  E12_null<-(sqrt(sigma2u_null*sigma2v)*rho_null*(alpha*phi+c)/c) + (sqrt(sigma2u_null)*(mv*rho_null*phi/c)) +
    (sqrt(sigma2v)*(mu_null*phi/c)) + mu_null*mv
  
  E12<-(sqrt(sigma2u*sigma2v)*rho*(alpha*phi+c)/c) + (sqrt(sigma2u)*(mv*rho*phi/c)) +
    (sqrt(sigma2v)*(mu*phi/c)) + mu*mv
  
  cov_null<-E12_null-E1_null*E2
  
  cov<-E12-E1*E2
  
  mw_th_null<-E1_null-theta*E2
  
  mw_th_alt<-E1-theta*E2
  
  vw_null<-V1_null+(theta^2)*V2-2*theta*cov_null
  
  vw_alt<-V1+(theta^2)*V2-2*theta*cov
  
  z<-qnorm(0.975,0,1)
  
  k<-mw_th_null + z*sqrt(vw_null)
  
  power<-1 - pnorm(((k-mw_th_alt)/sqrt(vw_alt)),0,1)
  
  return(power)
}


# sample size calculation #
samplesize_fn_freq<-function(r.alloc,parE,theta)
{
  a<-1:1000
  
  power<-NULL
  
  for(i in a)
  {
    power<-c(power,power_freq(r.alloc=r.alloc,n=i,parE=parE,theta=theta))
  }
  n<-a[min(which(power>=0.8))]
  
  return(n)
}

### Sample size calculation for power=0.8 ###

parR<-.7
parP<-.1

p_star<-0.975

### sample size table ###
n.table.freq<-NULL
r.alloc.list<-matrix(c(1,1,1,1,2,2,1,2,3),3,3)

for (i.r.alloc in c(1:3))
{
  r.alloc<-r.alloc.list[i.r.alloc,]
  for (theta in c(0.8,0.7))
  {
    for (parE in c(.9, .85, .8, .75, .7, .65))
    {
      n<-samplesize_fn_freq(r.alloc=r.alloc,parE=parE,theta=theta)
      n.table.freq<-rbind(n.table.freq,c(r.alloc=r.alloc,parE=parE,theta=theta,n.freq=n))
    }
  }
}

save(n.table.freq,file="n.table.freq.RData")
(n.table.freq)
