#### Approximate bayesian power ####

# internal function #
power_approxbayes<-function(r.alloc=r.alloc,n=i,parE=parE,theta=theta)
{
  nP=r.alloc[1]*n
  nR=r.alloc[2]*n
  nE=r.alloc[3]*n
  aE<-1; bE<-1; aR<-1; bR<-1; aP<-1; bP<-1
  me=aE/(aE+bE)
  mr=aR/(aR+bR)
  mp=aP/(aP+bP)
  
  sigma2e=(aE*bE)/(((aE+bE)^2)*(aE+bE+1))
  sigma2r=(aR*bR)/(((aR+bR)^2)*(aR+bR+1))
  sigma2p=(aP*bP)/(((aP+bP)^2)*(aP+bP+1))
  
  sigma2_T<-parE*(1-parE)/nE + (theta^2)*parR*(1-parR)/nR + ((1-theta)^2)*parP*(1-parP)/nP 
  
  
  m_etaEP=me-mp
  m_etaRP=mr-mp
  
  sigma2_etaEP<-sigma2e + sigma2p
  
  sigma2_etaRP<-sigma2r + sigma2p 
  
  rho<-sigma2p/(sqrt(sigma2_etaRP*sigma2_etaEP)) 
  
  alpha<-(-m_etaRP/sqrt(sigma2_etaRP))
  
  phi=dnorm(alpha,0,1)
  
  c<-1-pnorm(alpha,0,1)
  
  delta_alpha=(phi/c)*((phi/c)-alpha)
  
  E1<-m_etaEP+sqrt(sigma2_etaEP)*rho*phi/c
  E2<-m_etaRP+sqrt(sigma2_etaRP)*phi/c
  
  V1<-sigma2_etaEP*(1+(rho^2*alpha*phi/c)-(rho*phi/c)^2)
  V2<-sigma2_etaRP*(1-delta_alpha)
  
  E12<-(sqrt(sigma2_etaEP*sigma2_etaRP)*rho*(alpha*phi+c)/c) + (sqrt(sigma2_etaEP)*m_etaRP*rho*phi/c) +
    (sqrt(sigma2_etaRP)*m_etaEP*phi/c) + m_etaEP*m_etaRP
  
  cov<-E12-E1*E2
  
  mu_star_th<-E1-theta*E2
  
  sigma2_star_nu<-V1+(theta^2)*V2-2*theta*cov
  
  mu_T_theta<-parE - theta*parR -(1-theta)*parP
  
  z<-qnorm(1-p_star,0,1)
  
  point<-z*sqrt((1/sigma2_T) + (1/sigma2_star_nu))*sqrt(sigma2_T) + 
    (mu_star_th/sigma2_star_nu)*(sqrt(sigma2_T)) + (mu_T_theta/sqrt(sigma2_T))
  
  power<-pnorm(point,0,1)
  
  return(power)
}

# sample size calculation #
samplesize_fn_approxbayes<-function(r.alloc,parE,theta)
{
  a<-1:1000
  
  power<-NULL
  
  for(i in a)
  {
    power<-c(power,power_approxbayes(r.alloc=r.alloc,n=i,parE=parE,theta=theta))
  }
  power<-round(power,3)
  n<-a[min(which(power>=0.8))]
  
  return(n)
}

#################### Approx Bayes ###################

p_star<-0.975
parR<-.7
parP<-.1

n.table.ab<-NULL
r.alloc.list<-matrix(c(1,1,1,1,2,2,1,2,3),3,3)

for (i.r.alloc in c(1:3))
{
  r.alloc<-r.alloc.list[i.r.alloc,]
  for (theta in c(0.8,0.7))
  {
    for (parE in c(.9, .85, .8, .75, .7))
    {
      n<-samplesize_fn_approxbayes(r.alloc=r.alloc,parE=parE,theta=theta)
      n.table.ab<-rbind(n.table.ab,c(r.alloc=r.alloc,parE=parE,theta=theta,n.ab=n))
    }
  }
}

n.table.ab
