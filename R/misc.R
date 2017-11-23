## Convert cumulants to non-central moments
##    recursive formula produces as many cumulants as moments
##    References: Kenneth Lange: Numerical Analysis for Statisticians, 2nd ed. Page 15
cum2mnc <- function(kappa){
  N = length(kappa)+1
  mc = rep(0, N); mc[1] = 1
  for(k in 1:(N-1)){
    mc[k+1] = sum(choose(k-1, 0:(k-1))*kappa[k:1]*mc[1:k])
  }
  return(mc[-1])
}
## Convert non-central to central moments, uses recursive formula
##    optionally adjusts first moment to return mean
mnc2mc <- function(mnc){
  N = length(mnc)
  mc = rep(0,N); mc[1] = 0
  s1 = rep(c(1,-1), N)
  mnc = c(1,mnc)
  for(k in 1:(N-1)){
    mc[k+1] = sum( choose(k+1, 0:(k+1))*s1[(k+2):1]*mnc[1:(k+2)]*mnc[2]^((k+1):0) )
  }
  return(mc)
}
## non-central chisq \eqn{\chi_k^2(\delta)} cumulants 
## k: DF; lam: ncp 
chisq.cum <- function(k, delta, N){
  ik = 1:N
  2^(ik-1)*gamma(ik)*(k+ik*delta)
}
## 1-DF non-central chisq mix \eqn{\sum_k\lambda_k\chi_1^2(\delta_k)} cumulants
##   lam: mixing coef
##   delta: ncp
chisq1m.cum <- function(lam, delta, N){
  ## if(is.null(delta)) delta = rep(0,length(lam))
  ik = 1:N
  a1 = 2^(ik-1)*gamma(ik)
  cl = rep(0, N)
  for(k in 1:N) cl[k] = a1[k]*sum(lam^k*(1+k*delta))
  cl
}
## Liu non-central chisq approx
liu.params <- function(lam, delta){
  c1 = rep(0,4); for(i in 1:4){ c1[i] = sum(lam^i*(1+i*delta)) }
  muQ = c1[1];  sigmaQ = sqrt(2*c1[2])
  s1 = c1[3]/c1[2]^(3/2);  s2 = c1[4]/c1[2]^2
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2));  d = s1 *a^3 - a^2;  l = a^2 - 2*d
  } else {
    l = 1/s1^2;  a = sqrt(l);  d = 0
  }
  muX = l+d;  sigmaX = sqrt(2)*a
  list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
}
liu.pval1 <- function(Q.all, lambda,delta){
  param = liu.params(lambda,delta)
  Q.Norm = (Q.all - param$muQ)/param$sigmaQ
  Q.Norm1 = Q.Norm*param$sigmaX + param$muX
  pchisq(Q.Norm1, df = param$l,ncp=param$d, lower.tail=FALSE)
}
## Lee non-central chisq approx
lee.params <- function(lam,delta){
  c1 = rep(0,4); for(i in 1:4){ c1[i] = sum(lam^i*(1+i*delta)) }
  muQ = c1[1];  sigmaQ = sqrt(2 *c1[2])
  s1 = c1[3]/c1[2]^(3/2);  s2 = c1[4]/c1[2]^2
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2));  d = s1 *a^3 - a^2;  l = a^2 - 2*d
  } else {
    l = 1/s2;  a = sqrt(l);  d = 0
  }
  muX = l+d;  sigmaX = sqrt(2)*a
  list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
}
lee.pval1 <- function(Q.all, lambda,delta){
  param = lee.params(lambda,delta)
  Q.Norm = (Q.all - param$muQ)/param$sigmaQ
  Q.Norm1 = Q.Norm*param$sigmaX + param$muX
  pchisq(Q.Norm1, df = param$l,ncp=param$d, lower.tail=FALSE)
}
lee.qval <- function(pval, lambda,delta){
  param = lee.params(lambda,delta)
  q0  = qchisq(pval,df=param$l,ncp=param$d, lower.tail=FALSE)
  q1 = (q0-param$muX)/param$sigmaX*param$sigmaQ+param$muQ
  return(q1)
}
## Wu non-central chisq approx
wu.pval1 <- function(Q.all, lambda,delta, N=12){
  param = wu.params(lambda,delta,N)
  Q.Norm = (Q.all - param$muQ)/param$sigmaQ
  Q.Norm1 = Q.Norm*param$sigmaX + param$muX
  pchisq(Q.Norm1, df = param$l,ncp=param$d, lower.tail=FALSE)
}

## saddlepoint pval
saddle <- function(x,lambda,lower.tail=FALSE){
  d = max(lambda)
  lambda = lambda/d
  x = x/d
  k0 = function(zeta) -sum(log(1-2*zeta*lambda))/2
  kprime0 = function(zeta) sapply(zeta, function(zz) sum(lambda/(1-2*zz*lambda)))
  kpprime0 = function(zeta) 2*sum(lambda^2/(1-2*zeta*lambda)^2)
  n = length(lambda)
  if (any(lambda < 0)) {
    lmin = max(1/(2 * lambda[lambda < 0])) * 0.99999
  } else if (x>sum(lambda)){
    lmin = -0.01
  } else {
    lmin = -length(lambda)/(2*x)
  }
  lmax = min(1/(2*lambda[lambda>0]))*0.99999
  hatzeta = uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, upper = lmax, tol=1e-12,maxiter=1e4)$root
  w = sign(hatzeta)*sqrt(2*(hatzeta*x-k0(hatzeta)))
  v = hatzeta*sqrt(kpprime0(hatzeta))
  if(abs(hatzeta)<1e-4){
    return(NA)
  } else{
    return( pnorm(w+log(v/w)/w, lower.tail=lower.tail) ) ## surv func
  }
}
Sadd.pval <- function(Q.all,lambda,lower.tail=FALSE){
  sad = rep(1*(!lower.tail),length(Q.all))
  if(sum(Q.all>0)>0){
    sad[Q.all>0] = sapply(Q.all[Q.all>0],saddle,lambda=lambda,lower.tail=lower.tail)
  }
  id = which(is.na(sad))
  if(length(id)>0){
    sad[id] = lee.pval1(Q.all[id], lambda,0)
  }
  return(sad)
}
## SKAT p-value calc
##   mix of 1-DF central chisq
KAT.pval <- function(Q.all, lambda, acc=1e-12,lim=1e8){
  pval = rep(0, length(Q.all))
  i1 = which(is.finite(Q.all))
  for(i in i1){
    tmp = davies(Q.all[i],lambda,acc=acc,lim=lim); pval[i] = tmp$Qq
    if((tmp$ifault>0)|(pval[i]<0)|(pval[i]>1)){
      pval[i] = Sadd.pval(Q.all[i],lambda)
    }
  }
  return(pval)
}
## KAT tail prob calc
##  mix of 1-DF non-central chi-sq
KAT.tpr <- function(Q.all, lambda, delta, acc=1e-12,lim=1e8){
  Kdf = rep(1, length(lambda))
  pr = rep(0, length(Q.all))
  i1 = which(is.finite(Q.all))
  for(i in i1){
    tmp = davies(Q.all[i],lambda,Kdf,delta, acc=acc,lim=lim); pr[i] = tmp$Qq
    if((tmp$ifault>0)|(pr[i]<0)|(pr[i]>1)){
      pr[i] = wu.pval1(Q.all[i],lambda,delta,N=6)
    }
  }
  return(pr)
}
KAT.qdf <- function(pval, lam, lim=1e8,acc=1e-12){
  f0 = function(x) KAT.pval(x,lam, lim=lim,acc=acc) - pval
  q0 = lee.qval(pval, lam,0)
  a1 = q0*0.9; a2 = q0*1.5
  it = 1
  while( (it<20)&(f0(a1)*f0(a2)>0) ){
    a1 = a1*0.9; a2 = a2*1.5
    it = it+1
  }
  q1 = try(uniroot(f0, c(a1,a2), tol=pval*1e-5)$root)
  if(class(q1)=='try-error') return(q0)
  return(q1)
}
