#' Non-central chi-square approximation to mixture of 1-DF chi-square distributions
#'
#' Given the mixture of 1-DF chi-square distribution, \eqn{\sum_i\lambda_i\chi_1^2(\delta_i)},
#' the non-central chi-square distribution approximation is based on matching higher moments. See the ref.
#' @param  lambda  mixing coefficients for chi-square mixtures
#' @param  delta  vector of non-centrality parameters for chi-square mixtures
#' @param  N  match (N-1,N)-th moments
#' @return
#' \describe{
#'   \item{l}{ estimated DF }
#'   \item{d}{ estimated non-centrality parameter }
#'   \item{muQ,sigmaQ}{ first two moments based on the chi-square mixture distribution}
#'   \item{muX,sigmaX}{ first two moments based on the estimated non-central chi-square distribution}
#' }
#' 
#' @export
#' @references
#' Wu,B. and Pankow,J.S. (2016) On sample size and power calculation for variant set-based association tests. /AHG/, 80(2), 136-143.
#'
#' Wu,B. and Pankow,J.S. (2017) On computing the tail probability of non-negative definite quadratic forms in central normal variables. tech rep.
wu.params <- function(lambda,delta, N=12){
  cl = chisq1m.cum(lambda,delta, N)
  a0 = cum2mnc(cl)
  muQ = a0[1]; sigmaQ = sqrt(a0[2]-a0[1]^2)
  a1 = mnc2mc(a0)
  a1 = a1/sqrt(a1[2])^(1:N)  
  f1 = function(xpar){
    k = exp(xpar[1])
    v = xpar[2]
    a2 = mnc2mc(cum2mnc(chisq.cum(k,v,N)))
    a2 = a2/sqrt(a2[2])^(1:N)  
    (a1[N-1]-a2[N-1])^2 + (a1[N]-a2[N])^2
  }
  tmp = bobyqa(c(0,1), f1, lower=c(-Inf,0),upper=c(Inf,Inf))
  xpar = tmp$par
  l = exp(xpar[1])
  d = xpar[2]
  if(f1(c(xpar[1],0))<=tmp$fval){
    d=0
    f.1 = function(xpar) f1(c(xpar,0))
    l = exp(bobyqa(xpar[1], f.1)$par)
  }
  muX = l+d; sigmaX = sqrt(2*l+4*d) ## sqrt(chisq.cum(l,d,N=2)[2])
  list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
}


#' Compute power for SKAT (given chi-square mixture parameters) using non-central chi-square approximation
#'
#' Given the 1-DF chi-square mixture distribution, \eqn{\sum_i\lambda_i\chi_1^2(\delta_i)},
#' the non-central chi-square approximation is based on moment-matching (typically higher moment). See the reference.
#' @param  alpha significance level
#' @param  lambda  vector of mixing coefficients
#' @param  delta  vector of non-centrality parameters
#' @param  N0  match (N0-1,N0)-th moments for null p-value calculation
#' @param  N1  match (N1-1,N1)-th moments for power calculation
#' @return
#' \describe{
#'   \item{pwr}{ estimated power }
#' }
#' 
#' @export
#' @references
#' Wu,B. and Pankow,J.S. (2016) On sample size and power calculation for variant set-based association tests. /AHG/, 80(2), 136-143.
#'
#' Wu,B. and Pankow,J.S. (2017) On computing the tail probability of non-negative definite quadratic forms in central normal variables. tech rep.
wu.pwr <- function(alpha, lambda,delta,N0=12,N1=6){
  lampar = wu.params(lambda,delta=0,N0)
  lampar1 = wu.params(lambda,delta,N1)
  q0 = qchisq(alpha,df=lampar$l,ncp=lampar$d,lower=FALSE)
  Q0 = (q0-lampar$muX)/lampar$sigmaX*lampar$sigmaQ+lampar$muQ
  Q1 = (Q0-lampar1$muQ)/lampar1$sigmaQ*lampar1$sigmaX+lampar1$muX
  pwr = pchisq(Q1,df=lampar1$l,ncp=lampar1$d,lower=FALSE)
  return(pwr)
}
#' Compute power for SKAT (given chi-square mixture parameters) using non-central chi-square approximation
#'
#' The non-central chi-square approximation is based on the Liu method. 
#' @param  alpha significance level
#' @param  lambda  vector of mixing coefficients
#' @param  delta  vector of non-centrality parameters
#' @return estimated power
#' @export
#' @references
#' Wu,B. and Pankow,J.S. (2016) On sample size and power calculation for variant set-based association tests. /AHG/, 80(2), 136-143.
#'
#' Wu,B. and Pankow,J.S. (2017) On computing the tail probability of non-negative definite quadratic forms in central normal variables. tech rep.
#'
#' H. Liu, Y. Tang, H.H. Zhang (2009) A new chi-square approximation to the distribution of non-negative definite quadratic forms in non-central normal variables, Computational Statistics and Data Analysis, 53, 853-856.
liu.pwr <- function(alpha, lambda,delta){
  lampar = liu.params(lambda,delta=0)
  lampar1 = liu.params(lambda,delta)
  q0 = qchisq(alpha,df=lampar$l,ncp=lampar$d,lower=FALSE)
  Q0 = (q0-lampar$muX)/lampar$sigmaX*lampar$sigmaQ+lampar$muQ
  Q1 = (Q0-lampar1$muQ)/lampar1$sigmaQ*lampar1$sigmaX+lampar1$muX
  pchisq(Q1,df=lampar1$l,ncp=lampar1$d,lower=FALSE)
}
#' Compute power for SKAT (given chi-square mixture parameters) using non-central chi-square approximation
#'
#' The non-central chi-square approximation is based on a modified Liu method. 
#' @param  alpha significance level
#' @param  lambda  vector of mixing coefficients
#' @param  delta  vector of non-centrality parameters
#' @return estimated power
#' @export
#' @references
#' Wu,B. and Pankow,J.S. (2016) On sample size and power calculation for variant set-based association tests. /AHG/, 80(2), 136-143.
#'
#' Wu,B. and Pankow,J.S. (2017) On computing the tail probability of non-negative definite quadratic forms in central normal variables. tech rep.
#'
#' Lee, S., Wu, M. C., and Lin, X. (2012) Optimal tests for rare variant effects in sequencing association studies. Biostatistics, 13, 762-775. 
lee.pwr <- function(alpha, lambda,delta){
  lampar = lee.params(lambda,delta=0)
  lampar1 = lee.params(lambda,delta)
  q0 = qchisq(alpha,df=lampar$l,ncp=lampar$d,lower=FALSE)
  Q0 = (q0-lampar$muX)/lampar$sigmaX*lampar$sigmaQ+lampar$muQ
  Q1 = (Q0-lampar1$muQ)/lampar1$sigmaQ*lampar1$sigmaX+lampar1$muX
  pchisq(Q1,df=lampar1$l,ncp=lampar1$d,lower=FALSE)
}
#' Compute power for SKAT (given chi-square mixture parameters) using Davies method
#'
#' The Davies method is based on inverting characteristic function to compute the exact distribution function. 
#' @param  alpha significance level
#' @param  lambda  vector of mixing coefficients
#' @param  delta  vector of non-centrality parameters
#' @return estimated power
#' @export
#' @references
#' Wu,B. and Pankow,J.S. (2016) On sample size and power calculation for variant set-based association tests. /AHG/, 80(2), 136-143.
#'
#' Lafaye De Micheaux, P. (2013) CompQuadForm: Distribution function of quadratic forms in normal variables. R package version 1.4.1. 
davies.pwr <- function(alpha,lambda,delta){
  Q0 = KAT.qdf(alpha,lambda)
  KAT.tpr(Q0,lambda,delta, lim=1e8,acc=1e-12)
}


#' Compute parameters for (weighted linear) kernel association test of rare variant set
#'
#' Assume a linear model for continous traits, \eqn{Y=\alpha+G\beta+\epsilon}.
#' Input population genotype covariance matrix (VLD), variant weights (Wv), variant effect sizes (Ves, i.e., \eqn{\beta}), and sample size (n).
#' This function computes the asymptotic mixture of 1-DF chi-square distributions, \eqn{\sum_i\lambda_i\chi_1^2(\delta_i)}.
#' @param  n  sample size of a planned study
#' @param  VLD  population genotype covariance matrix for the variant set
#' @param  Ves  variant effect sizes
#' @param  MAF  minor allele freqs of variant set
#' @param  Wv  variant weights
#' @param  W.beta  beta dist parameters to compute variant weights based on MAF when Wv is not specified.
#' @return
#' \describe{
#'   \item{lambda}{ mixing cofficients of 1-DF chi-square mixtures }
#'   \item{delta}{ non-centrality parameters of 1-DF chi-square mixtures }
#' }
#' 
#' @export
#' @references
#' Wu,B. and Pankow,J.S. (2016) On sample size and power calculation for variant set-based association tests. /AHG/, 80(2), 136-143.
RVS.params <- function(n=2000, VLD, Ves, MAF, Wv=NULL, W.beta=c(1,25)){
  if(is.null(Wv)){
    W = dbeta(MAF,W.beta[1],W.beta[2])
    Wv = W/sum(W)*length(W)
  }
  lmf = log(MAF)
  pr1 =  1-exp(lmf*2*n)
  R = t(VLD*Wv)*Wv
  R1 = t(R*pr1)*pr1
  diag(R1) = diag(R1)/pr1
  a0 = svd(R1,nv=0)
  lam = a0$d
  dta = sqrt(lam)*colSums(a0$u*Ves/Wv)*sqrt(n)
  return(list(lambda=lam,delta=dta^2))
}

