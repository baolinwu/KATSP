# KATSP
Sample size and power calculation for variant-set based association tests

# Reference
Wu,B. and Pankow,J.S. (2016) On sample size and power calculation for variant set-based association tests. *AHG*, 80(2), 136-143.

# Sample R codes

```r
library(KATSP)
library(minqa)
library(CompQuadForm)
## simulate genotype and es for QT
## simulate 1e4 by 20 genotypes from MVN with pairwise corr=0.1
Z1 = matrix(rnorm(1e4*20,0,sqrt(0.9)),1e4,20) + rnorm(1e4,0,sqrt(0.1))
Z2 = matrix(rnorm(1e4*20,0,sqrt(0.9)),1e4,20) + rnorm(1e4,0,sqrt(0.1))
## population MAF U[0.0005,0.02]
maf = runif(20, 0.0005,0.02)
q0 = qnorm(maf, lower=FALSE)
G = t( I(t(Z1)>q0) + I(t(Z2)>q0) )
## simulate es
VLD = cov(G)
MAF = colMeans(G)/2
Ves = runif(20,-0.25,0.25)
a1 = RVS.params(n=5e3,VLD,Ves,MAF)
## power
alpha = 2.5e-6
davies.pwr(alpha,a1$lambda,a1$delta)
wu.pwr(alpha,a1$lambda,a1$delta)
lee.pwr(alpha,a1$lambda,a1$delta)
liu.pwr(alpha,a1$lambda,a1$delta)
## given lambda/delta
lam = c(0.012, 0.0076, 0.0061, 0.0035, 0.0018, 0.0013, 0.00092, 0.00055, 0.00037)
dta = c(0.028, -0.044, -0.033, 0.032, -0.015, 0.0073, 0.0036, -0.012, 0.011)
davies.pwr(alpha,lam,dta^2*5000)
wu.pwr(alpha,lam,dta^2*5000)
lee.pwr(alpha,lam,dta^2*5000)
liu.pwr(alpha,lam,dta^2*5000)
```
