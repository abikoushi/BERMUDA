# BERMUDA
BERMUDA is an R package for BERnoulli and MUltinomial Distribution-base latent Allocation

### Depends:

R(>=3.5.0)

Rcpp, RcppArmadillo

### Authors:

Ko Abe and Teppei Shimamura

Contact: ko.abe[at]med.nagoya-u.ac.jp and shimamura[at]med.nagoya-u.ac.jp

## Installation

Install the latest version of this package from Github by pasting in the following.

~~~R
devtools::install_github("abikoushi/BERMUDA")
~~~

## An example of synthetic data

~~~R
library(BERMUDA)
N <-100
L <-3
K <- 20
M <- 1000
rho <- c(3,5,7)/10
phi <- rep(1/L,L)
set.seed(1192)
unnorm <- matrix(rexp(L*5),L,5)
p <- unnorm/rowSums(unnorm)
z <- sample.int(L,N,replace = TRUE,prob = phi)
W <-t(sapply(1:N, function(i){rmultinom(1,M,p[z[i],])}))
y <-rbinom(N,1,rho[z])
fit <-doEM(y,W,L,
                phi = phi,rho = rho,
                p = p,gamma=0,beta = 1)
print(fit$rho)
print(rho)

doPred(W, L, phi = fit$phi, p = fit$p, rho = fit$rho)
~~~

## Genral overview
We propose new probabilistic model called BERMUDA (BERnoulli and MUltinomial Distribution-base latent Allocation).

BERMUDA enables us to describe the differences in bacteria composition and disease among samples.

In BERMUDA, abundances of each taxon can be viewed as a mixture of vari- ous groups where enables us to describe the differences in bacteria composition among samples.

For details, see Abe et al.(2018).

## Reference
Ko Abe, Masaaki Hirayama, Kinji Ohno, and Teppei Shimamura, A Latent Allocation Model for the Analysis of Microbial
Composition and Disease, submitted.
