\name{doEM}
\alias{doEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter estimation using AEM algorithm
}
\description{
Return EM algorithm output for mixtures of multinomial distributions.
}
\usage{
doEM(y, W, L, phi, p, rho, gamma, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Binary outcome of Bernoulli response.
}
  \item{W}{
A matrix of multinial count.
}
  \item{L}{
Number of clusters.
}
  \item{phi}{
Initial value of mixing ratio of each components.
}
  \item{p}{
Initial value of multinial paramaters.
}
  \item{rho}{
Initial value of Bernoulli paramaters.
}
  \item{gamma}{
Hyper parameter of Dirichlet prior.
}
  \item{beta}{
Inverse temparture for anealing EM algorithm.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{p}{The final multinomial parameters.}
\item{phi}{The final mixing proportions.}
\item{rho}{The final Bernoulli parameters.}
\item{z}{The probability which each samples belong to the cluster.}
\item{loglik}{The final value of log-likelihood.}
\item{iter}{Iteration count until convergence.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
