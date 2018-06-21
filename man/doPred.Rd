\name{doPred}
\alias{doPred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict binary response using multinomial count.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
doPred(W, L, phi, p, rho)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{
A matrix of multinial count.
}
  \item{L}{
Number of clusters.
}
  \item{phi}{
Mixing ratio parameters of each components.
}
  \item{p}{
The parameteres of multinomial distribution.
}
  \item{rho}{
The parameteres of Bernoulli distribution.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%
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
