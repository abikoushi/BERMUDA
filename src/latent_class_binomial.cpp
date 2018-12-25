#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

arma::rowvec colSums(const arma::mat & X){
  int nCols = X.n_cols;
  arma::rowvec out(nCols);
  for(int i = 0; i < nCols; i++){
    out[i] = sum(X.col(i));
  }
  return(out);
}

arma::vec colMeans(const arma::mat & X){
  int nCols = X.n_cols;
  arma::vec out(nCols);
  for(int i = 0; i < nCols; i++){
    out(i) = sum(X.col(i))/X.n_rows;
  }
  return(out);
}

arma::vec rowSums(const arma::mat & X){
  int nRows = X.n_rows;
  arma::vec out(nRows);
  for(int i = 0; i < nRows; i++){
    out(i) = sum(X.row(i));
  }
  return(out);
}

arma::rowvec softmax(const arma::rowvec & x){
  double tmp = max(x);
  arma::rowvec res = exp(x-tmp)/sum(exp(x-tmp));
  return res;
}

double logsumexp(const arma::rowvec & x){
  double tmp = max(x);
  return tmp + std::log(sum(exp(x-tmp)));
}


double multi_lpmf(const arma::rowvec & W, const double & M, const arma::rowvec & p){
  int K = W.n_cols;
  double lp = 0;
  for(int k=0; k<K; k++){
      if(!(W[k]==0 && p[k]==0)){
        lp = lp + W[k]*std::log(p[k])-std::lgamma(W[k]+1);
      }
    }
  lp = lp + std::lgamma(M+1);
  return lp;
}


// [[Rcpp::export]]
List doEM(arma::vec y, arma::mat W, int L, arma::vec phi, arma::mat p, arma::vec rho, double gamma, double beta){
  arma::vec M;
  M = rowSums(W);
  int K = W.n_cols;
  int N = y.n_rows;
  arma::mat z = arma::zeros(N,L);
  arma::rowvec unnorm(L);
  arma::rowvec unnorm2(L);
  arma::rowvec tmp;
  arma::vec tmp2;
  arma::vec phi2=phi;
  arma::mat p2=p;
  arma::vec rho2=rho;
  bool tol;
  bool nan;
  int iter = 0;
  double loglik;
  double loglik2;
  for(int i=0;i<10000;i++){
    iter = iter+1;
    loglik2 = 0;
  for(int n=0;n<N;n++){
    for(int l=0; l<L; l++){
      unnorm[l] = beta*(std::log(phi[l]) + y[n]*std::log(rho[l]) + (1-y[n])*log(1-rho[l]) + multi_lpmf(W.row(n),M[n],p.row(l)));
    }
    loglik2 += logsumexp(unnorm/beta);
    z.row(n) = softmax(unnorm);
  }
  phi2  = colMeans(z);
  for(int l=0; l<L; l++){
    rho2[l] = sum(z.col(l) % y) / sum(z.col(l));
    p2.row(l) = (colSums(W.each_col() % z.col(l)) + gamma)/(sum(M % z.col(l))+K*gamma);
  }
  tol = all(abs(phi2-phi)<1e-8) && all(all(abs(p2-p)<1e-8)) && all(abs(rho2-rho)<1e-8) && std::abs(loglik2-loglik)<1e-8;
  nan = phi2.has_nan()|p2.has_nan()|rho2.has_nan();
  if(tol|nan){
    break;
  }
  p = p2;
  phi = phi2;
  rho = rho2;
  loglik = loglik2;
  }
  return List::create(Named("p")=p2,_["phi"]=phi2,_["rho"]=rho2,_["z"]=z,_["loglik"]=loglik2,_["iter"]=iter);
}

// [[Rcpp::export]]
arma::vec doPred(arma::mat W, int L, arma::vec phi, arma::mat p, arma::vec rho){
  arma::vec M;
  M = rowSums(W);
  int N = W.n_rows;
  arma::vec y(N);
  arma::rowvec unnorm(L);
  arma::rowvec unnorm2(L);
  arma::rowvec tmp;
  arma::vec tmp2;
  arma::vec phi2=phi;
  arma::mat p2=p;
  arma::vec rho2=rho;
  for(int n=0;n<N;n++){
    for(int l=0; l<L; l++){
      unnorm[l] = std::log(phi[l]) + log(rho[l]) + multi_lpmf(W.row(n),M[n],p.row(l));
      unnorm2[l] = std::log(phi[l]) + multi_lpmf(W.row(n),M[n],p.row(l));
    }
    y[n] = std::exp(logsumexp(unnorm)-logsumexp(unnorm2));
  }
  return y;
}
