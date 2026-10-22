#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

struct EMState { mat means; cube cov; vec weights; mat prob; double ll; double det_margin; };

EMState em_mstep(const mat& X, const mat& R, double logdelta) {
  const int n=X.n_rows,d=X.n_cols,k=R.n_cols;
  EMState s; s.means.set_size(k,d); s.cov.set_size(d,d,k);s.det_margin=arma::datum::inf;
  s.weights=sum(R,0).t()/n;
  for(int j=0;j<k;j++) {
    double mass=accu(R.col(j));
    if(!(mass>0)) stop("Zero numerical component mass");
    mat shifted=X.each_row()-X.row(0);
    rowvec offset=(R.col(j).t()*shifted)/mass;
    rowvec mu=X.row(0)+offset; s.means.row(j)=mu;
    mat D=shifted.each_row()-offset; mat WD=D.each_col()%R.col(j);
    mat V=symmatu(D.t()*WD/mass); mat L;
    if(!chol(L,V,"lower")) stop("Singular weighted covariance: candidate rejected, no eigenvalue floor");
    double ld=2*accu(log(L.diag()));
    V*=std::exp(std::max(0.0,(logdelta-ld)/d));
    if(!chol(L,V,"lower"))stop("Numerical singularity after covariance scaling");
    double margin=2*accu(log(L.diag()))-logdelta;
    if(margin < -1e-8)stop("Numerical determinant constraint violation");
    s.det_margin=std::min(s.det_margin,margin);s.cov.slice(j)=V;
  }
  return s;
}

mat em_logjoint(const mat& X, const mat& means, const cube& cov, const vec& weights) {
  int n=X.n_rows,d=X.n_cols,k=weights.n_elem; mat scores(n,k);
  for(int j=0;j<k;j++) {
    mat L; if(!chol(L,cov.slice(j),"lower"))stop("Non-positive EM covariance");
    mat D=X.each_row()-means.row(j); mat z=solve(trimatl(L),D.t());
    scores.col(j)=std::log(weights(j))-.5*d*std::log(2*M_PI)-accu(log(L.diag()))-.5*sum(square(z),0).t();
  }
  return scores;
}
void em_estep(const mat& X, EMState& s) {
  mat scores=em_logjoint(X,s.means,s.cov,s.weights); vec mx=max(scores,1);
  mat e=exp(scores.each_col()-mx); vec den=sum(e,1);
  s.ll=mean(mx+log(den)); s.prob=e.each_col()/den;
  if(!std::isfinite(s.ll))stop("Nonfinite EM likelihood");
}
// [[Rcpp::export]]
List emmap_fit_cpp(const arma::mat& X, const IntegerVector& labels, double C,
                   int max_iter=1000, double tolerance=1e-8, double soften=.02) {
  int n=X.n_rows,d=X.n_cols;
  if(n<1 || d<1 || !X.is_finite() || labels.size()!=n ||
     !std::isfinite(C) || C<=0 || max_iter<1 || !std::isfinite(tolerance) || tolerance<0 ||
     !std::isfinite(soften) || soften<=0 || soften>=1) stop("Invalid EM input");
  for(int i=0;i<n;i++) if(IntegerVector::is_na(labels[i]) || labels[i]<1 || labels[i]>n) stop("Invalid EM labels");
  int k=max(labels); mat R(n,k,fill::value(soften/k));
  for(int i=0;i<n;i++)R(i,labels[i]-1)+=1-soften;
  double logdelta=-d*std::log(2*M_PI)-2*std::log(C);
  EMState s=em_mstep(X,R,logdelta);em_estep(X,s);
  std::vector<double> trace;trace.push_back(s.ll);std::string reason="EM_iteration_budget";
  int quiet=0;double min_margin=s.det_margin;
  for(int it=0;it<max_iter;it++) {
    EMState next;
    try { next=em_mstep(X,s.prob,logdelta);em_estep(X,next); }
    catch(std::exception &ex) { reason="EM_singular_or_numerical_limit";break; }
    double gain=next.ll-s.ll,scale=std::max(1.0,std::abs(s.ll));
    if(gain < -1e-9*scale){reason="EM_numerical_decrease_rejected";break;}
    if(gain<0){reason="EM_numerical_plateau";break;}
    s=next; trace.push_back(s.ll);min_margin=std::min(min_margin,s.det_margin);
    quiet=gain<=tolerance*scale ? quiet+1 : 0;
    if(quiet>=3){reason="EM_likelihood_stable";break;}
    if(it%50==0)checkUserInterrupt();
  }
  IntegerVector lab(n);for(int i=0;i<n;i++)lab[i]=s.prob.row(i).index_max()+1;
  return List::create(_["labels"]=lab,_["means"]=s.means,_["cov"]=s.cov,
    _["weights"]=s.weights,_["loglik_mean_trace"]=trace,_["stop"]=reason,
    _["converged"]=reason=="EM_likelihood_stable" || reason=="EM_numerical_plateau",_["soft_K"]=k,
    _["minimum_logdet_constraint_margin"]=min_margin);
}
// [[Rcpp::export]]
IntegerVector emmap_predict_cpp(const arma::mat& X, const arma::mat& means,
                                const arma::cube& cov, const arma::vec& weights) {
  if(!X.is_finite() || means.n_cols!=X.n_cols || means.n_rows!=weights.n_elem ||
     cov.n_slices!=weights.n_elem || cov.n_rows!=X.n_cols || cov.n_cols!=X.n_cols ||
     weights.n_elem<1 || !weights.is_finite() || any(weights<=0)) stop("Invalid EM prediction input");
  mat scores=em_logjoint(X,means,cov,weights);IntegerVector labels(X.n_rows);
  for(unsigned int i=0;i<X.n_rows;i++)labels[i]=scores.row(i).index_max()+1;
  return labels;
}
