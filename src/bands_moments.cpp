#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

// Finite-C raw CEM core adapted from CEClustNext, Copyright Thierry Dumont 2026, MIT.
// The unbounded branch explicitly rejects singular covariance states.
// [[Rcpp::export]]
Rcpp::List bands_state_cpp(const arma::mat& X,const Rcpp::IntegerVector& labels,
                         double C,double lambda) {
  using namespace Rcpp;
  int n=X.n_rows,d=X.n_cols,K=0;
  if(n<2||d<1||labels.size()!=n||!X.is_finite())stop("Invalid data or labels.");
  if(!(C>0)||R_IsNA(C)||R_IsNaN(C)||!R_finite(lambda)||lambda<0)
    stop("C must be positive (possibly Inf); lambda must be finite and nonnegative.");
  for(int label:labels)if(label==NA_INTEGER||label<1||label>n)stop("Invalid labels.");
  const bool unbounded=!R_finite(C);
  for(int i=0;i<n;++i)K=std::max(K,labels[i]);
  IntegerVector count(K),first(K,-1),cluster(K),ranks(K);
  arma::mat sum(K,d,arma::fill::zeros),means(K,d,arma::fill::zeros);
  std::vector<std::vector<int>> groups(K);
  for(int i=0;i<n;++i) {
    int k=labels[i]-1;
    if(k<0||k>=K)stop("Invalid labels.");
    if(first[k]<0)first[k]=i;
    ++count[k];groups[k].push_back(i);
    for(int j=0;j<d;++j)sum(k,j)+=X(i,j)-X(first[k],j);
  }
  for(int k=0;k<K;++k) {
    if(!count[k])stop("Labels must be compressed.");
    means.row(k)=X.row(first[k])+sum.row(k)/count[k];cluster[k]=k+1;
  }
  NumericVector mass(K),ldet(K),emin(K),emax(K);
  List covs(K);std::vector<int> bounded;
  arma::mat logdens(n,K);
  const double eps=std::numeric_limits<double>::epsilon();
  const double log2pi=std::log(2*M_PI);
  const double scale_min=1.0/(std::sqrt(2*M_PI)*std::pow(C,1.0/d));
  const double logdelta=unbounded ? R_NegInf : 2*d*std::log(scale_min);
  for(int k=0;k<K;++k) {
    arma::mat S(d,d,arma::fill::zeros);
    std::vector<double> offset(d),diff(d);
    for(int j=0;j<d;++j)offset[j]=sum(k,j)/count[k];
    for(int i:groups[k]) {
      for(int j=0;j<d;++j)diff[j]=(X(i,j)-X(first[k],j))-offset[j];
      for(int a=0;a<d;++a)for(int b=0;b<d;++b)S(a,b)+=diff[a]*diff[b];
    }
    S/=count[k];S=(S+S.t())/2;
    arma::vec ev;arma::mat vectors;
    if(!arma::eig_sym(ev,vectors,S))stop("Covariance eigendecomposition failed.");
    double tol=8*d*eps*arma::abs(ev).max();
    if(ev.min() < -tol)stop("Invalid covariance.");
    if(tol>0)for(int j=0;j<d;++j)
      if(ev[j]>.25*tol && ev[j]<4*tol)stop("Numerically ambiguous rank; use reference evaluation.");
    int rank=arma::accu(ev>tol);ranks[k]=rank;emin[k]=ev.min();emax[k]=ev.max();
    mass[k]=double(count[k])/n;
    if(rank<d) {
      if(unbounded)return List::create(_["valid"]=false,
        _["reason"]="unbounded_singular_covariance",_["group"]=k+1);
      ldet[k]=R_NegInf;bounded.push_back(k+1);
      arma::uvec ids=arma::find(ev<=tol);
      // Column order has no effect on distance to the null subspace.
      arma::mat basis=vectors.cols(ids);
      double support_tol=std::sqrt(tol);
      for(int i:groups[k])support_tol=std::max(support_tol,arma::norm((X.row(i)-means.row(k))*basis,2));
      double logC=-.5*d*log2pi-.5*logdelta;
      if(std::abs(logC)<=32*eps)logC=0;
      NumericMatrix stored=wrap(S);
      stored.attr("CEC_singular_support")=List::create(_["rank"]=rank,_["null_basis"]=basis,
        _["rank_tolerance"]=tol,_["support_tolerance"]=support_tol,_["log_C"]=logC);
      covs[k]=stored;
      double mean_scale=arma::abs(means.row(k)).max();
      for(int i=0;i<n;++i) {
        double residual=arma::norm((X.row(i)-means.row(k))*basis,2);
        double allowance=support_tol+64*d*eps*std::max(arma::abs(X.row(i)).max(),mean_scale);
        logdens(i,k)=residual<=allowance ? logC : R_NegInf;
      }
    } else {
      arma::mat L;
      if(!arma::chol(L,S,"lower"))stop("Covariance factorization failed; no jitter is added.");
      double ld=2*arma::sum(arma::log(L.diag()));ldet[k]=ld;
      double inflation=std::max(0.0,(logdelta-ld)/d);
      arma::mat fitted=S*std::exp(inflation);
      if(!fitted.is_finite()||!arma::chol(L,fitted,"lower"))stop("Fitted covariance factorization failed.");
      covs[k]=fitted;
      if(!unbounded && ld<=logdelta+1e-10*std::max(1.0,std::abs(logdelta)))bounded.push_back(k+1);
      double fitted_ld=2*arma::sum(arma::log(L.diag()));
      arma::mat centered=X.each_row()-means.row(k);
      arma::mat solved=arma::solve(arma::trimatl(L),centered.t());
      logdens.col(k)=-.5*(d*log2pi+fitted_ld)-.5*arma::sum(arma::square(solved),0).t();
    }
  }
  NumericMatrix mean_matrix=wrap(means);
  List params=List::create(_["states"]=cluster,_["nu"]=mass,_["m"]=mean_matrix,_["Sigma"]=covs,
    _["lambda"]=lambda,_["C"]=C,_["familyType"]="gaussVector",_["functionBoundReached"]=wrap(bounded));
  List stats=List::create(_["cluster"]=cluster,_["n"]=count,_["mass"]=mass,_["rank"]=ranks,
    _["logdet"]=ldet,_["eigmin"]=emin,_["eigmax"]=emax);
  stats.attr("class")="data.frame";
  stats.attr("row.names")=IntegerVector::create(NA_INTEGER,-K);
  return List::create(_["valid"]=true,_["stats"]=stats,
    _["params"]=params,_["logdens"]=logdens);
}

double bands_profile_cpp(const Rcpp::List& stats,double C,double lambda,int d,double beta) {
 Rcpp::NumericVector mass=stats["mass"],ld=stats["logdet"];
 Rcpp::IntegerVector rank=stats["rank"];
 double value=beta*(mass.size()-1),bound=-d*std::log(2*M_PI)-2*std::log(C);
 for(int k=0;k<mass.size();++k) {
  double ce=rank[k]<d ? -std::log(C) :
   .5*d*std::log(2*M_PI)+.5*std::max(ld[k],bound)+.5*d*std::exp(-std::max(0.,(bound-ld[k])/d));
  value+=mass[k]*(ce-lambda*std::log(mass[k]));
 }
 return value;
}

Rcpp::IntegerVector bands_compress_cpp(const Rcpp::IntegerVector& labels) {
 Rcpp::IntegerVector result(labels.size());std::vector<int> map(Rcpp::max(labels)+1,0);int K=0;
 for(int i=0;i<labels.size();++i) {int k=labels[i];if(!map[k])map[k]=++K;result[i]=map[k];}
 return result;
}

Rcpp::IntegerVector bands_propose_cpp(const Rcpp::List& state,double lambda,double floor) {
 Rcpp::NumericMatrix density=state["logdens"];Rcpp::List params=state["params"];
 Rcpp::NumericVector mass=params["nu"];Rcpp::IntegerVector labels(density.nrow());
 for(int i=0;i<density.nrow();++i) {
  int chosen=0;double best=std::max(floor,density(i,0))+lambda*std::log(mass[0]);
  for(int k=1;k<density.ncol();++k) {
   double value=std::max(floor,density(i,k))+lambda*std::log(mass[k]);
   if(value>best){best=value;chosen=k;}
  }
  if(!R_finite(best))Rcpp::stop("No supported CEM assignment");
  labels[i]=chosen+1;
 }
 return bands_compress_cpp(labels);
}

// [[Rcpp::export]]
Rcpp::List bands_cem_cpp(const arma::mat& X,Rcpp::IntegerVector labels,double C,double lambda,
                     int max_iter=1000,double tolerance=1e-8,double beta=0) {
 if(labels.size()!=X.n_rows||labels.size()==0)Rcpp::stop("Invalid label length.");
 for(int label:labels)if(label==NA_INTEGER||label<1||label>labels.size())Rcpp::stop("Invalid labels.");
 if(max_iter<1||!R_finite(tolerance)||tolerance<0||!R_finite(beta)||beta<0||
    !R_finite(lambda)||lambda<=0)Rcpp::stop("Invalid CEM controls.");
 labels=bands_compress_cpp(labels);Rcpp::List current=bands_state_cpp(X,labels,C,lambda);
 if(!Rcpp::as<bool>(current["valid"]))return Rcpp::List::create(
  Rcpp::_["labels"]=labels,Rcpp::_["J"]=NA_REAL,Rcpp::_["trace"]=Rcpp::NumericVector(),
  Rcpp::_["stop"]="unbounded_singular_initial",Rcpp::_["valid"]=false);
 double value=bands_profile_cpp(current["stats"],C,lambda,X.n_cols,beta);
 if(!R_finite(value))Rcpp::stop("Nonfinite initial CEM criterion.");
 std::vector<double> trace{value};std::string stop="iteration_budget";
 for(int iteration=0;iteration<max_iter;++iteration) {
  if(iteration%20==0)Rcpp::checkUserInterrupt();
  Rcpp::IntegerVector next=bands_propose_cpp(current,lambda,R_NegInf);
  if(Rcpp::is_true(Rcpp::all(next==labels))){stop="unchanged_partition";break;}
  Rcpp::List candidate=bands_state_cpp(X,next,C,lambda);
  if(!Rcpp::as<bool>(candidate["valid"])) {stop="unbounded_singular_proposal";break;}
  double after=bands_profile_cpp(candidate["stats"],C,lambda,X.n_cols,beta);
  if(!R_finite(after))Rcpp::stop("Nonfinite proposed CEM criterion.");
  if(after>=value-tolerance*std::max(1.,std::abs(value))){stop="no_improving_proposal";break;}
  current=candidate;labels=next;value=after;trace.push_back(value);
 }
 return Rcpp::List::create(Rcpp::_["labels"]=labels,Rcpp::_["valid"]=true,
  Rcpp::_["J"]=value,Rcpp::_["trace"]=trace,Rcpp::_["stop"]=stop,
  Rcpp::_["proposal"]="raw");
}

