#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include <chrono>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

// Shared sufficient statistics, determinant-only covariance fitting and raw
// reference log densities. No change to initialization or search stopping.
// [[Rcpp::export]]
Rcpp::List cec_cpp_cem_evaluate(const arma::mat& X,const Rcpp::IntegerVector& labels,
                         double C,double lambda) {
  using namespace Rcpp;
  int n=X.n_rows,d=X.n_cols,K=0;
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
  const double logdelta=2*d*std::log(scale_min);
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
      if(ld<=logdelta+1e-10*std::max(1.0,std::abs(logdelta)))bounded.push_back(k+1);
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
  return List::create(_["stats"]=stats,_["params"]=params,_["logdens"]=logdens);
}

// [[Rcpp::export]]
Rcpp::List cec_cpp_cem_descent(const arma::mat& X,double C,double lambda,
                           Rcpp::IntegerVector labels,int max_iter,double tolerance,
                           double beta,Rcpp::Function profile) {
  using namespace Rcpp;
  const int n=X.n_rows,d=X.n_cols;
  const double logfloor=std::min(-std::log(C),std::log(C)-std::log(1000.0));
  auto clock=[]() {return std::chrono::steady_clock::now();};
  auto seconds=[](std::chrono::steady_clock::time_point a) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now()-a).count();
  };
  double refit_time=0,proposal_time=0,deletion_time=0;
  std::vector<int> stepK;
  std::vector<double> raw_scores,penalties,scores,full_penalties,loss_trace;
  std::vector<std::string> moves;
  std::vector<List> deletion_plans;
  int deletion_count=0;
  auto better=[&](double a,double b) {return R_finite(a)&&(!R_finite(b)||a<b-tolerance*std::max(1.0,std::abs(b)));};
  auto compress=[&](const IntegerVector& a,int oldK) {
    IntegerVector out(n);std::vector<int> mapping(oldK+1,0);int next=0;
    for(int i=0;i<n;++i) {int k=a[i];if(!mapping[k])mapping[k]=++next;out[i]=mapping[k];}
    return out;
  };
  List current;
  List part;
  arma::mat density;
  arma::vec mass,log_weights;
  auto evaluate=[&](const IntegerVector& lab) {
    auto begin=clock();
    List e=cec_cpp_cem_evaluate(X,lab,C,lambda);
    List stats=e["stats"];
    List p=profile(stats,C,lambda,d);
    p.push_front(stats,"stats");p.push_front(lab,"labels");
    e["part"]=p;
    arma::mat raw=as<arma::mat>(e["logdens"]);
    long double own=0;
    for(int i=0;i<n;++i)own+=raw(i,lab[i]-1);
    double from_density=lambda*as<double>(p["h"])-double(own/n);
    if(std::abs(from_density-as<double>(p["J_raw"]))>=1e-6)
      stop("Fused density/profile consistency check failed.");
    refit_time+=seconds(begin);return e;
  };
  auto activate=[&](List e) {
    current=e;part=as<List>(e["part"]);density=as<arma::mat>(e["logdens"]);
    List params=e["params"];mass=as<arma::vec>(params["nu"]);
    log_weights=lambda*arma::log(mass);
  };
  auto append=[&](const std::string& move,double loss=NA_REAL) {
    stepK.push_back(as<int>(part["K"]));raw_scores.push_back(as<double>(part["J_raw"]));
    penalties.push_back(as<double>(part["BIC_penalty_relative"]));scores.push_back(as<double>(part["J"]));
    full_penalties.push_back(as<double>(part["BIC_penalty"]));loss_trace.push_back(loss);moves.push_back(move);
  };
  activate(evaluate(labels));double initial=as<double>(part["J"]);append("initial");
  std::string reason;
  for(;;) {
    std::string cem_stop="iteration_limit";
    for(int iter=0;iter<max_iter;++iter) {
      if((iter&15)==0)Rcpp::checkUserInterrupt();
      auto begin=clock();
      int K=mass.n_elem;IntegerVector next(n);
      std::vector<double> best(n);
      const double* raw=density.memptr();
      for(int i=0;i<n;++i) {best[i]=std::max(logfloor,raw[i])+log_weights[0];next[i]=1;}
      for(int k=1;k<K;++k) {
        const double* column=raw+static_cast<size_t>(k)*n;double weight=log_weights[k];
        for(int i=0;i<n;++i) {
          double candidate=std::max(logfloor,column[i])+weight;
          if(candidate>best[i]) {best[i]=candidate;next[i]=k+1;}
        }
      }
      next=compress(next,K);
      bool unchanged=true;for(int i=0;i<n;++i)if(next[i]!=labels[i]) {unchanged=false;break;}
      proposal_time+=seconds(begin);
      if(unchanged) {cem_stop="unchanged_partition";break;}
      List candidate=evaluate(next);List p=candidate["part"];
      if(as<int>(p["K"])>K)stop("CEM proposal increased order.");
      if(!better(as<double>(p["J"]),as<double>(part["J"]))) {cem_stop="no_penalized_CEM_improvement";break;}
      labels=next;activate(candidate);append("CEM");
    }
    if(cem_stop=="iteration_limit") {reason="CEM_iteration_budget";break;}
    auto begin=clock();
    int K=mass.n_elem;NumericVector costs(K);IntegerVector alternative(n);
    int chosen=-1;double loss=R_PosInf;
    if(K>1) {
      for(int i=0;i<n;++i) {
        int own=labels[i]-1,best=-1;double value=R_NegInf;
        for(int k=0;k<K;++k)if(k!=own) {
          double score=density(i,k)+log_weights[k];
          if(best<0||score>value) {best=k;value=score;}
        }
        alternative[i]=best+1;
        costs[own]+=(density(i,own)+log_weights[own]-value)/n;
      }
      int candidate=0;for(int k=1;k<K;++k)if(costs[k]<costs[candidate])candidate=k;
      double threshold=beta-tolerance*std::max(1.0,std::abs(as<double>(part["J"])));
      if(R_finite(costs[candidate]) && costs[candidate]<threshold) {chosen=candidate;loss=costs[candidate];}
    } else costs=NumericVector();
    deletion_plans.push_back(List::create(_["K"]=K,_["costs"]=costs,
      _["selected"]=chosen<0?NA_INTEGER:chosen+1,_["loss"]=loss,_["beta"]=beta,_["CEM_stop"]=cem_stop));
    deletion_time+=seconds(begin);
    if(chosen<0) {reason="no_improving_frozen_deletion";break;}
    IntegerVector next=clone(labels);
    for(int i=0;i<n;++i)if(labels[i]==chosen+1)next[i]=alternative[i];
    next=compress(next,K);List e=evaluate(next);List p=e["part"];
    if(as<int>(p["K"])!=K-1 || !better(as<double>(p["J"]),as<double>(part["J"])) ||
       as<double>(p["J"])-as<double>(part["J"])>loss-beta+1e-6)stop("Deletion acceptance check failed.");
    labels=next;activate(e);append("delete_component",loss);++deletion_count;
  }
  int steps=scores.size();IntegerVector index(steps);NumericVector full_score(steps);
  double max_increase=0;
  for(int i=0;i<steps;++i) {
    index[i]=i;full_score[i]=raw_scores[i]+full_penalties[i];
    if(i) {
      double change=scores[i]-scores[i-1];if(change>=0)stop("Nondecreasing accepted trace.");
      if(i==1||change>max_increase)max_increase=change;
    }
  }
  DataFrame trace=DataFrame::create(_["step"]=index,_["move"]=moves,_["K"]=stepK,_["CEC"]=raw_scores,
    _["penalty_relative"]=penalties,_["J"]=scores,_["full_BIC_penalty"]=full_penalties,
    _["full_J"]=full_score,_["frozen_deletion_loss"]=loss_trace);
  return List::create(_["part"]=part,_["valid"]=true,_["error"]="",_["initial_J"]=initial,
    _["steps"]=steps-1,_["stop"]=reason,_["descent_trace"]=trace,_["deletion_plans"]=deletion_plans,
    _["deletion_count"]=deletion_count,_["max_increase"]=max_increase,
    _["stage_timings"]=NumericVector::create(_["refit"]=refit_time,_["CEM_proposal"]=proposal_time,_["deletion"]=deletion_time));
}
