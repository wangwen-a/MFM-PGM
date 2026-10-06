// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#ifndef ORDINAL_MFM_KERNELS_H
#define ORDINAL_MFM_KERNELS_H
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
namespace ordinal_mfm {
  using arma::mat;
  using arma::vec;
  const double LOG2PI=std::log(2.0*3.14159265358979323846);
  inline double logu(){
    return std::log(R::runif(0.0,1.0));
  }
  inline void check_graph(const mat&G) {
    if(G.n_rows!=G.n_cols || !G.is_finite()) Rcpp::stop("Invalid graph.");
    for(unsigned i=0;i<G.n_rows;++i) for(unsigned j=0;j<G.n_cols;++j)
    if((G(i,j)!=0 && G(i,j)!=1) || G(i,j)!=G(j,i) || (i==j && G(i,j)!=0))
    Rcpp::stop("Graph must be symmetric, binary, with zero diagonal.");
  }
  // Complete non-free Cholesky entries so that missing edges imply zero precision entries.
  inline mat complete(const mat&T,const mat&G) {
    const int p=G.n_rows;
    mat F=arma::trimatu(T);
    for(int i=0;i<p;++i) {
      if(!(F(i,i)>0) || !std::isfinite(F(i,i))) Rcpp::stop("Nonpositive Cholesky diagonal.");
      for(int j=i+1;j<p;++j) if(G(i,j)==0) {
        double s=0;
        for(int k=0;k<i;++k)s+=F(k,i)*F(k,j);
        F(i,j)=-s/F(i,i);
      }
    }
    return F;
  }
  // Unnormalized G-Wishart log density in free Cholesky coordinates, including its Jacobian.
  // The graph-independent factor 2^p cancels from all ratios.
  inline double log_h(const mat&F,const mat&G,double b,const mat&D) {
    if(!F.is_finite())return -INFINITY;
    double ans=0;
    for(unsigned i=0;i<G.n_rows;++i){
      int nu=0;
      for(unsigned j=i+1;j<G.n_rows;++j)nu+=(G(i,j)!=0);
      if(!(F(i,i)>0))return -INFINITY;
      ans+=(b+nu-1.0)*std::log(F(i,i));
    }
    mat O=F.t()*F;
    ans-=0.5*arma::accu(D%O);
    return std::isfinite(ans)?ans:-INFINITY;
  }
  struct WishartDraw {
    mat omega;
    mat phi;
    int attempts;
  };
  // Envelope rejection draw. The released sampler calls this only with b=3 and D=I.
  // Exhausting the attempt budget aborts; it never returns a substitute distribution.
  inline WishartDraw gwish(const mat&G,double b,const mat&D,int max_attempts=1000000) {
    check_graph(G);
    const int p=G.n_rows;
    if(!(b>2) || D.n_rows!=G.n_rows || D.n_cols!=G.n_cols || !D.is_finite())
    Rcpp::stop("G-Wishart rejection requires b>2 and finite compatible SPD D.");
    mat invD,R;
    if(!arma::inv_sympd(invD,D)||!arma::chol(R,invD))Rcpp::stop("D must be SPD.");
    for(int attempt=1;attempt<=max_attempts;++attempt) {
      if(attempt%4096==0)Rcpp::checkUserInterrupt();
      mat P(p,p,arma::fill::zeros),F(p,p,arma::fill::zeros);
      double nonfree_sq=0;
      for(int i=0;i<p;++i) {
        int nu=0;
        for(int j=i+1;j<p;++j)nu+=(G(i,j)!=0);
        P(i,i)=std::sqrt(R::rchisq(b+nu));
        F(i,i)=P(i,i)*R(i,i);
        for(int j=i+1;j<p;++j) {
          double previous=0;
          for(int k=i;k<j;++k)previous+=P(i,k)*R(k,j);
          if(G(i,j)!=0) {
            P(i,j)=R::rnorm(0,1);
            F(i,j)=previous+P(i,j)*R(j,j);
          } else {
            double s=0;
            for(int k=0;k<i;++k)s+=F(k,i)*F(k,j);
            F(i,j)=-s/F(i,i);
            P(i,j)=(F(i,j)-previous)/R(j,j);
            nonfree_sq+=P(i,j)*P(i,j);
          }
        }
      }
      if(std::isfinite(nonfree_sq)&&logu() < -0.5*nonfree_sq)
      return {
        F.t()*F,F,attempt
      };
    }
    Rcpp::stop("Exact G-Wishart rejection budget exhausted; run aborted, no fallback sample.");
    return {
    };
  }
  // Stable normal rejection: exponential tail envelope, uniform central envelope, or normal proposals.
  // No clipping of CDF probabilities and no deterministic fallback in a stochastic update.
  inline double std_tn(double a,double b) {
    if(std::isnan(a)||std::isnan(b)||!(a<b))Rcpp::stop("Empty truncated-normal interval.");
    if(std::isfinite(a)&&std::isfinite(b)&&std::nextafter(a,b)>=b)
    Rcpp::stop("Truncated-normal interval has no representable interior double.");
    if(b<=0)return -std_tn(-b,-a);
    for(int n=0;n<1000000;++n) {
      double z;
      if(a>=0) {
        double alpha=0.5*a+0.5*std::hypot(a,2.0);
        double mass=std::isfinite(b)?-std::expm1(-alpha*(b-a)):1.0;
        z=a-std::log1p(-R::runif(0,1)*mass)/alpha;
        if(z>a&&z<b&&logu() < -0.5*(z-alpha)*(z-alpha))return z;
      } else if(std::isfinite(a)&&std::isfinite(b)&&b-a<2.0) {
        z=a+(b-a)*R::runif(0,1);
        if(z>a&&z<b&&logu() < -0.5*z*z)return z;
      } else {
        z=R::rnorm(0,1);
        if(z>a&&z<b)return z;
      }
    }
    Rcpp::stop("Truncated-normal rejection budget exhausted; no deterministic fallback.");
    return NA_REAL;
  }
  inline double tn(double mean,double sd,double lo,double hi) {
    if(!std::isfinite(mean)||!std::isfinite(sd)||!(sd>0)||!(lo<hi))
    Rcpp::stop("Invalid truncated-normal parameters.");
    const double a=(lo-mean)/sd,b=(hi-mean)/sd;
    for(int k=0;k<1000;++k) {
      double z=mean+sd*std_tn(a,b);
      if(std::isfinite(z)&&z>lo&&z<hi)return z;
    }
    Rcpp::stop("Truncated-normal rounding cannot represent a valid draw.");
    return NA_REAL;
  }
  inline mat augment(const mat&F,const mat&G,const vec&scale) {
    mat T=F;
    for(unsigned i=0;i<G.n_rows;++i)for(unsigned j=i+1;j<G.n_rows;++j)
    if(G(i,j)==0)T(i,j)=F(i,j)+R::rnorm(0,scale[j]);
    return T;
  }
  inline double pseudo_log(const mat&T,const mat&F,const mat&G,const vec&scale) {
    double x=0;
    for(unsigned i=0;i<G.n_rows;++i)for(unsigned j=i+1;j<G.n_rows;++j)
    if(G(i,j)==0){
      double v=(T(i,j)-F(i,j))/scale[j];
      x+=-0.5*v*v-std::log(scale[j])-0.5*LOG2PI;
    }
    return x;
  }
  struct GraphDraw {
    mat graph;
    mat omega;
    int graph_accept;
    int precision_accept;
    int auxiliary_attempts;
  };
  // Validate G and Omega BEFORE completion or RNG use. A mismatch is an error, not a projection.
  // Standardized nonedge tolerance permits floating-point roundoff only.
  inline mat checked_precision_cholesky(const mat&G,const mat&O,double tol=1e-10) {
    check_graph(G);
    const unsigned p=G.n_rows;
    if(p==0 || O.n_rows!=p || O.n_cols!=p || !O.is_finite())
    Rcpp::stop("Input precision must be finite and match graph dimensions.");
    if(!(tol>0) || !std::isfinite(tol))Rcpp::stop("Invalid graph-precision tolerance.");
    for(unsigned i=0;i<p;++i)
    if(!(O(i,i)>0))Rcpp::stop("Input precision must have positive diagonal entries.");
    for(unsigned i=0;i<p;++i)for(unsigned j=i+1;j<p;++j) {
      const double upper=O(i,j)/std::sqrt(O(i,i))/std::sqrt(O(j,j));
      const double lower=O(j,i)/std::sqrt(O(i,i))/std::sqrt(O(j,j));
      if(!std::isfinite(upper)||!std::isfinite(lower)||std::abs(upper-lower)>tol)
      Rcpp::stop("Input precision must be symmetric before graph update.");
      if(G(i,j)==0 && std::max(std::abs(upper),std::abs(lower))>tol)
      Rcpp::stop("Graph-precision mismatch at (%d,%d): standardized nonedge magnitude %.6g exceeds %.6g; no projection performed.",
      (int)i+1,(int)j+1,std::max(std::abs(upper),std::abs(lower)),tol);
    }
    mat F;
    if(!arma::chol(F,O))Rcpp::stop("Input precision must be SPD.");
    return F;
  }
  // Discrete-time graph/precision transition. Return the terminal pair, not an inner average.
  // Each step updates one free precision coordinate and proposes one edge toggle.
  // Auxiliary prior draws and normalized pseudo-priors cancel unknown graph normalizers.
  inline GraphDraw graph_step(mat G,mat O,double b,double bs,const mat&D,const mat&Ds,
  double q,int steps) {
    check_graph(G);
    if(!(q>0&&q<1)||steps<1)Rcpp::stop("Invalid graph-kernel controls.");
    const int p=G.n_rows;
    // Validate the exported entry before using scales or drawing any random numbers.
    // This sampler implements the b>2 G-Wishart parameter domain.
    if(!std::isfinite(b)||!std::isfinite(bs)||b<=2.0||bs<=2.0)
      Rcpp::stop("G-Wishart b and b_star must be finite and greater than 2.");
    auto check_scale=[p](const mat& A,const char* name){
      if(A.n_rows!=(unsigned)p||A.n_cols!=(unsigned)p||!A.is_finite())
        Rcpp::stop("%s must be a finite p-by-p matrix.",name);
      const double tol=1e-12*std::max(1.0,arma::abs(A).max());
      if(arma::abs(A-A.t()).max()>tol)
        Rcpp::stop("%s must be symmetric.",name);
      mat C;
      if(!arma::chol(C,A))Rcpp::stop("%s must be positive definite.",name);
    };
    check_scale(D,"D");check_scale(Ds,"D_star");
    mat F=checked_precision_cholesky(G,O);
    F=complete(F,G);
    vec scale=1.5/arma::sqrt(Ds.diag());
    int ga=0,pa=0,attempts=0;
    for(int it=0;it<steps;++it) {
      if(it%256==0)Rcpp::checkUserInterrupt();
      std::vector<std::pair<int,int>> free;
      for(int i=0;i<p;++i){
        free.push_back({i,i});
        for(int j=i+1;j<p;++j)if(G(i,j))free.push_back({i,j});
      }
      auto ij=free[std::min((int)free.size()-1,(int)(R::runif(0,1)*free.size()))];
      int i=ij.first,j=ij.second;
      mat candidate=F;
      double jac=0;
      // A log-diagonal random walk needs its proposal Jacobian; off-diagonal moves are symmetric.
      if(i==j){
        double change=R::rnorm(0,1.0/std::sqrt(bs));
        candidate(i,i)*=std::exp(change);
        jac=change;
      }
      else candidate(i,j)+=R::rnorm(0,2.0/std::sqrt(Ds(j,j)));
      candidate=complete(candidate,G);
      if(logu()<log_h(candidate,G,bs,Ds)-log_h(F,G,bs,Ds)+jac){
        F=candidate;
        ++pa;
      }
      if(p>1) {
        int e=(int)(R::runif(0,1)*(p*(p-1)/2)),ix=0;
        i=0;
        j=1;
        for(int ii=0;ii<p-1;++ii)for(int jj=ii+1;jj<p;++jj){
          if(ix==e){
            i=ii;
            j=jj;
          }++ix;
        }
        mat Gp=G;
        Gp(i,j)=Gp(j,i)=1-G(i,j);
        mat T=augment(F,G,scale),Fp=complete(T,Gp);
        WishartDraw aux=gwish(Gp,b,D);
        attempts+=aux.attempts;
        mat U=augment(aux.phi,Gp,scale),Uold=complete(U,G);
        const double prior_graph=(Gp(i,j)-G(i,j))*std::log(q/(1-q));
        // Exchange acceptance ratio on the common augmented coordinate space.
        double ratio=prior_graph+
        log_h(Fp,Gp,bs,Ds)+pseudo_log(T,Fp,Gp,scale)-log_h(F,G,bs,Ds)-pseudo_log(T,F,G,scale)
        +log_h(Uold,G,b,D)+pseudo_log(U,Uold,G,scale)-log_h(aux.phi,Gp,b,D)-pseudo_log(U,aux.phi,Gp,scale);
        if(std::isnan(ratio))Rcpp::stop("Non-finite exchange acceptance ratio.");
        if(logu()<ratio){
          G=Gp;
          F=Fp;
          ++ga;
        }
      }
    }
    return {
      G,F.t()*F,ga,pa,attempts
    };
  }
  struct Component{
    mat graph;
    mat omega;
    vec mean;
  };
  // Draw a complete component from p(G) p(Omega|G) p(mu|Omega); no data enter here.
  inline Component prior_component(int p,double q,double a) {
    mat G(p,p,arma::fill::zeros),O(p,p,arma::fill::eye),F;
    {
      for(int i=0;i<p;++i)for(int j=i+1;j<p;++j)G(i,j)=G(j,i)=R::rbinom(1,q);
      WishartDraw x=gwish(G,3.0,arma::eye(p,p));
      O=x.omega;
      F=x.phi;
    }
    vec z(p);
    for(int j=0;j<p;++j)z[j]=R::rnorm(0,1);
    vec mu=arma::solve(arma::trimatu(F),z)/std::sqrt(a);
    return {
      G,O,mu
    };
  }
  inline double loglike(const vec&z,const Component&x) {
    mat R;
    if(!arma::chol(R,x.omega))Rcpp::stop("Component precision is not SPD.");
    vec d=z-x.mean;
    return arma::sum(arma::log(R.diag()))-0.5*z.n_elem*LOG2PI-0.5*arma::dot(d,x.omega*d);
  }
}
#endif

using namespace Rcpp;
// [[Rcpp::export]]
double draw_truncated_normal(double mean,double sd,double lower,double upper) {
  return ordinal_mfm::tn(mean,sd,lower,upper);
}
// [[Rcpp::export]]
List draw_gwishart(int n,NumericMatrix graph,double b,NumericMatrix D,int max_attempts=1000000) {
  List values(n);
  IntegerVector attempts(n);
  for(int i=0;i<n;++i){
    auto x=ordinal_mfm::gwish(as<arma::mat>(graph),b,as<arma::mat>(D),max_attempts);
    values[i]=x.omega;
    attempts[i]=x.attempts;
  }
  return List::create(Named("draws")=values,Named("attempts")=attempts);
}
// [[Rcpp::export]]
List update_graph_precision(NumericMatrix graph,NumericMatrix precision,double b,double b_star,
NumericMatrix D,NumericMatrix D_star,double graph_prior,
int iterations=500) {
  auto x=ordinal_mfm::graph_step(as<arma::mat>(graph),as<arma::mat>(precision),b,b_star,
  as<arma::mat>(D),as<arma::mat>(D_star),graph_prior,iterations);
  return List::create(Named("G")=x.graph,Named("K")=x.omega,
  Named("graph_accepted")=x.graph_accept,
  Named("precision_accepted")=x.precision_accept,Named("auxiliary_attempts")=x.auxiliary_attempts);
}
// Auxiliary-parameter allocation: refresh candidates for each observation.
// [[Rcpp::export]]
List update_allocations(NumericMatrix latent,double g,NumericVector VN,
IntegerVector z,IntegerVector K,IntegerVector N_k,
List graphs,List precision_matrices,List mean_vectors,
int n_grid,double g_prior,double a,int mc_samples=40,
bool random_scan=true) {
  using ordinal_mfm::Component;
  const int p=latent.ncol(),M=std::max(1,mc_samples),initial=K[0];
  if(n_grid!=latent.nrow()||n_grid<1||!(g>0&&a>0&&g_prior>0&&g_prior<1))
  stop("Invalid auxiliary-allocation parameters.");
  IntegerVector labels=clone(z);
  std::vector<Component> components;
  std::vector<int> counts;
  for(int k=0;k<initial;++k) {
    components.push_back({as<arma::mat>(graphs[k]),as<arma::mat>(precision_matrices[k]),as<arma::vec>(mean_vectors[k])});
    counts.push_back(N_k[k]);
  }
  IntegerVector order=seq(0,n_grid-1);
  if(random_scan)for(int i=n_grid-1;i>0;--i){
    int j=(int)R::runif(0,i+1);
    std::swap(order[i],order[j]);
  }
  arma::mat Z=as<arma::mat>(latent);
  int removed=0,born=0,total_candidates=0;
  double probsum=0,probmax=0;
  for(int ii=0;ii<n_grid;++ii) {
    int i=order[ii],old=labels[i]-1;
    Component previous=components[old];
    bool singleton=(counts[old]==1);
    labels[i]=0;
    --counts[old];
    if(singleton) {
      components.erase(components.begin()+old);
      counts.erase(counts.begin()+old);
      ++removed;
      for(int j=0;j<n_grid;++j)if(labels[j]>old+1)--labels[j];
    }
    const int t=components.size();
    std::vector<Component> candidates;
    // A singleton must retain its previous parameter in a uniformly chosen candidate slot.
    int recycle=singleton?(int)R::runif(0,M):-1;
    for(int h=0;h<M;++h){
      if(h==recycle)candidates.push_back(previous);
      else candidates.push_back(ordinal_mfm::prior_component(p,g_prior,a));
      ++total_candidates;
    }
    arma::vec zi=Z.row(i).t();
    std::vector<double> logw(t+candidates.size());
    for(int k=0;k<t;++k)logw[k]=std::log(counts[k]+g)+ordinal_mfm::loglike(zi,components[k]);
    double newbase=t?std::log(g)+VN[t]-VN[t-1]-std::log((double)M):-std::log((double)M);
    for(unsigned h=0;h<candidates.size();++h)
    logw[t+h]=newbase+ordinal_mfm::loglike(zi,candidates[h]);
    double maxw=*std::max_element(logw.begin(),logw.end()),sum=0,newsum=0;
    std::vector<double>w(logw.size());
    for(unsigned k=0;k<w.size();++k){
      w[k]=std::exp(logw[k]-maxw);
      sum+=w[k];
      if((int)k>=t)newsum+=w[k];
    }
    if(!std::isfinite(sum)||sum<=0)stop("Non-finite auxiliary allocation weights.");
    double np=newsum/sum;
    probsum+=np;
    probmax=std::max(probmax,np);
    double u=R::runif(0,sum),running=0;
    unsigned choice=w.size()-1;
    for(unsigned k=0;k<w.size();++k){
      running+=w[k];
      if(u<running){
        choice=k;
        break;
      }
    }
    if((int)choice<t){
      labels[i]=choice+1;
      ++counts[choice];
    }
    else {
      components.push_back(candidates[choice-t]);
      counts.push_back(1);
      labels[i]=t+1;
      ++born;
    }
  }
  int t=components.size();
  List gout(n_grid),oout(n_grid),mout(n_grid);
  IntegerVector nk(n_grid);
  for(int k=0;k<n_grid;++k){
    arma::mat G(p,p,arma::fill::zeros),O(p,p,arma::fill::eye);
    arma::vec mu(p,arma::fill::zeros);
    if(k<t){
      G=components[k].graph;
      O=components[k].omega;
      mu=components[k].mean;
      nk[k]=counts[k];
    }
    gout[k]=G;
    oout[k]=O;
    mout[k]=mu;
  }
  return List::create(Named("K")=t,Named("z")=labels,Named("N_k")=nk,
  Named("graphs")=gout,Named("precision_matrices")=oout,Named("mean_vectors")=mout,
  Named("initial_K")=initial,Named("removed_clusters")=removed,
  Named("new_cluster_draws")=born,Named("mean_new_cluster_prob")=probsum/n_grid,
  Named("max_new_cluster_prob")=probmax,
  Named("candidates_evaluated")=total_candidates);
}
// Sum V_N in log space until the relative omitted-tail bound is below tol.
// [[Rcpp::export]]
NumericVector log_mfm_coefficients(int n_grid,double gamma,double lambda,double tol=1e-12) {
  if(n_grid<1||!(gamma>0&&lambda>0&&tol>0&&tol<1))stop("Invalid V_N parameters.");
  NumericVector ans(n_grid+10),bounds(n_grid+10);
  IntegerVector limits(n_grid+10);
  int global_max=0;
  for(int t=1;t<=n_grid+10;++t){
    double acc=R_NegInf;
    bool done=false;
    for(int k=t;k<=100000;++k){
      double term=R::lgammafn(k+1.)-R::lgammafn(k-t+1.)+R::dpois(k-1.,lambda,true);
      for(int j=0;j<n_grid;++j)term-=std::log(k*gamma+j);
      double m=std::max(acc,term);
      acc=m+std::log(std::exp(acc-m)+std::exp(term-m));
      double ratio=lambda*(k+1.)/(k*(k+1.-t));
      if(ratio<1){
        double logbound=term+std::log(ratio)-std::log1p(-ratio)-acc;
        if(logbound<std::log(tol)){
          ans[t-1]=acc;
          bounds[t-1]=std::exp(logbound);
          limits[t-1]=k;
          global_max=std::max(global_max,k);
          done=true;
          break;
        }
      }
    }
    if(!done)stop("V_N relative tail bound did not reach requested tolerance.");
  }
  ans.attr("relative_tail_bounds")=bounds;
  ans.attr("per_t_Kmax")=limits;
  ans.attr("truncation_Kmax")=global_max;
  ans.attr("omitted_poisson_mass")=R::ppois(global_max-1.,lambda,false,false);
  ans.attr("tail_bound_tolerance")=tol;
  return ans;
}
// [[Rcpp::export]]
NumericMatrix update_latent(NumericMatrix current,NumericMatrix X,NumericMatrix theta,
IntegerVector labels,List means,List precisions,IntegerVector levels,int n_sweeps=3) {
  const int n=current.nrow(),p=current.ncol(),K=means.size(),L=theta.ncol()-1;
  if(X.nrow()!=n||X.ncol()!=p||theta.nrow()!=p||L<2||labels.size()!=n||
  precisions.size()!=K||levels.size()!=p||n_sweeps<1)stop("Invalid batched latent dimensions.");
  std::vector<NumericMatrix> Os;
  std::vector<NumericVector> mus;
  for(int k=0;k<K;++k){
    NumericMatrix O=as<NumericMatrix>(precisions[k]);
    NumericVector mu=as<NumericVector>(means[k]);
    if(O.nrow()!=p||O.ncol()!=p||mu.size()!=p)stop("Invalid cluster dimensions.");
    arma::mat M=as<arma::mat>(O),C;
    if(!M.is_finite()||arma::abs(M-M.t()).max()>1e-10||!arma::chol(C,M))stop("Invalid latent precision.");
    for(int j=0;j<p;++j)if(!std::isfinite(mu[j]))stop("Invalid latent mean.");
    Os.push_back(O);
    mus.push_back(mu);
  }
  NumericMatrix out=clone(current);
  for(int i=0;i<n;++i){
    if(i%32==0)checkUserInterrupt();
    int k=labels[i]-1;
    if(k<0||k>=K)stop("Invalid cluster label.");
    std::vector<double> lo(p),hi(p);
    for(int j=0;j<p;++j){
      double x=X(i,j);
      if(!std::isfinite(x)||x!=std::floor(x)||x<1||levels[j]<3||levels[j]>L||x>levels[j])stop("Invalid ordinal category.");
      lo[j]=theta(j,(int)x-1);
      hi[j]=theta(j,(int)x);
      if(!std::isfinite(out(i,j))||!(out(i,j)>lo[j]&&out(i,j)<hi[j]))stop("Latent input violates ordinal support.");
    }
    NumericMatrix O=Os[k];
    NumericVector mu=mus[k];
    for(int s=0;s<n_sweeps;++s)for(int j=0;j<p;++j){
      // Preserve the reference R arithmetic: round products to double before extended accumulation.
      long double acc=0;
      for(int l=0;l<p;++l)if(l!=j){
        volatile double difference=out(i,l)-mu[l];
        volatile double product=O(j,l)*difference;
        acc+=(long double)product;
      }
      double cm=mu[j]-(double)acc/O(j,j);
      out(i,j)=ordinal_mfm::tn(cm,std::sqrt(1.0/O(j,j)),lo[j],hi[j]);
    }
  }
  return out;
}
