#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//calculate VN
// [[Rcpp::export]]
NumericVector fn_calc_VN(int n_grid, double gamma, double lambda) {
  NumericVector VN(n_grid + 10L, 0);
  double b, r, m;
  for (int t = 1; t < n_grid + 11; t++) {
    
    
    r = R_NegInf;
    for (int k = t; k < n_grid + 101; k++) {
      
      NumericVector l1 = log(Range(k - t + 1, k));
      NumericVector l2 = log(Range(k * gamma, k * gamma + n_grid - 1));
      b = sum(l1) - sum(l2) + dpois(NumericVector::create(k - 1), lambda, true)[0];
      
      m = max(NumericVector::create(b, r));
      r = log(exp(r - m) + exp(b - m)) + m;
    }
    VN[t - 1] = r;
  }
  return(VN);
}



// [[Rcpp::export]]
void fn_remove_grid(int id_grid, IntegerVector z, IntegerVector K, IntegerVector N_k,
                    List graphs, List precision_matrices, List mean_vectors,
                    int n_grid){

  int cur_z, i;
  
  cur_z = z[id_grid];
  z[id_grid] = (-1);
  
  --N_k[cur_z - 1]; 
  
  
  //If a grid is clustered by itself
  if (N_k[cur_z - 1] == 0) {
    for (i = 0; i < n_grid; i++) {
      if (z[i] > cur_z) {
        --z[i];
      }
    }
    
    for (i = cur_z; i < n_grid; i++) {
      N_k[i - 1] = N_k[i];
      graphs[i - 1] = graphs[i];
      precision_matrices[i - 1] = precision_matrices[i];
      mean_vectors[i - 1] = mean_vectors[i];
      
    }
    
    
    int p=as<arma::mat>(precision_matrices[0]).n_cols;
    vec zero_vector=zeros(p);
    mat zero_matrix = zeros(p,p);
    mat diag_matrix;
    diag_matrix.eye(p, p);
    
    
    N_k[(n_grid - 1)] = 0;
    
    graphs[(n_grid - 1)] = wrap(zero_matrix);
    precision_matrices[(n_grid - 1)] = wrap(diag_matrix); // Set to zero matrix
    mean_vectors[(n_grid - 1)] = wrap(zero_vector); // Set to zero vector
    
    --K[0];
  }	
}



// [[Rcpp::export]]
void fn_z_update_mod(NumericMatrix latent, double g, NumericVector VN,
                     IntegerVector z, IntegerVector K, IntegerVector N_k,
                     List graphs, List precision_matrices, List mean_vectors,
                     int n_grid,double g_prior,double a){
 
  NumericVector prob(n_grid + 0L, 0);
  int p=as<arma::mat>(precision_matrices[0]).n_cols;
  
  
  NumericVector mean_k;
  NumericMatrix precision_k;
  NumericVector data_i;
  int i, k;
  arma::mat latent_arma=as<arma::mat>(latent);
  Function dmvnorm("dmvnorm");//package mvtnorm
  Function gnorm("gnorm");
  Function rgwish("rgwish");
  for (i = 0; i < n_grid; i++)
  {
    
    fn_remove_grid(i, z, K, N_k, graphs,precision_matrices,mean_vectors,n_grid);
    std::fill(prob.begin(), prob.end(), 0);
    data_i = wrap((latent_arma).row(i));
    
    
    for (k = 0; k < K[0]; k++) {
      mean_k = mean_vectors[k];
      arma::mat precision_k = as<arma::mat>(precision_matrices[k]);
      arma::mat covariance_k = arma::inv(precision_k);
      NumericMatrix covariance_k_r = wrap(covariance_k);
      
      
      NumericVector log_density =(dmvnorm(data_i, Named("mean", mean_k), Named("sigma", covariance_k_r),Named("log", true)));
      
      prob[k] = log(g + N_k[k]) + log_density[0];
      
    }
    
    
    
    int b=3;
    arma::mat D = arma::eye(p, p);
    
    arma::mat D_tilde =D + (arma::vec(data_i))*(arma::vec(data_i)).t()*(1-1/(a+1));
   
    double mc_estimate = 0.0;
    
    int mc_samples=15;//minus
    
    int valid_samples = 0; 
    
    
    
   
    for (int iii = 0; iii < mc_samples; iii++) {
      arma::mat adj_matrix = arma::zeros<arma::mat>(p, p);
     
      for (int row = 0; row < p; row++) {
        for (int col = row + 1; col < p; col++) {
          adj_matrix(row, col) = R::rbinom(1, g_prior);
        }
      }
      //Function gnorm("gnorm");
      double G_norm1 = as<double>(gnorm(Named("adj", adj_matrix), Named("b", b),Named("D", D)));           // b = 3.0 as per your R code
      double G_norm2 = as<double>(gnorm(adj_matrix, b+1, D_tilde));
      
      
      
      if (std::isnan(G_norm2)) {
        continue; 
      }
      
      
      mc_estimate += std::exp(G_norm2 - G_norm1);
      valid_samples++;
    }
    mc_estimate /= valid_samples;
    prob[K[0]] = log(g) + VN[K[0]] - VN[K[0] - 1]
    - (p / 2.0) * std::log(2 * M_PI) + std::log(mc_estimate)+(p / 2.0) * std::log(a/(a+1));
    
    IntegerVector range = seq(0, K[0]);
    NumericVector tprob = prob[range];
    
    tprob = exp(tprob);
    
    IntegerVector ind = Rcpp::RcppArmadillo::sample(range, 1, true, tprob);
    
    k = ind[0];
    z[i] = k + 1;
    
    
    if (k < K[0]) { //Join one existing cluster
      ++N_k[k];
    } else { //Start a new cluster
      N_k[k] = 1;
      
      arma::mat G = arma::randu<arma::mat>(p, p);
      G.transform( [&](double val) { return R::rbinom(1, g_prior); } ); 
      G = arma::symmatu(G);
      G.diag().zeros();  
      
      arma::mat D = arma::eye(p, p);  
      
      //Function rgwish("rgwish");
      arma::mat precision_matrix = as<arma::mat>(rgwish(1, G, b, D));
      arma::vec mean_vector = arma::zeros(p);
      graphs[k]= wrap(G);
      precision_matrices[k]=wrap(precision_matrix); 
      mean_vectors[k]=wrap(arma::trans(mean_vector));
      ++K[0];
    }
    
    
  }}
