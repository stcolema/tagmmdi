# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// --- Hyperparameters ---------------------------------------------------------

// Returns a variable involved in updating the scale, it is similar to sample 
// covariance
arma::mat S_n(arma::mat data,
              arma::vec sample_mean,
              arma::uword n,
              arma::uword n_col){
  
  arma::mat sample_covariance = arma::zeros<arma::mat>(n_col, n_col);
  // sample_covariance.zeros();
  
  // If n > 0 (as this would crash for empty clusters)
  if(n > 0){
    data.each_row() -= arma::trans(sample_mean);
    sample_covariance = arma::trans(data) * data;
  }
  return sample_covariance;
}

// Returns a vector of the mean of a cluster of size n
arma::vec mean_n(double lambda_0,
                 arma::vec mu_0,
                 arma::uword n,
                 arma::uword n_col,
                 arma::vec sample_mean){
  
  arma::vec mu_n(n_col);
  
  mu_n = (lambda_0 * mu_0 + n * sample_mean) / (double)(lambda_0 + n);
  return mu_n;
}

// Returns the matrix for the scale hyperparameter for n observations
arma::mat calc_scale_n(arma::mat scale_0, 
                       arma::mat sample_cov, 
                       double lambda, 
                       int n, 
                       arma::vec sample_mean, 
                       arma::vec mu_0,
                       int d){
  
  arma::mat dist_from_prior(d, d);
  arma::mat scale_n(d, d);
  
  // Calculate the distance of the sample mean from the prior
  dist_from_prior = (sample_mean - mu_0) * arma::trans(sample_mean - mu_0);
  
  // Update the scale hyperparameter
  scale_n = scale_0 + sample_cov + ((lambda * n) / (double)(lambda + n)) * dist_from_prior;
  
  return scale_n;
}


// ---- This is not use --------------------------------------------------------

// Returns the matrix for the scale hyperparameter for n observations
arma::mat scale_n(arma::mat scale_0,
                  arma::vec mu_0,
                  double lambda_0,
                  arma::mat sample_covariance,
                  arma::uword n,
                  arma::vec sample_mean){
  
  arma::mat scale_n;
  double lambda_n = ((lambda_0 * n) / (double)(lambda_0 + n));
  
  if(n > 0){
    // The posterior of the scale parameter
    scale_n = (scale_0
                 + sample_covariance
                 + lambda_n * (sample_mean - mu_0) * arma::trans(sample_mean - mu_0)
    );
    
    return scale_n;
  }
  
  // If there's no data present we retain the prior belief undiluted
  scale_n = scale_0; 
  return scale_out;
}

// --- Variance and mean posteriors --------------------------------------------

arma::mat variance_posterior(int df_0,
                             arma::mat scale_0,
                             double lambda,
                             arma::vec mu_0,
                             arma::mat data
){
  
  // The sample size
  int n = data.n_rows;
  
  // The dimension of the data
  int d = data.n_cols;

  // Declare the degrees of freedom and scale parameters, initially setting to
  // the prior values (in case n = 0 and there's no update to our belief about 
  // them)
  double df_n = df_0;
  arma::mat scale_n(d, d);
  arma::vec sample_mean(d);
  arma::mat sample_cov(d, d);
  
  scale_n = scale_0;
  
  // if any points belong to this data update the priors
  if(n > 0){
    // The posterior degrees of freedom
    df_n = df_0 + n;
    
    if(n > 1){
      // Calculate the sample mean (the 0 means the mean of the columns is calculated)
      sample_mean = arma::trans(arma::mean(data, 0));
      
    } else {
      // If only one point in cluster the mean is the entry
      sample_mean = arma::trans(data.row(0));
    }
    
    // Sample covariance without division by sample size
    sample_cov = S_n(data, sample_mean, n, d);
    
    // Calculate the scale variable posterior
    scale_n = calc_scale_n(scale_0, sample_cov, lambda, n, sample_mean, mu_0, d);
  }
  
  // Draw the current variance from the inverse wishart
  arma::mat variance(d, d);
  variance = arma::iwishrnd(scale_n, df_n);
  
  return variance;
}

// sample from the mean posterior
arma::vec mean_posterior(arma::vec mu_0,
                         arma::mat variance,
                         double lambda_0,
                         arma::mat data){
  
  arma::uword n_col = data.n_cols;
  arma::vec sample_mean = arma::zeros<arma::vec>(n_col);
  arma::vec mu_out = arma::zeros<arma::vec>(n_col);
  arma::vec mu_n = arma::zeros<arma::vec>(n_col);
  arma::mat variance_n = arma::zeros<arma::mat>(n_col, n_col);
  arma::uword n = data.n_rows;
  
  if (n > 0){
    sample_mean = trans(arma::mean(data, 0));
  }
  
  // The lambda hyperparameter posterior value
  double lambda_n = lambda_0 + n;
  
  // Update mean hyperparameter
  mu_n = mean_n(lambda_0, mu_0, n, n_col, sample_mean);
  
  // Update variance hyperparameter
  variance_n = variance / lambda_n;
  
  // Draw from multivariate Normal distribution
  arma::mat x = arma::mvnrnd(mu_n, variance_n, 1);
  mu_out = x.col(0);
  
  return mu_out;
  
}