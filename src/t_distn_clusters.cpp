# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;




// Returns the t-distribution log likelihood for a given point
double CalcTdistnLikelihood(arma::vec point,
                            arma::vec mu,
                            arma::mat variance,
                            arma::uword n_col,
                            double df
){
  
  double log_det = 0.0;
  double exponent = 0.0;
  double log_likelihood = 0.0;
  
  exponent = arma::as_scalar(
    arma::trans(point - mu) 
    * arma::inv(variance)
    * (point - mu)
  );
  
  log_det = arma::log_det(variance).real();
  
  log_likelihood = lgamma((df + n_col)/2.0) 
    - lgamma(df/2.0) 
    - (double)n_col/2.0 * log(df * M_PI) 
    - 0.5 * log_det 
    - ((df + (double)n_col)/2.0) * log(1.0 + (1.0/df) * exponent);
    
    return log_likelihood;
}


// Samples from a Beta distribution based on the idea that two independent gamma
// functions can be used. 
// See https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions.
// Used in calculating the weights for allocation to outlier
double SampleBetaDistn(double a, double b, double theta = 1.0){
  double X = arma::randg( arma::distr_param(a, 1/theta) );
  double Y = arma::randg( arma::distr_param(b, 1/theta) );
  double beta = X / (double)(X + Y);
  return beta;
}


// Samples if the point is an outlier or not (comparing to assigned class)
double SampleOutlier(arma::vec point,
                      arma::mat data,
                      double outlier_weight,
                      arma::vec global_mean,
                      arma::mat global_variance,
                      double t_df = 4.0){
  
  arma::uword n = data.n_rows;
  arma::uword d = data.n_cols;
  double log_likelihood = 0.0;
  double prob = 0.0;
  
  // The probability of belonging to the outlier class
  log_likelihood = t_likelihood(point, global_mean, global_variance, d, t_df);
  
  prob = log(outlier_weight) + log_likelihood;
  
  return prob;
}

// Calculate the probability of being an outlier
arma::vec CalculateOutlierProb(arma::vec point,
                                 arma::vec global_mean,
                                 arma::mat global_variance,
                                 arma::uword n_col,
                                 double t_df,
                                 double outlier_weight,
                                 arma::col mu,
                                 arma::mat var){
  double out_likelihood = 0.0;
  double in_likelihood = 0.0;
  double non_outlier_weight = 1.0 - outlier_weight;
  arma::vec outlier_prob(2);
  outlier_prob.zeros();
  
  // Calculate outlier likelihood
  out_likelihood = t_likelihood(point, global_mean,  global_variance, n_col, t_df);
  out_likelihood += log(outlier_weight);
  
  outlier_prob(1) = out_likelihood;
  
  // Calculate likelihood of belonging to current cluster
  in_likelihood = normal_likelihood(point, mu, var, n_col);
  
  // Update with the non-outlier weight
  in_likelihood += log(non_outlier_weight);
  
  outlier_prob(0) = in_likelihood;
  
  return outlier_prob;
}