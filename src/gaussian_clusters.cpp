# include <RcppArmadillo.h>
// # include <iostream>
// # include <fstream>
# include "common_functions.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

//' Returns the normal distribution log-likelihood.
//' 
//' @param point Column vector describing the measurements for a sample of data.
//' @param mu The mean for the data (same dimensionality as point).
//' @param variance The covariance matrix.
//' @param d The dimensionality of point.
//' 
//' @return The log likelihood of the multivariate normal distribution for an
//' observation.
double CalcNormalLikelihood(arma::vec point,
                            arma::vec mu,
                            arma::mat variance,
                            arma::uword d
) {
  double log_det = 0.0;
  double exponent = 0.0;
  double log_likelihood = 0.0;
  
  // Exponent in normal PDF
  exponent = arma::as_scalar( (point - mu).t()
                             * arma::inv(variance)
                             * (point - mu) );
                               
  // Log determinant of the variance
  log_det = arma::log_det(variance).real();
  
  // Normal log likelihood
  log_likelihood = -0.5 *(log_det + exponent + (double) d * log(2.0 * M_PI)); 
  
  return log_likelihood;
}

//' For k classes returns a k-vector of probabilities for point to belong to said
//' classes.
//' 
//' @param point The column vector of measurements for a given sample.
//' @param data A matrix of the entire dataset.
//' @param k The number of clusters present.
//' @param class_weights The k-vector of cluster weights.
//' @param mu d x k matrix of mean vectors for each cluster.
//' @param variance k-cube of d x d covariance matrices for each cluster.
arma::vec SampleGaussianMembership(arma::vec point,
                                   arma::mat data,
                                   arma::uword k,
                                   arma::vec class_weights,
                                   arma::mat mu,
                                   arma::cube variance
) {
  
  double curr_weight = 0.0;
  double log_likelihood = 0.0;
  
  arma::vec prob_vec(k);
  prob_vec.zeros();
  
  arma::uword d = data.n_cols;
  
  // Calculate the likelihood for each class
  for(arma::uword i = 0; i < k; i++){
    curr_weight = log(class_weights(i));
    
    log_likelihood = CalcNormalLikelihood(point, 
                                          mu.col(i),
                                          variance.slice(i),
                                          d);
    
    prob_vec(i) = curr_weight + log_likelihood;
  } 
  return prob_vec;
}

//' Predicts cluster membership.
//' 
//' @param my_vec Vector of cluster assignment probabilities.
//' 
//' @return Label representing predicted cluster membership.
arma::uword PredictIndex(arma::vec my_vec) {
  
  double u = 0.0;
  arma::uword pred_ind = 0;
  
  // Overflow handling and convert from logs
  my_vec = exp(my_vec - max(my_vec));
  
  // Normalise the vector
  my_vec = my_vec / sum(my_vec);
  
  // Sample from uniform distribution to select which value to use
  u = arma::randu<double>( );
  
  // include + 1 if labels begin at 1
  pred_ind = sum(u > cumsum(my_vec));
  
  return pred_ind;
  
}

//' Converts a vector of log, unnormalised probabilities to probabilities.
//' @param my_log_vec A vector of log, unnormalised porbabilities.
//' 
//' @return The exponentiated, normalised transform of my_log_vec.
arma::vec HandleOverflow(arma::vec my_log_vec) {
  
  arma::uword p = my_log_vec.n_cols;
  arma::vec my_vec(p);
  
  // Overflow handling and convert from logs
  my_vec = exp(my_log_vec - max(my_log_vec) );
  
  // Normalise the vector
  my_vec = my_vec / sum(my_vec);
  
  return(my_vec);
}

// Predicts the cluster assignments based on a vector of probabilities using
// the rejection method
// Old name: cluster_predictor
// [[Rcpp::export]]
arma::uword PredictClusterMembership(arma::vec probabilities) {
  double u;
  arma::uword pred;
  u = arma::randu<double>( );
  
  // include + 1 if labels centred on 1
  pred = 1 + sum(u > cumsum(probabilities)); 
  return pred;
}