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
double calcNormalLikelihood(arma::vec point,
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
arma::vec calcGaussianMembership(arma::vec point,
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
    
    log_likelihood = calcNormalLikelihood(point, 
                                          mu.col(i),
                                          variance.slice(i),
                                          d);
    
    prob_vec(i) = curr_weight + log_likelihood;
  } 
  return prob_vec;
}
