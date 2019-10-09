# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// Returns the normal distribution log-likelihood
double CalcNormalLikelihood(arma::vec point,
                         arma::vec mu,
                         arma::mat variance,
                         arma::uword d){
  double log_det = 0.0;
  double exponent = 0.0;
  double log_likelihood = 0.0;
  
  // Exponent in normal PDF
  exponent = arma::as_scalar(arma::trans(point - mu) 
                               * arma::inv(variance)
                               * (point - mu));
                               
  // Log determinant of the variance
  log_det = arma::log_det(variance).real();
  
  // Normal log likelihood
  log_likelihood = -0.5 *(log_det + exponent + (double) d * log(2.0 * M_PI)); 
  
  return log_likelihood;
}

// For k classes returns a k-vector of probabilities for point to belong to said
// classes
arma::vec SampleGaussianMembership(arma::vec point,
                                  arma::mat data,
                                  arma::uword k,
                                  arma::vec class_weights,
                                  arma::mat mu,
                                  arma::cube variance){
  
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

arma::uword PredictIndex(arma::vec my_vec){
  
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