# include <RcppArmadillo.h>
# include "CommonFunctions.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

//' Converts a vector of log, unnormalised probabilities to probabilities.
//' @param my_log_vec A vector of log, unnormalised porbabilities.
//' 
//' @return The exponentiated, normalised transform of my_log_vec.
arma::vec makeProbabilities(arma::vec my_log_vec) {
  
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
//' @param prob_vec Vector of cluster assignment probabilities.
//' @param start_label The value clusters are counted from (typcially 0 in C++ 
//' and 1 in R).
arma::uword predictIndex(arma::vec probabilities, arma::uword start_label) {
  double u;
  arma::uword pred;
  u = arma::randu<double>( );
  
  // include + 1 if labels centred on 1
  pred = sum(u > cumsum(probabilities)) + start_label; 
  return pred;
}

//' Predicts cluster membership.
//' 
//' @param prob_vec Vector of cluster assignment scores (unnormalised, log 
//' probabilities).
//' @param start_label The value clusters are counted from (typcially 0 in C++ 
//' and 1 in R).
//' @return Label representing predicted cluster membership.
arma::uword predictCluster(arma::vec my_vec, arma::uword start_label) {
  
  double u = 0.0;
  arma::uword pred = 0;
  
  // Overflow handling and convert from assignment scores to probabilities
  my_vec = makeProbabilities(my_vec);
  
  // Sample from uniform distribution to select which value to use
  u = arma::randu<double>( );
  
  // Predict the component membership
  pred = predictIndex(my_vec, start_label);
  
  return pred;
}