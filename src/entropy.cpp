# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// Calculate the entropy for the current cluster weights
// [[Rcpp::export]]
double CalcEntropy(arma::vec class_weights){
  
  // Declare objects
  arma::uword n = class_weights.n_elem;
  arma::vec entropy_components(n);
  double entropy_out = 0.0;
  
  // Calculate the entropy
  entropy_out = - sum(class_weights.t() * log(class_weights));

  return entropy_out;
}