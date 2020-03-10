# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

//' Calculate the entropy for the current cluster weights.
//' 
//' @param class_weights The vector of the cluster weights from a mixture model.
//' 
//' @return The information entropy for the current cluster weights.
// [[Rcpp::export]]
double calcEntropy(arma::vec class_weights) {
  
  // Declare objects
  arma::uword n = class_weights.n_elem;
  arma::vec entropy_components(n);
  double entropy_out = 0.0;
  
  // Calculate the entropy
  entropy_out = - sum(class_weights.t() * log(class_weights));

  return entropy_out;
}