# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// Compares how similar two points are with regards to their clustering across 
// all iterations.
// Works for unsupervised methods (i.e. allows label flipping)
double point_similarity(arma::uword point, 
                        arma::uword comparison_point,
                        arma::umat cluster_record,
                        arma::uword num_iter) {
  
  // Declare objects
  double out = 0.0;
  arma::urowvec ind_1(num_iter);
  arma::urowvec ind_2(num_iter);
  arma::umat out_1(1, num_iter);
  
  // The relvant cluster vectors
  ind_1 = cluster_record.row(point);
  ind_2 = cluster_record.row(comparison_point);
  
  // Compare vector of allocations element-wise
  out_1.row(0) = (ind_1 == ind_2);
  
  // Similarity is the sum of the above divided by the number of entries
  // Convert the sum to a double as otherwise is integer divison and does not 
  // work
  out = (double)arma::sum(out_1.row(0)) / (double)num_iter;
  
  return out;
}

// Constructs a similarity matrix comparing all points clustering across the 
// iterations
// [[Rcpp::export]]
arma::mat similarity_mat(arma::umat cluster_record){
  
  // Sample size
  arma::uword n = cluster_record.n_rows;
  
  // Number of iterations from MCMC
  arma::uword n_iter = cluster_record.n_cols;
  
  // Output matrix (n x n similarity matrix)
  // arma::mat out(n, n);
  arma::mat out = arma::ones<arma::mat>(n,n);
  
  // Compare every entry to every other entry. As symmetric and diagonal is I
  // do not need to compare points with self and only need to calcualte (i, j) 
  // entry
  for (arma::uword i = 0; i < n - 1; i++){ 
    for (arma::uword j = i + 1; j < n; j++){
      out(i, j) = point_similarity(i, j, cluster_record, n_iter);
      out(j, i) = out(i, j);
    }
  }
  return out;
}
