# include <RcppArmadillo.h>
# include "CommonFunctions.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

//' Compares how similar two points are with regards to their clustering across 
//' all iterations. Works for unsupervised methods (i.e. allows label flipping).
//' 
//' @param point The index row within cluster_record for the first object.
//' @param comparison_point The comparison row index from cluster_record.
//' @param cluster_record A matrix of combined membership column vectors where 
//' the ijth entry is the membership assigned to the ith person in the jth 
//' iteration of sampling.
//' @param num_iter The number of samples recorded.
//' 
//' @return A score between 0 and 1 of the fraction of iterations for which the 
//' objects denoted by the point and comparison_point rows are assigned the same
//' label.
double calcPointSimilarity(arma::uword point, 
                           arma::uword comparison_point,
                           arma::umat cluster_record,
                           arma::uword num_iter
) {
  
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

//' Constructs a similarity matrix comparing all points clustering across the 
//' iterations.
//' 
//' @param cluster_record Matrix of label assignment for data across iterations.
//' 
//' @return A symmetric n x n matrix (for n rows in cluster record) describing 
//' the fraction of iterations for which each pairwise combination of points are
//' assigned the same label.
//' @export
// [[Rcpp::export]]
arma::mat createSimilarityMat(arma::umat cluster_record){
  
  // Sample size
  arma::uword n = cluster_record.n_rows;
  
  // Number of iterations from MCMC
  arma::uword n_iter = cluster_record.n_cols;
  
  // Output matrix (n x n similarity matrix)
  arma::mat out = arma::ones<arma::mat>(n,n);
  
  // Current value
  double entry = 0.0;
  
  // Compare every entry to every other entry. As symmetric and diagonal is I
  // do not need to compare points with self and only need to calcualte (i, j) 
  // entry
  for (arma::uword i = 0; i < n - 1; i++){ 
    for (arma::uword j = i + 1; j < n; j++){
      entry = calcPointSimilarity(i, j, cluster_record, n_iter);
      out(i, j) = entry;
      out(j, i) = entry;
    }
  }
  return out;
}
