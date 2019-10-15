// # include <RcppArmadillo.h>
// # include <iostream>
// # include <fstream>
// 
// // [[Rcpp::depends(RcppArmadillo)]]
// 
// using namespace Rcpp ;
// 
// // count unique entries in a vector (Armadillo has this already - unnecssary)
// arma::uword unique_counter(arma::uvec v){
//   std::sort(v.begin(), v.end());
//   arma::uword unique_count = std::unique(v.begin(), v.end()) - v.begin();
//   return unique_count;
// }
// 
// // returns a vector of the number of unqiue values in each column
// // [[Rcpp::export]]
// arma::uvec cat_counter(arma::umat data){
//   arma::uword num_cols = data.n_cols;
//   arma::uvec num_categories(num_cols);
//   for(arma::uword i = 0; i < num_cols; i++){
//     num_categories(i) = unique_counter(data.col(i));
//   }
//   return num_categories;
// }
// 
// // find the number of categories in each covariate and declare the appropriate
// // matrix to record the associated probabilties for each cluster
// // [[Rcpp::export]]
// arma::field<arma::mat> declare_class_probs_field(arma::uvec cat_count,
//                                                  arma::uword num_cols,
//                                                  arma::uword num_clusters){
//   
//   arma::field<arma::mat> class_probabilities(num_cols);
//   for(arma::uword i = 0; i < num_cols; i++){
//     arma::mat phi_j(num_clusters, cat_count(i));
//     phi_j.zeros();
//     class_probabilities(i) = phi_j;
//   }
//   return class_probabilities;
// }
// 
// // Sample the probabilities for each category across all clusters for each covariate
// // [[Rcpp::export]]
// arma::field<arma::mat> sample_class_probabilities(arma::umat data,
//                                                   arma::field<arma::mat> class_probs,
//                                                   arma::field<arma::vec> phi_prior,
//                                                   arma::uvec cluster_labels,
//                                                   arma::uvec cat_count,
//                                                   arma::uword num_clusters,
//                                                   arma::uword num_cols
// ){
//   
//   arma::umat cluster_data;
//   for(arma::uword k = 1; k < num_clusters + 1; k++){
//     
//     cluster_data = data.rows(find(cluster_labels == k));
//     
//     for(arma::uword j = 0; j < num_cols; j++){
//       
//       
//       class_probs(j).row(k - 1) = arma::trans(dirichlet_posterior_class(phi_prior(j),
//                                         cluster_data.col(j),
//                                         cat_count(j)
//       )
//       );
//       
//     }
//   }
//   return class_probs;
// }