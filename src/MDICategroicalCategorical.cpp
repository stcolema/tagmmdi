# include <RcppArmadillo.h>
// # include <iostream>
// # include <fstream>
# include "common_functions.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// Refers to categorical_parameters.cpp, mdi_parameters.cpp



// MDI clustering for two categorical datasets
// [[Rcpp::export]]
Rcpp::List mdiCatCat(arma::umat data_1,
                     arma::umat data_2,
                     arma::field<arma::vec> class_dist_prior_1,
                     arma::field<arma::vec> class_dist_prior_2,
                     arma::vec clust_weight_priors_1,
                     arma::vec clust_weight_priors_2,
                     arma::uvec clust_labels_1,
                     arma::uvec clust_labels_2,
                     arma::uword n_clust_1,
                     arma::uword n_clust_2,
                     arma::uvec fix_vec_1,
                     arma::uvec fix_vec_2,
                     double a_0,
                     double b_0,
                     arma::uword num_iter,
                     arma::uword burn,
                     arma::uword thinning
) {
  // As we tend to iterate to less than num_iter, add 1 to itself
  // Done here as many objects are sized based on this.
  num_iter++;
  
  // Declare the sample size and dimensionality of the continuous and 
  // categorical data
  arma::uword n = data_1.n_rows;
  arma::uword n_cols_1 = data_1.n_cols;
  arma::uword n_cols_2 = data_2.n_cols;
  arma::uword min_n_clust = std::min(n_clust_1, n_clust_2); // for comparing across datasets
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  arma::uword v_a_0 = 1;
  arma::uword v_b_0 = 1;
  arma::uword cpp_class_start = 0;
  arma::uword r_class_start = 1;
  double v = 0.0; // strategic latent variable
  double phi = arma::randg(arma::distr_param(a_0, 1/b_0) ); // Context similarity
  double Z = 1.0; // Normalising constant
  arma::uvec cat_count_1(n_cols_1);
  arma::uvec cat_count_2(n_cols_2);
  arma::vec clust_weights_1 = arma::zeros<arma::vec>(n_clust_1);
  arma::vec clust_weights_2 = arma::zeros<arma::vec>(n_clust_2);
  arma::vec curr_prob_vec_1(n_clust_1);
  arma::vec curr_prob_vec_2(n_clust_2);
  arma::vec phi_record(eff_count);
  arma::vec labels_weights_phi(n + n_clust_2 + 1);
  arma::vec entropy_cw(num_iter);
  arma::vec rate_0_1(n_clust_1);
  arma::vec rate_0_2(n_clust_2);
  arma::umat record_1(n, eff_count);
  arma::umat record_2(n, eff_count);
  arma::mat sim_1(n, n); 
  arma::mat sim_2(n, n);
  arma::field<arma::mat> class_prob_1(n_cols_1);
  arma::field<arma::mat> class_prob_2(n_cols_2);
  
  // Declare the field for the phi variable for the categorical data
  cat_count_1 = countCategories(data_1);
  
  class_prob_1 = declareClassProbsField(cat_count_1,
                                        n_cols_1,
                                        n_clust_1);
  
  
  
  cat_count_2 = countCategories(data_2);
  class_prob_2 = declareClassProbsField(cat_count_2,
                                        n_cols_2,
                                        n_clust_2);
  
  // Placeholder prior
  rate_0_1.fill(1);
  rate_0_2.fill(1);
  
  // Initialise v based on the prior
  v = arma::randg( arma::distr_param(v_a_0, 1.0/(v_b_0) ) );
  
  record_1.zeros();
  record_2.zeros();
  
  for(arma::uword i = 0; i < num_iter; i++){
  
    // sample cluster weights for the two datasets
    clust_weights_1 = sampleMDIClusterWeights(clust_weight_priors_1,
                                              rate_0_1,
                                              v,
                                              n_clust_1,
                                              n_clust_2,
                                              clust_weights_2,
                                              clust_labels_1,
                                              clust_labels_2,
                                              phi);
    
    clust_weights_2 = sampleMDIClusterWeights(clust_weight_priors_2,
                                          rate_0_2,
                                          v,
                                          n_clust_2,
                                          n_clust_1,
                                          clust_weights_1,
                                          clust_labels_2,
                                          clust_labels_1,
                                          phi);
    
    // Entropy for graphing convergence
    entropy_cw(i) = calcEntropy(clust_weights_1);
    
    // For the categorical data, sample the probabilities for each class
    class_prob_1 = sampleCategoryProbabilities(data_1,
                                              class_prob_1,
                                              class_dist_prior_1,
                                              clust_labels_1,
                                              cat_count_1,
                                              n_clust_1,
                                              n_cols_1,
                                              r_class_start);
    
    class_prob_2 = sampleCategoryProbabilities(data_2,
                                              class_prob_2,
                                              class_dist_prior_2,
                                              clust_labels_2,
                                              cat_count_2,
                                              n_clust_2,
                                              n_cols_2,
                                              r_class_start);
    
    
    // Calculate the current normalising constant (consider being more clever 
    // about this) 
    Z = calcNormalisingConst(clust_weights_1,
                             clust_weights_2,
                             phi,
                             n_clust_1,
                             n_clust_2);
    
    // sample the strategic latent variable, v
    v = arma::randg( arma::distr_param(v_a_0 + n, 1.0/(v_b_0 + Z) ) );
    
    // sample the context similarity parameter (as only two contexts this is a
    // single double - number not a drink)
    phi = sampleMDIPhi(clust_labels_1,
                       clust_labels_2,
                       clust_weights_1,
                       clust_weights_2,
                       v,
                       n,
                       min_n_clust,
                       a_0,
                       b_0);
    
    // sample 
    for(arma::uword j = 0; j < n; j++){
      
      // for each point create the vector of probabilities associated with 
      // assignment to each cluster
      
      curr_prob_vec_1 = sampleMDICatClustProb(j, 
                                              data_1,
                                              class_prob_1,
                                              n_clust_1,
                                              n_cols_1,
                                              phi,
                                              clust_weights_1,
                                              clust_labels_1,
                                              clust_labels_2);
      
      curr_prob_vec_2 = sampleMDICatClustProb(j, 
                                              data_2,
                                              class_prob_2,
                                              n_clust_2,
                                              n_cols_2,
                                              phi,
                                              clust_weights_2,
                                              clust_labels_2,
                                              clust_labels_1);
      
      // update labels - in gaussian data this is only if the current point is 
      // not fixed
      if(fix_vec_1[j] == 0){
        clust_labels_1(j) = predictCluster(curr_prob_vec_1, 1); // 1 moves to R labelling
      }
      
      if(fix_vec_2[j] == 0){
        clust_labels_2(j) = predictCluster(curr_prob_vec_2, 1);
      }
    }
    
    // Update cluster labels in second dataset
    // Return the new labels, weights and similarity in a single vector
    labels_weights_phi = updateClusterLabels(clust_labels_1,
                                             clust_labels_2,
                                             clust_weights_1,
                                             clust_weights_2,
                                             n_clust_1,
                                             n_clust_2,
                                             phi,
                                             min_n_clust,
                                             v,
                                             n,
                                             a_0,
                                             b_0,
                                             Z);
    
    // Separate the output into the relevant components
    clust_labels_2 = arma::conv_to<arma::uvec>::from(labels_weights_phi.subvec(0, n - 1));
    clust_weights_2 = labels_weights_phi.subvec(n, n + n_clust_2 - 1);
    phi = arma::as_scalar(labels_weights_phi(n + n_clust_2));
    
    // if current iteration is a recorded iteration, save the labels
    if (i >= burn && (i - burn) % thinning == 0) {
      record_1.col((i - burn) / thinning) = clust_labels_1;
      record_2.col((i - burn) / thinning) = clust_labels_2;
      phi_record((i - burn) / thinning) = phi;
    }
  }
  
  // construct similarity matrix
  sim_1 = createSimilarityMat(record_1);
  sim_2 = createSimilarityMat(record_2);
  
  // std::cout << "Context similarity: " << context_similarity << "\n";
  
  return List::create(Named("similarity_1") = sim_1,
                      Named("similarity_2") = sim_2,
                      Named("class_record_1") = record_1,
                      Named("class_record_2") = record_2,
                      Named("context_similarity") = phi_record,
                      Named("entropy") = entropy_cw);
  
}
