# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// The actual categorical clustering all wrapped up in one function
// [[Rcpp::export]]
// Rcpp::List
arma::umat  categorical_clustering(arma::umat data,
                                  arma::field<arma::vec> phi_prior,
                                  arma::uvec cluster_labels,
                                  arma::uvec fix_vec,
                                  arma::vec cluster_weight_priors,
                                  arma::uword num_clusters,
                                  arma::uword num_iter,
                                  arma::uword burn,
                                  arma::uword thinning){
  
  // To allow using < and keeping in line with object sizes
  num_iter++;
  
  // Find the dimensions of the dataset
  arma::uword n = data.n_rows;
  arma::uword num_cols = data.n_cols;
  
  // Find the number of categories in each column
  arma::uvec cat_count(num_cols);
  cat_count = cat_counter(data);
  
  // Declare the object that will hold the pobability of assignment to each 
  // class for each column
  arma::field<arma::mat> class_probabilities(num_cols);
  class_probabilities = declare_class_probs_field(cat_count,
                                                  num_cols,
                                                  num_clusters);
  
  // Cluster weights
  arma::vec cluster_weights(num_clusters);
  
  // Vector to hold probability of assignment for current individual
  arma::vec curr_cluster_probs(num_clusters);
  
  // Number of samples recorded
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  
  // The matrix of labels to record
  arma::umat record(n, eff_count);
  record.zeros();
  
  // // The posterior similarity matrix (I think this should be outside this function)
  arma::mat sim(n, n);
  sim.zeros()
  
  // Iterate over the desired number of iterations
  for(arma::uword i = 0; i < num_iter; i++){
    
    // Sample cluster weights from a Dirichlet distribution
    cluster_weights = dirichlet_posterior(cluster_weight_priors,
                                          cluster_labels,
                                          num_clusters);
    
    
    // sample probabilities for each class for each cluster (categorical)
    class_probabilities = sample_class_probabilities(data,
                                                     class_probabilities,
                                                     phi_prior,
                                                     cluster_labels,
                                                     cat_count,
                                                     num_clusters,
                                                     num_cols
                                                       
    );
    
    // For each individual sample allocation probability
    for(arma::uword j = 0; j < n; j++){
      
      // sample cluster for each point here
      curr_cluster_probs = categorical_cluster_probabilities(data.row(j),
                                                             data,
                                                             class_probabilities,
                                                             cluster_weights,
                                                             num_clusters,
                                                             num_cols);
      
      // If not fixed (i.e. a known label) update the current label
      // use + 1 as R is counting from 1 not 0
      if(fix_vec(j) == 0){
        cluster_labels(j) = predict_ind(curr_cluster_probs) + 1; 
        // cluster_labels(j) = cluster_predictor(curr_cluster_probs);
      }
    }
    
    // Record outputs of interest
    if (i >= burn && (i - burn) % thinning == 0) {
      record.col((i - burn) / thinning) = cluster_labels;
    }
  }
  
  sim = similarity_mat(record);
  
  return List::create(Named("similarity") = sim,
                      Named("class_record") = record);
}