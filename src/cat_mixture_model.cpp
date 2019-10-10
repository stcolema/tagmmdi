# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// uses categorical_parameters.cpp for defining classes and sampling
// uses gaussian_clusters.cpp for PredictIndex

// The actual categorical clustering all wrapped up in one function
// [[Rcpp::export]]
// Rcpp::List
arma::umat  CategoricalClustering(arma::umat data,
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
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning); // Number of samples recorded
  arma::uword cpp_class_start = 0;
  arma::uword r_class_start = 1;
  arma::uvec cat_count(num_cols);
  arma::vec cluster_weights(num_clusters);
  arma::vec curr_cluster_probs(num_clusters);
  arma::umat record(n, eff_count);
  arma::mat sim(n, n);
  arma::field<arma::mat> class_probabilities(num_cols);
  
  // Find the number of categories in each column
  cat_count = CountCatergories(data);
  
  // This is the object that will hold the pobability of assignment to each 
  // class for each column. This function defines the matrices within the field,
  // ensuring they are of the correct dimensionality.
  // This object should be a field of matrices (a matrix for each variable / 
  // column) with the ith matrix having dimeionsality of K (the number of clusters)
  // rows and P_i columns where P_i is the number of categories for the ith 
  // variable
  class_probabilities = DeclareClassProbsField(cat_count,
                                               num_cols,
                                               num_clusters);
  
  // Cluster weights
  cluster_weights.zeros();
  
  // Vector to hold probability of assignment for current individual
  curr_cluster_probs.zeros();

  // The matrix of labels to record
  record.zeros();
  
  // The posterior similarity matrix (I think this should be outside this function)
  sim.zeros()
  
  // Iterate over the desired number of iterations
  for(arma::uword i = 0; i < num_iter; i++){
    
    // Sample cluster weights from a Dirichlet distribution
    cluster_weights = SampleDirichletPosterior(cluster_weight_priors,
                                               cluster_labels,
                                               num_clusters,
                                               r_class_start);
    
    
    // sample probabilities for each class for each cluster (categorical)
    class_probabilities = SampleCategoryProbabilities(data,
                                                      class_probabilities,
                                                      phi_prior,
                                                      cluster_labels,
                                                      cat_count,
                                                      num_clusters,
                                                      num_cols,
                                                      cpp_class_start
                                                       
    );
    
    // For each individual sample allocation probability
    for(arma::uword j = 0; j < n; j++){
      
      // sample cluster for each point here
      curr_cluster_probs = SampleCategoricalDistn(data.row(j),
                                                  data,
                                                  class_probabilities,
                                                  cluster_weights,
                                                  num_clusters,
                                                  num_cols);
      
      // If not fixed (i.e. a known label) update the current label
      // use + 1 as R is counting from 1 not 0
      if(fix_vec(j) == 0){
        cluster_labels(j) = PredictIndex(curr_cluster_probs) + 1; 
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