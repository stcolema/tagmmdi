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

// Calculate the entropy for the current cluster weights
// [[Rcpp::export]]
double entropy(arma::vec class_weights){
  
  // Declare objects
  arma::uword n = class_weights.n_elem;
  arma::vec entropy_components(n);
  double entropy_out = 0.0;
  
  // Calculate the entropy
  entropy_out = - sum(class_weights.t() * log(class_weights));
  
  return entropy_out;
}

// update the concentration parameter in the Dirichlet distribution
arma::vec concentration_n(arma::vec concentration_0,
                          arma::uvec cluster_labels,
                          arma::uword n_cat){
  arma::uvec class_members;
  class_members.zeros();
  arma::uword class_count = 0;
  arma::vec concentration(n_cat);
  concentration.zeros();
  
  // Iterate between 1 and n_cat as the input is based upon R
  for (arma::uword i = 1; i < n_cat + 1; i++) {
    
    // // The number of members in the class is set to 0
    // class_count = 0;
    
    // Find how many labels have the value
    class_members = find(cluster_labels == i);
    class_count = class_members.n_elem;
    
    // Update the concentration parameter using the prior and membership
    concentration(i - 1) = arma::as_scalar(concentration_0(i - 1)) + class_count;
  }
  return concentration;
}

// THIS SEEMS WEIRD
arma::vec concentration_n_class(arma::vec concentration_0,
                                arma::uvec cluster_labels,
                                arma::uword n_cat){
  
  arma::uword class_count = 0;
  arma::vec concentration(n_cat);
  concentration.zeros();
  arma::uvec class_members;
  class_members.zeros();
  
  for (arma::uword i = 0; i < n_cat; i++) {
    class_members = find(cluster_labels == i);
    class_count = class_members.n_elem;
    
    concentration(i) = arma::as_scalar(concentration_0(i)) + class_count;
  }
  return concentration;
}


// sample parameters for a dirichlet distribution (normally for the clusters)
// [[Rcpp::export]]
arma::vec dirichlet_posterior(arma::vec concentration_0,
                              arma::uvec cluster_labels,
                              arma::uword n_clusters){
  
  // Initialise the cluster weight vector
  arma::vec cluster_weight = arma::zeros<arma::vec>(n_clusters);
  
  // Calculate the concentration parameter
  arma::vec concentration(n_clusters);
  concentration = concentration_n(concentration_0,
                                  cluster_labels,
                                  n_clusters);
  
  // Loop over the clusters updating the weight 
  for (arma::uword i = 0; i < n_clusters; i++) {
    
    // Update weights by sampling from a Gamma distribution
    cluster_weight(i) = arma::randg( arma::distr_param(arma::as_scalar(concentration(i)), 1.0) );
    
  }
  
  // Find the sum of the cluster weights (separate object as is a doouble - this is not necessary)
  double total_cluster_weight = sum(cluster_weight);
  
  // Convert the cluster weights (previously gamma distributed) to Beta distributed
  cluster_weight = cluster_weight / total_cluster_weight;
  return cluster_weight;
}

// Dirichlet posterior for class weights (difference of base 0 compared to vanilla
// dirichlet_posterior function
arma::vec dirichlet_posterior_class(arma::vec concentration_0,
                                    arma::uvec cluster_labels,
                                    arma::uword n_clusters){
  // Initialise the cluster weight vector
  arma::vec cluster_weight = arma::zeros<arma::vec>(n_clusters);
  
  // Calculate the concentration parameter
  arma::vec concentration(n_clusters);
  concentration = concentration_n_class(concentration_0,
                                        cluster_labels,
                                        n_clusters);
  
  // Loop over the clusters updating the weight 
  for (arma::uword i = 0; i < n_clusters; i++) {
    
    // Update weights by sampling from a Gamma distribution
    cluster_weight(i) = arma::randg(arma::distr_param(arma::as_scalar(concentration(i)), 1.0) );
  }
  
  // Find the sum of the cluster weights (separate object as is a double - this is not necessary)
  double total_cluster_weight = sum(cluster_weight);
  
  // Convert the cluster weights (previously gamma distributed) to Beta distributed
  cluster_weight = cluster_weight / total_cluster_weight;
  return cluster_weight;
}

// count unique entries in a vector (Armadillo has this already - unnecssary)
arma::uword unique_counter(arma::uvec v){
  std::sort(v.begin(), v.end());
  arma::uword unique_count = std::unique(v.begin(), v.end()) - v.begin();
  return unique_count;
}

// returns a vector of the number of unqiue values in each column
// [[Rcpp::export]]
arma::uvec cat_counter(arma::umat data){
  arma::uword num_cols = data.n_cols;
  arma::uvec num_categories(num_cols);
  for(arma::uword i = 0; i < num_cols; i++){
    num_categories(i) = unique_counter(data.col(i));
  }
  return num_categories;
}

// find the number of categories in each covariate and declare the appropriate
// matrix to record the associated probabilties for each cluster
// [[Rcpp::export]]
arma::field<arma::mat> declare_class_probs_field(arma::uvec cat_count,
                                                 arma::uword num_cols,
                                                 arma::uword num_clusters){
  
  arma::field<arma::mat> class_probabilities(num_cols);
  for(arma::uword i = 0; i < num_cols; i++){
    arma::mat phi_j(num_clusters, cat_count(i));
    phi_j.zeros();
    class_probabilities(i) = phi_j;
  }
  return class_probabilities;
}

// Sample the probabilities for each category across all clusters for each covariate
// [[Rcpp::export]]
arma::field<arma::mat> sample_class_probabilities(arma::umat data,
                                                  arma::field<arma::mat> class_probs,
                                                  arma::field<arma::vec> phi_prior,
                                                  arma::uvec cluster_labels,
                                                  arma::uvec cat_count,
                                                  arma::uword num_clusters,
                                                  arma::uword num_cols
){
  
  arma::umat cluster_data;
  for(arma::uword k = 1; k < num_clusters + 1; k++){
    
    cluster_data = data.rows(find(cluster_labels == k));
    
    for(arma::uword j = 0; j < num_cols; j++){
      
      
      class_probs(j).row(k - 1) = arma::trans(dirichlet_posterior_class(phi_prior(j),
                                        cluster_data.col(j),
                                        cat_count(j)
      )
      );
      
    }
  }
  return class_probs;
}


// Sample the cluster membership of point
// [[Rcpp::export]]
arma::vec categorical_cluster_probabilities(arma::urowvec point,
                                            arma::umat data,
                                            arma::field<arma::mat> class_probabilities,
                                            arma::vec cluster_weights,
                                            arma::uword num_clusters,
                                            arma::uword num_cols){
  
  arma::vec probabilities = arma::zeros<arma::vec>(num_clusters);
  
  double curr_weight = 0.0;
  
  for(arma::uword i = 0; i < num_clusters; i++){
    curr_weight = log(cluster_weights(i));
    for(arma::uword j = 0; j < num_cols; j++){
      
      probabilities(i) = probabilities(i) + std::log(class_probabilities(j)(i, point(j)));
    }
    probabilities(i) = probabilities(i) + curr_weight;
  }
  
  // Overflow handling and convert from logs
  probabilities = exp(probabilities - max(probabilities));
  
  // Normalise
  probabilities = probabilities / sum(probabilities);
  
  return probabilities;
}

// Predicts the cluster assignments based on a vector of probabilities using
// the rejection method
// [[Rcpp::export]]
arma::uword cluster_predictor(arma::vec probabilities){
  double u;
  arma::uword pred;
  u = arma::randu<double>( );
  
  // include + 1 if labels centred on 1
  pred = 1 + sum(u > cumsum(probabilities)); 
  return pred;
}


arma::uword predict_ind(arma::vec my_vec){
  
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
  sim.zeros();
    
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