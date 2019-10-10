# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;


// update the concentration parameter in the Dirichlet distribution
arma::vec CalcConcentrationn(arma::vec concentration_0,
                             arma::uvec cluster_labels,
                             arma::uword n_cat,
                             arma::uword count_from){

  arma::uword class_count = 0;
  arma::uvec class_members;
  arma::vec concentration(n_cat);
  class_members.zeros();
  concentration.zeros();
  
  // Iterate between 1 and n_cat as the input is based upon R
  for (arma::uword i = count_from; i < n_cat + count_from; i++) {
    
    // Find how many labels have the value
    class_members = find(cluster_labels == i);
    class_count = class_members.n_elem;
    
    // Update the concentration parameter using the prior and membership
    concentration(i - 1) = arma::as_scalar(concentration_0(i - 1)) + class_count;
  }
  return concentration;
}

// THIS SEEMS WEIRD
arma::vec CalcConcentrationnClass(arma::vec concentration_0,
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
arma::vec SampleDirichletPosterior(arma::vec concentration_0,
                              arma::uvec cluster_labels,
                              arma::uword n_clusters,
                              arma::uword count_from){
  
  // Initialise the cluster weight vector
  arma::vec cluster_weight = arma::zeros<arma::vec>(n_clusters);
  
  // Calculate the concentration parameter
  arma::vec concentration(n_clusters);
  concentration = CalcConcentrationn(concentration_0,
                                     cluster_labels,
                                     n_clusters,
                                     count_from);
  
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
arma::vec SampleDirichletPosteriorClass(arma::vec concentration_0,
                                        arma::uvec cluster_labels,
                                        arma::uword n_clusters){
  // Initialise the cluster weight vector
  arma::vec cluster_weight = arma::zeros<arma::vec>(n_clusters);
  
  // Calculate the concentration parameter
  arma::vec concentration(n_clusters);
  concentration = CalcConcentrationn(concentration_0,
                                     cluster_labels,
                                     n_clusters,
                                     1);
  
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
arma::uword CountUniqueCases(arma::uvec v){
  std::sort(v.begin(), v.end());
  arma::uword unique_count = std::unique(v.begin(), v.end()) - v.begin();
  return unique_count;
}

// returns a vector of the number of unqiue values in each column
// [[Rcpp::export]]
arma::uvec CountCatergories(arma::umat data){
  arma::uword num_cols = data.n_cols;
  arma::uvec num_categories(num_cols);
  for(arma::uword i = 0; i < num_cols; i++){
    num_categories(i) = CountUniqueCases(data.col(i));
  }
  return num_categories;
}

// find the number of categories in each covariate and declare the appropriate
// matrix to record the associated probabilties for each cluster
// [[Rcpp::export]]
arma::field<arma::mat> DeclareClassProbsField(arma::uvec cat_count,
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
arma::field<arma::mat> SampleCategoryProbabilities(arma::umat data,
                                                   arma::field<arma::mat> class_probs,
                                                   arma::field<arma::vec> phi_prior,
                                                   arma::uvec cluster_labels,
                                                   arma::uvec cat_count,
                                                   arma::uword num_clusters,
                                                   arma::uword num_cols,
                                                   arma::uword class_start
){
  
  arma::umat cluster_data;
  for(arma::uword k = 1; k < num_clusters + 1; k++){
    
    cluster_data = data.rows(find(cluster_labels == k));
    
    for(arma::uword j = 0; j < num_cols; j++){
      
      
      class_probs(j).row(k - 1) = SampleDirichletPosterior(phi_prior(j),
                                                           cluster_data.col(j),
                                                           cat_count(j),
                                                           class_start
      ).t();

      
    }
  }
  return class_probs;
}

// Sample the cluster membership of point
// Old name: categorical_cluster_probabilities
// [[Rcpp::export]]
arma::vec SampleCategoricalDistn(arma::urowvec point,
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
