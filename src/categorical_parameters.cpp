# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;


// update the concentration parameter in the Dirichlet distribution
arma::vec CalcConcentrationn(arma::vec concentration_0,
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
                              arma::uword n_clusters){
  
  // Initialise the cluster weight vector
  arma::vec cluster_weight = arma::zeros<arma::vec>(n_clusters);
  
  // Calculate the concentration parameter
  arma::vec concentration(n_clusters);
  concentration = CalcConcentrationn(concentration_0,
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
arma::vec SampleDirichletPosteriorClass(arma::vec concentration_0,
                                    arma::uvec cluster_labels,
                                    arma::uword n_clusters){
  // Initialise the cluster weight vector
  arma::vec cluster_weight = arma::zeros<arma::vec>(n_clusters);
  
  // Calculate the concentration parameter
  arma::vec concentration(n_clusters);
  concentration = CalcConcentrationnClass(concentration_0,
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


