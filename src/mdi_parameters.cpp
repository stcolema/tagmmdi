# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

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

// The posterior for a gamma distribution
arma::vec gamma_posterior(arma::vec concentration_0,
                          arma::uvec cluster_labels,
                          arma::uword num_clusters,
                          double rate
){
  
  // Initialise the cluster weights
  arma::vec cluster_weight = arma::zeros<arma::vec>(num_clusters);
  
  // Calculate the concentration parameter
  arma::vec concentration(num_clusters);
  concentration = concentration_n(concentration_0,
                                  cluster_labels,
                                  num_clusters);
  
  // Sample the cluster weights from a Gamma distribution
  for (arma::uword i = 0; i < num_clusters; i++) {
    cluster_weight(i) = arma::randg( arma::distr_param(arma::as_scalar(concentration(i)), 1.0) );
  }
  return cluster_weight;
}

// Calculate the rate for the gamma distribution for the class weights for MDI
// This is defined as the sume of the cluster weights (upweighted by the 
// correlation parameter, phi, if the labels match) multiplied by the variable v
double mdi_cluster_rate(double v,
                        arma::uword n_clust,
                        arma::uword cluster_index,
                        arma::vec cluster_weights,
                        double phi){
  // Initialise b, the rate
  double b = 0.0;

  if(cluster_index < n_clust){
    b = v * (arma::sum(cluster_weights) + phi * cluster_weights(cluster_index));
  } else {
    b = v * (arma::sum(cluster_weights));
  }
  
  // // Loop over the number of clusters
  // for(arma::uword i = 0; i < n_clust; i++){
  //   b += cluster_weights(i) * (1 + phi * (cluster_index == i));
  // }
  // 
  // b = b * v;

  return b;
}



// Returns the cluster weights for MDI
arma::vec mdi_cluster_weights(arma::vec shape_0,
                              arma::vec rate_0,
                              double v,
                              arma::uword n_clust,
                              arma::uword n_clust_comp,
                              arma::vec cluster_weights_comp,
                              arma::uvec cluster_labels,
                              arma::uvec cluster_labels_comp,
                              double phi){
  
  // The number of clusters relevant to MDI is the minimum of the number of
  // clusters present in either dataset
  arma::uword n_rel = min(n_clust, n_clust_comp);
  
  // Initialise the cluster weights, the rate and shape
  arma::vec cluster_weight = arma::zeros<arma::vec>(n_rel);
  double b = 0.0;
  double b_n = 0.0;
  arma::vec shape_n(n_rel);
  
  // Calculate the concentration parameter for each cluster
  shape_n = concentration_n(shape_0,
                            cluster_labels,
                            n_rel) + 1;
  
  for (arma::uword i = 0; i < n_rel; i++) {
  
    // Calculate the rate based upon the current clustering
    b = mdi_cluster_rate(v,
                         n_rel,
                         i,
                         cluster_weights_comp,
                         phi);
    
    // Update the prior
    b_n = b + rate_0(i);
    
    // Sample the weights from a gamma distribution
    cluster_weight(i) = arma::randg(arma::distr_param(shape_n(i), 1 / b_n));

  }
  return cluster_weight;
}


// returns the count of the number of points with the same label in both contexts
arma::uword count_common_cluster(arma::uvec cl_1,
                                 arma::uvec cl_2,
                                 arma::uword n){
  // Initialise the output
  arma::uword common_cluster = 0;
  
  // Initialise the vector recording common membership
  arma::uvec common_membership(n);
  
  // Compare the clusterings
  common_membership = cl_1 == cl_2;

  // The output is the count of common entries
  common_cluster = arma::sum(common_membership);
  
  
  // for(arma::uword i = 0; i < n; i++){
  //   if(cl_1(i) == cl_2(i)){
  //     common_cluster++;
  //   }
  // }
  return common_cluster;
}


// Calculate the factorial
int factorial(arma::uword n)
{
  if(n <= 1){
    return 1;
  }
  return n * factorial(n - 1);
  
  // return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}



// Returns the log of the factorial of n
double log_factorial(arma::uword n){
  if(n <= 1){
    return 0;
  }
  return log(n) + log_factorial(n - 1);
}


// Calculates the rate based on the data
double observed_rate(double v,
                     arma::uword min_num_clusters,
                     arma::vec cluster_weights_gaussian,
                     arma::vec cluster_weights_categorical){
  // declare the rate
  double b = 0.0;
  
  // the rate is proportional to the sum of the product of cluster weights 
  // across contexts (i.e. the product of cluster_weight_i in context 1 times
  // cluster_weight_i in context 2)
  for(arma::uword i = 0; i < min_num_clusters; i++){
    b = b + cluster_weights_gaussian(i) * cluster_weights_categorical(i);
  }
  
  // the rate is equal to this sum times the strategic latent variable, v
  b = v * b;
  
  return b;
}

// Calculate the rate for the gamma distribution for the class weights for MDI
// This is defined as the sume of the cluster weights (upweighted by the 
// correlation parameter, phi, if the labels match) multiplied by the variable v
double mdi_phi_rate(double v,
                    arma::uword n_clust,
                    arma::vec cluster_weights_1,
                    arma::vec cluster_weights_2){
  // Initialise b, the rate
  double b = 0.0;
  
  // Find the subset of weights to use in calculating b
  arma::vec sub_weights_1(n_clust);
  arma::vec sub_weights_2(n_clust);
  
  sub_weights_1 = cluster_weights_1(arma::span(0, n_clust - 1) );
  sub_weights_2 = cluster_weights_2(arma::span(0, n_clust - 1) );
  
  // Calculate b
  b = v * (sub_weights_1 % sub_weights_2);
  
  double b1 = 0.0;
  
  // Loop over the number of clusters
  for(arma::uword i = 0; i < n_clust; i++){
    b1 += cluster_weights_1(i) * cluster_weights_2(i);
  }
  
  b1 = b1 * v;
  
  if(b1 != b){
    return "x";
  }
  
  return b;
}



// Calculate the weights for each Gamma distribution ohi may be sampled from
arma::vec phi_weights(arma::uword n,
                      double a_0,
                      double b_n){
  arma::vec phi_weights(n);
  phi_weights.zeros();
  
  for(arma::uword i = 0; i < n; i++){
    // this is the weight of which gamma to sample for the phi
    phi_weights(i) = log_factorial(n)
    - log_factorial(i)
    - log_factorial(n - i)
    + log_factorial(i + a_0 - 1)
    - (i + a_0)*log(b_n);
  }
  return phi_weights;
}

// samples a gamma distn for the current iterations context similarity parameter
// (phi in the original 2012 paper).
double sample_phi(arma::uvec cl_1,
                  arma::uvec cl_2,
                  arma::vec cl_wgts_1,
                  arma::vec cl_wgts_2,
                  double v,
                  arma::uword n,
                  arma::uword min_n_clust,
                  double a_0,
                  double b_0){
  
  // The predicted index of the weighted sum to use
  arma::uword pred_ind = 0;
  
  // calculate the shape of the relevant gamma function
  arma::uword count_same_cluster = 0;
  count_same_cluster = count_common_cluster(cl_1, cl_2, n);
  
  arma::vec prob_vec(count_same_cluster);
  prob_vec.zeros();
  
  // calculate the rate
  double b = 0.0;
  // b = observed_rate(v, min_n_clust, cl_wgts_1, cl_wgts_2) + b_0;
  b = mdi_phi_rate(v, min_n_clust, cl_wgts_1, cl_wgts_2) + b_0;
  
  // Initialise the output, the clustering correlation parameter (phi)
  double phi = 0.0;

  // phi is a weighted sum of gammas; see section 1.5 and 1.51 from: 
  // https://github.com/stcolema/tagmmdi_notes/blob/master/notes/mdi_olly.pdf
  if(count_same_cluster > 0){
    
    // Calculated the weight for each Gamma distribution
    prob_vec = phi_weights(count_same_cluster, a_0, b);
    
    // Predict which to use based on prob_vec
    pred_ind = predict_ind(prob_vec);
    
    // Sample phi
    phi = arma::randg( arma::distr_param(pred_ind + a_0, 1.0/b) );
    
  } else {
    // Sample from the prior
    phi = arma::randg( arma::distr_param(a0, 1.0/b) );
  }
  
  return phi;
}


// Calculates the normalising constant for the posterior
double calculate_normalising_constant(arma::vec cl_wgts_1,
                                      arma::vec cl_wgts_2,
                                      double phi,
                                      arma::uword n_clust_1,
                                      arma::uword n_clust_2){
  double Z = 0.0;
  
  for(arma::uword i = 0; i < n_clust_1; i++){
    Z += arma::sum(cl_wgts_1(i) * cl_wgts_2);
    if(i < n_clust_2){
      Z += cl_wgts_1(i) * cl_wgts_2(i) * phi;
    }
  }
  return Z;
}
