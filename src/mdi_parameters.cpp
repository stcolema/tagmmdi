# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// Uses gaussian_clusters.cpp for PredictIndex
// uses categorical_parameters.cpp for calculating the concentration

// The posterior for a gamma distribution
// Old name: gamma_posterior
arma::vec SampleGammaPosterior(arma::vec concentration_0,
                               arma::uvec cluster_labels,
                               arma::uword num_clusters,
                               double rate
){
  
  // Initialise the cluster weights and the concentration parameter
  arma::uword r_class_start = 1;
  arma::vec cluster_weight = arma::zeros<arma::vec>(num_clusters);
  arma::vec concentration = arma::zeros<arma::vec>(num_clusters);
  
  // Calculate the concentration parameter
  concentration = CalcConcentrationn(concentration_0,
                                     cluster_labels,
                                     num_clusters,
                                     r_class_start);
  
  // Sample the cluster weights from a Gamma distribution
  for (arma::uword i = 0; i < num_clusters; i++) {
    cluster_weight(i) = arma::randg( arma::distr_param(arma::as_scalar(concentration(i)), 1.0) );
  }
  return cluster_weight;
}

// Calculate the rate for the gamma distribution for the class weights for MDI
// This is defined as the sume of the cluster weights (upweighted by the 
// correlation parameter, phi, if the labels match) multiplied by the variable v
// Old name: mdi_cluster_rate
double CalcRateMDIClassWeights(double v,
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
  return b;
}



// Returns the cluster weights for MDI
// Old name: mdi_cluster_weights
arma::vec SampleMDIClusterWeights(arma::vec shape_0,
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
  // Initialise the cluster weights, the rate and shape
  arma::uword n_rel = min(n_clust, n_clust_comp);
  arma::uword r_class_start = 1;
  double b = 0.0;
  double b_n = 0.0;
  arma::vec cluster_weight = arma::zeros<arma::vec>(n_rel);
  arma::vec shape_n = arma::zeros<arma::vec>(n_rel);
  
  // Calculate the concentration parameter for each cluster
  shape_n = CalcConcentrationn(shape_0,
                            cluster_labels,
                            n_rel,
                            r_class_start) + 1;
  
  for (arma::uword i = 0; i < n_rel; i++) {
  
    // Calculate the rate based upon the current clustering
    b = CalcRateMDIClassWeights(v,
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
// Old name: count_common_cluster
arma::uword CountCommonLabel(arma::uvec cl_1,
                             arma::uvec cl_2,
                             arma::uword n){
  
  arma::uword common_cluster = 0; // output
  arma::uvec common_membership(n); // vector recording common membership
  
  // Compare the clusterings
  common_membership = cl_1 == cl_2;

  // The output is the count of common entries
  common_cluster = arma::sum(common_membership);
  
  return common_cluster;
}


// Calculate the factorial
int Factorial(arma::uword n)
{
  if(n <= 1){
    return 1;
  }
  return n * Factorial(n - 1);
  
  // return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}



// Returns the log of the factorial of n
double LogFactorial(arma::uword n){
  if(n <= 1){
    return 0;
  }
  return log(n) + LogFactorial(n - 1);
}


// Calculates the rate based on the data
// Old name: observed_rate
double CalcObservedRate(double v,
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
// Old name: mdi_phi_rate
double CalcRateMDIPhi(double v,
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
  
  // double b1 = 0.0;
  // 
  // // Loop over the number of clusters
  // for(arma::uword i = 0; i < n_clust; i++){
  //   b1 += cluster_weights_1(i) * cluster_weights_2(i);
  // }
  // 
  // b1 = b1 * v;
  // 
  // if(b1 != b){
  //   return "x";
  // }
  // 
  
  return b;
}



// Calculate the weights for each Gamma distribution ohi may be sampled from
// Old name: phi_weights
arma::vec CalcMDIPhiWeights(arma::uword n,
                      double a_0,
                      double b_n){
  arma::vec phi_weights(n);
  phi_weights.zeros();
  
  for(arma::uword i = 0; i < n; i++){
    // this is the weight of which gamma to sample for the phi
    phi_weights(i) = LogFactorial(n)
    - LogFactorial(i)
    - LogFactorial(n - i)
    + LogFactorial(i + a_0 - 1)
    - (i + a_0)*log(b_n);
  }
  return phi_weights;
}

// samples a gamma distn for the current iterations context similarity parameter
// (phi in the original 2012 paper).
// Old name: sample_phi
double SampleMDIPhi(arma::uvec cl_1,
                    arma::uvec cl_2,
                    arma::vec cl_wgts_1,
                    arma::vec cl_wgts_2,
                    double v,
                    arma::uword n,
                    arma::uword min_n_clust,
                    double a_0,
                    double b_0){
  
  // The predicted index of the weighted sum to use
  arma::uword count_same_cluster = 0;
  arma::uword pred_ind = 0;
  double b = 0.0; // the rate
  double phi = 0.0; //output, the clustering correlation parameter (phi)
  arma::vec prob_vec(count_same_cluster);
  
  // calculate the shape of the relevant gamma function (this is the count of
  // the points with a common label across contexts)
  count_same_cluster = CountCommonLabel(cl_1, cl_2, n);
  prob_vec.zeros();
  
  // calculate the rate
  b = CalcRateMDIPhi(v, min_n_clust, cl_wgts_1, cl_wgts_2) + b_0;
  
  // phi is a weighted sum of gammas; see section 1.5 and 1.51 from: 
  // https://github.com/stcolema/tagmmdi_notes/blob/master/notes/mdi_olly.pdf
  if(count_same_cluster > 0){
    
    // Calculated the weight for each Gamma distribution
    prob_vec = CalcMDIPhiWeights(count_same_cluster, a_0, b);
    
    // Predict which to use based on prob_vec
    pred_ind = PredictIndex(prob_vec);
    
    // Sample phi
    phi = arma::randg( arma::distr_param(pred_ind + a_0, 1.0/b) );
    
  } else {
    // Sample from the prior
    phi = arma::randg( arma::distr_param(a_0, 1.0/b) );
  }
  
  return phi;
}


// Calculates the normalising constant for the posterior
// Old name: calculate_normalising_constant
double CalcNormalisingConst(arma::vec cl_wgts_1,
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


// Sample the cluster membership of a categorical sample for MDI
// Old name: mdi_cat_clust_prob
arma::vec SampleMDICatClustProb(arma::uword row_index,
                             arma::umat data,
                             arma::field<arma::mat> class_probs,
                             arma::uword num_clusters,
                             arma::uword num_cols_cat,
                             double phi,
                             arma::vec cluster_weights_categorical,
                             arma::uvec clust_labels,
                             arma::uvec clust_labels_comp){
  
  // cluster_labels_comparison is the labels of the data in the other context
  arma::uword common_cluster = 0;
  double curr_weight = 0.0;
  double similarity_upweight = 0.0; // Upweight for similarity of contexts
  arma::urowvec point = data.row(row_index);
  arma::vec prob_vec = arma::zeros<arma::vec>(num_clusters);
  
  for(arma::uword i = 0; i < num_clusters; i++){
    
    // calculate the log-weights for the context specific cluster and the across
    // context similarity
    // pretty much this is the product of probabilities possibly up-weighted by
    // being in the same cluster in a different context and weighted by the cluster
    // weight in the current context
    curr_weight = log(cluster_weights_categorical(i));
    
    // Check if in the same cluster in both contexts
    common_cluster = 1 * (clust_labels_comp(row_index) == clust_labels(row_index));
    
    similarity_upweight = log(1 + phi * common_cluster);
    
    for(arma::uword j = 0; j < num_cols_cat; j++){
      
      prob_vec(i) = prob_vec(i) + std::log(class_probs(j)(i, point(j)));
    }
    
    // As logs can sum rather than multiply the components
    prob_vec(i) = curr_weight + prob_vec(i) + similarity_upweight;
  }
  
  // // to handle overflowing
  // prob_vec = exp(prob_vec - max(prob_vec));
  // 
  // // normalise
  // prob_vec = prob_vec / sum(prob_vec);
  
  return prob_vec;
}

// sample calculate probabilties for cluster allocation in the gaussian data in 
// MDI (hence the presence of the context similarity parameter)
// Old name: mdi_gauss_clust_probs
arma::vec SampleMDIGaussClustProbs(arma::uword row_index,
                                arma::mat data,
                                arma::uword k,
                                arma::mat mu,
                                arma::cube variance,
                                double context_similarity,
                                arma::vec cluster_weights,
                                arma::uvec cluster_labels,
                                arma::uvec cluster_labels_comp){
  
  arma::uword d = data.n_cols;
  arma::uword common_cluster = 0;
  double curr_weight = 0.0;
  double log_likelihood = 0.0;
  double similarity_upweight = 0.0; // Upweight for similarity of contexts
  arma::uvec count_probs;
  arma::vec prob_vec(k);
  arma::vec point = arma::trans(data.row(row_index));
  
  common_cluster = 1 * (cluster_labels(row_index) == cluster_labels_comp(row_index));
  
  for(arma::uword i = 0; i < k ; i++){
    curr_weight = log(cluster_weights(i));
    
    // If in the cluster that the point is in in the comparison context, upweight
    if((i + 1) == cluster_labels(row_index)){
      similarity_upweight = log(1 + context_similarity * common_cluster);
    }
    
    log_likelihood = CalcNormalLikelihood(point, mu.col(i), variance.slice(i), d);
    
    prob_vec(i) = curr_weight + log_likelihood + similarity_upweight;
    
    // This is reset to 0 for the next iteration
    similarity_upweight = 0.0;
    
  }
  
  return prob_vec;
}



// In a vector changes all values of ``label_1'' to ``label_2'' and vice versa
// Old name: swap_labels
arma::uvec SwapLabels(arma::uvec cluster_labels, 
                      arma::uword label_1, 
                      arma::uword label_2){
  
  arma::uvec label_1_ind = find(cluster_labels == label_1);
  arma::uvec label_2_ind = find(cluster_labels == label_2);
  
  cluster_labels.elem(label_1_ind).fill(label_2);
  cluster_labels.elem(label_2_ind).fill(label_1);
  return cluster_labels;
}

// Swap cluster weights
// Old name: swap_cluster_weights
arma::vec SwapClusterWeights(arma::vec cluster_weights,
                               arma::uword label_1, 
                               arma::uword label_2){
  
  arma::vec decoy_weights = cluster_weights;
  
  cluster_weights(label_1) = decoy_weights(label_2);
  cluster_weights(label_2) = decoy_weights(label_1);
  
  return cluster_weights;
}


// This calculates a score for a labelling. This is to allow comparison of the
// current labelling to an updated labelling (i.e. to enable label flipping).
// Old name: comp
double CalcLabellingScore(arma::uword n,
                          double context_similarity,
                          arma::uvec cluster_labels_1,
                          arma::uvec cluster_labels_2){
  
  double score = 0.0;
  
  for(arma::uword i = 0; i < n; i++){
    score += log(1 + context_similarity * (cluster_labels_1(i) == cluster_labels_2(i)));
  }
  
  return score;
}


//  Have to create a function for label swapping
// This will involve comparing likelihoods with and without swap and then 
// a rejection method

// will have to re-order cluster weights vector for dataset 2; easiest to record
// which classes flip and hence the index of the gammas which have to swap
// to generate a vector of the indices use: 
// std::vector<int> v(100) ; // vector with 100 ints.
// std::iota (std::begin(v), std::end(v), 0);
// 
// Will compare improvement in context similarity if cat_cluster_label(j) changed
// to associating with cont_cluster_label(i) and then record swapping of (j) and 
// (i) in the cluster labels in context 2 and the change in the gamma vector
// Old name: cluster_label_update
arma::vec UpdateClusterLabels(arma::uvec cluster_labels_1,
                               arma::uvec cluster_labels_2,
                               arma::vec cluster_weights_1,
                               arma::vec cluster_weights_2,
                               arma::uword num_clusters_1,
                               arma::uword num_clusters_2,
                               double phi,
                               arma::uword min_num_clusters,
                               double v,
                               arma::uword n,
                               double a0,
                               double b0,
                               double Z){
  
  arma::uword new_pos = 0;
  double log_accept = 0.0;
  double accept = 0.0;
  double old_score = 0.0;
  double new_score = 0.0;
  arma::uvec new_labels(n);
  arma::vec new_weights(num_clusters_2);

  old_score = CalcLabellingScore(n, phi, cluster_labels_1, cluster_labels_2);
  
  
  // Should this be bound as the min_num_clusters or min_num_clusters - 1?
  for(arma::uword i = 0; i < num_clusters_2; i++){
    
    // Generate a new position not equal to i
    // multiply a random number from [0,1] by the upper bound of i less 1
    // this less 1 is used so that if we equal i we can add 1 and treat all
    // cases similarly (otherwise we would be biasing the result towards i + 1)
    // Should it be "- 2" rather than "- 1"? Due to C++'s method of accessing
    // elements and the +1 to avoid selecting the same position as i, I think we
    // need to use "- 2".
    new_pos = floor(arma::randu<double>( ) * (num_clusters_2 - 1));
    
    if(new_pos >= i){
      new_pos++;
    }
    
    new_labels = SwapLabels(cluster_labels_2, i + 1, new_pos + 1);
    
    new_weights = SwapClusterWeights(cluster_weights_2, i, new_pos);
    
    new_score = CalcLabellingScore(n, phi, cluster_labels_1, new_labels);
    
    log_accept = new_score - old_score;
    
    accept = 1.0;
    
    if(log_accept < 0){
      accept = exp(log_accept);
    }
    
    if(arma::randu<double>( ) < accept){
      cluster_labels_2 = new_labels;
      cluster_weights_2 = new_weights;
      
      old_score = CalcLabellingScore(n, phi, cluster_labels_1, cluster_labels_2);
      
    }
  }
  
  arma::vec output = arma::zeros<arma::vec>(n + num_clusters_2 + 1);
  output.subvec(0, n - 1) = arma::conv_to<arma::vec>::from(cluster_labels_2);
  
  output.subvec(n, n + num_clusters_2 - 1) = cluster_weights_2;
  output(n + num_clusters_2) = phi;
  return output;
}

// Predicts the cluster assignments based on a vector of probabilities using
// the rejection method
// Old name: cluster_predictor
// [[Rcpp::export]]
arma::uword PredictClusterMembership(arma::vec probabilities){
  double u;
  arma::uword pred;
  u = arma::randu<double>( );
  
  // include + 1 if labels centred on 1
  pred = 1 + sum(u > cumsum(probabilities)); 
  return pred;
}