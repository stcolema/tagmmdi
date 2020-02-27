# include <RcppArmadillo.h>
// # include <iostream>
// # include <fstream>
# include "common_functions.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// This file uses the entropy.cpp, similarity_matrix.cpp,
// the categorical_parameters.cpp (for Dirichlet sampling), 
// gaussian_parameters.cpp (for calculating the parameters for the MVN) 
// and gaussian_clusters.cpp files (for cluster assignment)

// [[Rcpp::export]]
Rcpp::List gaussianClustering(arma::uword num_iter,
                              arma::vec concentration_0,
                              arma::mat scale_0,
                              arma::uvec class_labels,
                              arma::uvec fix_vec,
                              arma::vec mu_0,
                              double lambda_0,
                              arma::mat data,
                              int nu_0,
                              arma::uword k,
                              arma::uword burn,
                              arma::uword thinning,
                              bool outlier = false,
                              double t_df = 4.0,
                              bool normalise = false,
                              double u = 2,
                              double v = 10
) {
  
  // To allow using < and keeping in line with object sizes
  num_iter++;
  
  arma::uword N = data.n_rows;
  arma::uword num_cols = data.n_cols;
  arma::uword eff_count = ceil((double) (num_iter - burn) / (double) thinning);
  arma::uword record_ind = 0;
  arma::uword predicted_outlier = 0;
  arma::uword predicted_class = 0;
  arma::uword cpp_class_start = 0;
  arma::uword r_class_start = 1;
  double b = 0.0;
  double outlier_weight = 0.0;
  double non_outlier_weight = 0.0;
  arma::uvec outlier_vec(N);
  arma::uvec rel_labels(N);
  arma::vec entropy_cw(num_iter);
  arma::vec global_mean(num_cols);
  arma::vec point;
  arma::vec class_weights(k);
  arma::vec cluster_prob(k);
  arma::vec curr_prob_vec(k);
  arma::vec outlier_weights(2);
  arma::vec outlier_prob(2);
  arma::umat record(N, eff_count);
  arma::umat outlier_probs_saved(N, eff_count);
  arma::mat global_variance(num_cols, num_cols);
  arma::mat sim(N, N);
  arma::mat cluster_data;
  arma::mat mu_n(num_cols, k);
  arma::mat alloc_prob_curr(N, k);
  arma::mat alloc_prob(N, k);
  arma::cube variance_n(num_cols, num_cols, k);

  //  Normalise the continuous data
  if(normalise){
    data = arma::normalise(data);
  }
  
  // for use in the outlier distribution
  global_variance = 0.5 * arma::cov(data); 
  global_mean = arma::trans(arma::mean(data, 0));
  
  // Some output matrixces (the PSM and record of cluster assignment)
  record.zeros();
  sim.zeros();
  
  // Cluster specific data
  cluster_data.zeros();
  
  
  // These are the lists recording the posterior mean and 
  // variance for each class for each recorded iteration
  class_weights.zeros();
  outlier_weights.zeros();
  
  // for the calculation of the a, b parameters in the posterior of the outlier 
  // class weight (epsilon in Crook et al 2018) 
  // double b = 0.0; // a will be N - b, no need to declare
  
  // Vector of 0 and 1 for points assigned to outlier group or not
  outlier_vec.zeros();
  
  // Begin with all non-fixed points (i.e. unlabelled) to outlier component
  outlier_vec = 1 - fix_vec;
  

  // Class labels of points not currently assigned as outliers
  rel_labels.zeros();
  
  // the current iterations probabilities of being an outlier (non-outlier prob
  // is 1 - outlier_prob)
  outlier_prob.zeros();
  
  // Record allocation probabilities
  alloc_prob.zeros();
  
  for(arma::uword ii = 0; ii < num_iter; ii++) {
    
    // To see which points are relevant for defining component parameters
    // use pairwise multiplication between the current label and the outlier
    rel_labels = class_labels % (1 - outlier_vec);
    
    // Update class weights
    class_weights = sampleDirichletPosterior(concentration_0,
                                             class_labels, 
                                             k, 
                                             r_class_start);
    
    // Calculate the entropy of the current weights
    entropy_cw(ii) = calcEntropy(class_weights);
    
    // Sample cluster specific parameters
    variance_n = sampleClusterVariance(data,
                                       rel_labels,
                                       k,
                                       nu_0,
                                       num_cols,
                                       scale_0,
                                       lambda_0,
                                       mu_0);
    
    mu_n = sampleClusterMeans(data,
                              rel_labels,
                              k,
                              num_cols,
                              variance_n,
                              lambda_0,
                              mu_0);
    
    // If outliers are allowed calculate the outlier component weight
    if(outlier) {
      b = (double) sum(outlier_vec);
      outlier_weight = sampleBetaDistn(b + u, N + v - b);
      non_outlier_weight = 1.0 - outlier_weight;
    }
    
    
    for (arma::uword jj = 0; jj < N; jj++) {
      // sample class allocation
      point = arma::trans(data.row(jj));
      
      cluster_prob = sampleGaussianMembership(point, 
                                              data,
                                              k, 
                                              class_weights, 
                                              mu_n,
                                              variance_n
      );
      
      // Predict the label (+1 due to R using a 1:n system rather than 0:(n-1))
      predicted_class = predictCluster(curr_prob_vec, 1); 
      
      if(outlier) {
        
        // Calculate the probabilities of being an outlier or not
        outlier_prob = calculateOutlierProb(point, 
                                            global_mean,
                                            global_variance,
                                            num_cols,
                                            t_df,
                                            outlier_weight,
                                            mu_n.col(predicted_class - 1),
                                            variance_n.slice(predicted_class - 1)
                                            );
        
        // Predict membership in the outlier set or not
        predicted_outlier = predictCluster(outlier_prob);
      }
      
      if (ii >= burn && (ii - burn) % thinning == 0) {
        alloc_prob_curr.row(jj) = arma::trans(curr_prob_vec);
      }
      
      // If point is not known, update allocation
      if(fix_vec[jj] == 0) {
        class_labels(jj) = predicted_class;
        outlier_vec(jj) = predicted_outlier;
      }
    }
    
    
    // Record output
    if (ii >= burn && (ii - burn) % thinning == 0) {
      record_ind = (ii - burn) / thinning;
      record.col(record_ind) = class_labels;
      outlier_probs_saved.col(record_ind) = outlier_vec;
      alloc_prob += alloc_prob_curr;
      
    }
  }
  
  // Calculate similarity matrix
  sim = createSimilarityMat(record);
  
  // Convert sum of recorded allocaiton probabilities to average probability
  alloc_prob = alloc_prob / eff_count; 
  
  return List::create(Named("similarity") = sim,
                      Named("class_record") = record,
                      Named("entropy") = entropy_cw,
                      Named("outliers") = outlier_probs_saved);
}
