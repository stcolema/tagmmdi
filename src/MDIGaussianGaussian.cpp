# include <RcppArmadillo.h>
// # include <iostream>
// # include <fstream>
# include "common_functions.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// MDI clustering for two gaussian datasets
// [[Rcpp::export]]
Rcpp::List mdiGaussGauss(arma::mat data_1,
                         arma::mat data_2,
                         arma::vec mu_0_1,
                         double lambda_0_1,
                         arma::mat scale_0_1,
                         int df_0_1,
                         arma::vec mu_0_2,
                         double lambda_0_2,
                         arma::mat scale_0_2,
                         int df_0_2,
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
                         arma::uword thinning,
                         bool outlier_1 = false,
                         double t_df_1 = 4.0,
                         bool outlier_2 = false,
                         double t_df_2 = 4.0,
                         bool normalise_1 = false,
                         bool normalise_2 = false,
                         double u_1 = 2,
                         double v_1 = 10,
                         double u_2 = 2,
                         double v_2 = 10,
                         arma::uword rate_1_0 = 1,
                         arma::uword rate_2_0 = 1
) {
  
  
  
  // As we tend to iterate to less than num_iter, add 1 to itself
  // Done here as many objects are sized based on this.
  num_iter++;
  
  arma::uword n = data_1.n_rows; // Declare the sample size and dimensionality
  arma::uword n_cols_1 = data_1.n_cols;
  arma::uword n_cols_2 = data_2.n_cols;
  arma::uword min_n_clust = std::min(n_clust_1, n_clust_2); // Frequently will compare clusters across contexts so this is a useful bound
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  arma::uword predicted_class_1 = 0;
  arma::uword predicted_class_2 = 0;
  arma::uword record_iter = 0; // Holds the current count of recorded iterations
  arma::uword predicted_outlier_1 = 0;
  arma::uword predicted_outlier_2 = 0;
  arma::uword start_iter = 0;
  arma::uword phi_count = 0;
  double v = 1.0; // strategic latent variable
  double phi = arma::randg(arma::distr_param(a_0, 1/b_0) ); // Context similarity - sample prior
  double Z = 1.0; // Declare normalising constant
  double predicted_norm_likelihood_1 = 0.0; // The likelihood of being non-outlier
  double predicted_norm_likelihood_2 = 0.0;// The likelihood of being non-outlier
  double b_1 = 0.0; // the number of items in the outlier component
  double b_2 = 0.0; // the number of items in the outlier component
  double outlier_likelihood_1 = 0.0;
  double outlier_weight_1 = 1 - sampleBetaDistn(u_1, v_1);
  double outlier_likelihood_2 = 0.0;
  double outlier_weight_2 = 1 - sampleBetaDistn(u_2, v_2);
  arma::uvec outlier_vec_1(n);
  arma::uvec outlier_vec_2(n);
  arma::uvec relevant_labels_1(n); // Class labels of points not currently assigned as outliers
  arma::uvec relevant_labels_2(n);
  arma::vec curr_prob_vec_1(n_clust_1);
  arma::vec curr_prob_vec_2(n_clust_2);
  arma::vec clust_weights_1(n_clust_1);
  arma::vec clust_weights_2(n_clust_2);
  arma::vec global_mean_1(n_cols_1);
  arma::vec global_mean_2(n_cols_2);
  arma::vec labels_weights_phi(n + n_clust_2 + 1); // To hold output of label flipping function
  arma::vec phi_record(eff_count); // To record the context similarity parameter
  arma::vec entropy_cw(num_iter); // record the entropy at each iteration
  arma::vec rate_0_1(n_clust_1); // Rate priors for weight smapling from Gamma distn
  arma::vec rate_0_2(n_clust_2); // Rate priors for weight smapling from Gamma distn
  arma::vec curr_class_probs_1(n_clust_1);
  arma::vec norm_likelihoods_1(n_clust_1);
  arma::vec curr_class_probs_2(n_clust_2);
  arma::vec norm_likelihoods_2(n_clust_2);
  arma::vec curr_outlier_prob_1(2);
  arma::vec curr_outlier_prob_2(2);
  arma::umat outlier_probs_saved_1(n, eff_count); // To save the outlier labels
  arma::umat outlier_probs_saved_2(n, eff_count); // To save the outlier labels
  arma::umat record_1(n, eff_count);
  arma::umat record_2(n, eff_count);
  arma::mat global_variance_1(n_cols_1, n_cols_1); // For outlier assignmnet
  arma::mat global_variance_2(n_cols_2, n_cols_2);
  arma::mat alloc_prob_1(n, n_clust_1); // The allocation probabilities for each class
  arma::mat alloc_prob_2(n, n_clust_2); // The allocation probabilities for each class
  arma::mat alloc_prob_curr_1(n, n_clust_1); // allocation probability matrix
  arma::mat alloc_prob_curr_2(n, n_clust_2); // allocation probability matrix
  arma::mat mu_n_1(n_cols_1, n_clust_1);
  arma::mat mu_n_2(n_cols_2, n_clust_2);
  arma::mat sim_1(n, n); 
  arma::mat sim_2(n, n);
  arma::cube variance_n_1(n_cols_1, n_cols_1, n_clust_1);
  arma::cube variance_n_2(n_cols_2, n_cols_2, n_clust_2);
  // Local 
  
  // Placeholder prior
  rate_0_1.fill(rate_1_0);
  rate_0_2.fill(rate_2_0);
  
  //  Normalise the continuous data
  if(normalise_1){
    data_1 = arma::normalise(data_1);
  }
  
  if(normalise_2){
    data_2 = arma::normalise(data_2);
  }
  
  // Cluster weights for each dataset
  clust_weights_1.zeros();
  clust_weights_2.zeros();
  
  // These are the local cubes of posterior mean and variance overwritten each
  // iteration
  mu_n_1.zeros();
  variance_n_1.zeros();
  
  mu_n_2.zeros();
  variance_n_2.zeros();

  // Record of clustering across MCMC iterations
  record_1.zeros();
  record_2.zeros();
  
  // OUTLIER COMPONENT
  // Variables to handle outlier component from tagm
  // Vector of 0 and 1 for points assigned to outlier group or not
  outlier_vec_1.zeros();
  outlier_vec_2.zeros();
  
  // Begin with all non-fixed points (i.e. unlabelled) to outlier component
  if(outlier_1 && arma::any(fix_vec_1)){
    outlier_vec_1 = 1 - fix_vec_1;
  }
  
  if(outlier_2 && arma::any(fix_vec_2)){
    outlier_vec_2 = 1 - fix_vec_2;
  }
  
  // Declare global variance and mean - used in outlier t-distribution 
  // (if included)
  global_variance_1 = 0.5 * arma::cov(data_1); // Olly's rec
  global_mean_1 = arma::trans(arma::mean(data_1, 0));
  global_variance_2 = 0.5 * arma::cov(data_2); // Olly's rec
  global_mean_2 = arma::trans(arma::mean(data_2, 0));
  
  // Set all entries to zero as we will be using += rather than = (i.e will 
  // never over-write the initial value held, only add to it)
  alloc_prob_1.zeros();
  alloc_prob_2.zeros();
  
  for(arma::uword i = 0; i < num_iter; i++){
    
    // Consider only the labels of points not considered to be outliers
    relevant_labels_1 = clust_labels_1 % (1 - outlier_vec_1);
    relevant_labels_2 = clust_labels_2 % (1 - outlier_vec_2);
    
    // Entropy for graphing convergence
    entropy_cw(i) = calcEntropy(clust_weights_1);
    
    // ## Sample the parameters for the two datasets ##
    
    // Sample the posterior mean and variance for the first dataset
    variance_n_1 = sampleClusterVariance(data_1,
                                         relevant_labels_1,
                                         n_clust_1,
                                         df_0_1,
                                         n_cols_1,
                                         scale_0_1,
                                         lambda_0_1,
                                         mu_0_1);
    
    mu_n_1 = sampleClusterMeans(data_1,
                                relevant_labels_1,
                                n_clust_1,
                                n_cols_1,
                                variance_n_1,
                                lambda_0_1,
                                mu_0_1);
    
    // Sample the posterior mean and variance for the second dataset
    variance_n_2 = sampleClusterVariance(data_2,
                                         clust_labels_2,
                                         n_clust_2,
                                         df_0_2,
                                         n_cols_2,
                                         scale_0_2,
                                         lambda_0_2,
                                         mu_0_2);
    
    mu_n_2 = sampleClusterMeans(data_2,
                                clust_labels_2,
                                n_clust_2,
                                n_cols_2,
                                variance_n_2,
                                lambda_0_2,
                                mu_0_2);
    
    // ## Sample cluster weights for the two datasets ##
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
    
    // Calculate the current normalising constant (consider being more clever
    // about this)
    Z = calcNormalisingConst(clust_weights_1,
                             clust_weights_2,
                             phi,
                             n_clust_1,
                             n_clust_2);
    
    // sample the strategic latent variable, v
    v = arma::randg( arma::distr_param(n, 1.0/Z ) );
    
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
    
    if(outlier_1){
      b_1 = sum(outlier_vec_1);
      outlier_weight_1 = sampleBetaDistn(b_1 + u_1, n + v_1 - b_1);
    }
    
    if(outlier_2){
      b_2 = sum(outlier_vec_2);
      outlier_weight_2 = sampleBetaDistn(b_2 + u_2, n + v_2 - b_2);
    }
    
    // sample class for each observation
    for(arma::uword j = 0; j < n; j++){
      
      // for each point create the vector of probabilities associated with 
      // assignment to each cluster within the gaussian dataset
      // Calculate the log-likelihoods for each cluster in each context
      norm_likelihoods_1 = sampleMDIGaussClustProbs(j,
                                                    data_1,
                                                    n_clust_1,
                                                    mu_n_1,
                                                    variance_n_1,
                                                    phi,
                                                    clust_weights_1,
                                                    relevant_labels_1,
                                                    clust_labels_2);
      
      norm_likelihoods_2 = sampleMDIGaussClustProbs(j,
                                                    data_2,
                                                    n_clust_2,
                                                    mu_n_2,
                                                    variance_n_2,
                                                    phi,
                                                    clust_weights_2,
                                                    relevant_labels_2,
                                                    clust_labels_1);

      // Predict the component to which the obersation belongs
      predicted_class_1 = predictCluster(curr_prob_vec_1, 1);
      predicted_class_2 = predictCluster(curr_prob_vec_2, 1);
      
      // The various probabilities to determine if the observation is considered 
      // an outlier or not (if using TAGM rather than traditional mixtures)
      if(outlier_1){
        
        // Calculate the likelihood associated with the outlier component
        outlier_likelihood_1 =  calcTdistnLikelihood(data_1.row(j).t(), 
                                                     global_mean_1, 
                                                     global_variance_1, 
                                                     n_cols_1,
                                                     t_df_1);
        
        
        outlier_likelihood_1 += log(outlier_weight_1);
        
        curr_outlier_prob_1(1) = outlier_likelihood_1;
        
        // The normal likelihood for the current class allocation
        predicted_norm_likelihood_1 = calcNormalLikelihood(data_1.row(j).t(),
                                                           mu_n_1.col(predicted_class_1 - 1),
                                                           variance_n_1.slice(predicted_class_1 - 1),
                                                           n_cols_1);
        
        predicted_norm_likelihood_1 += log(1 - outlier_weight_1);
        curr_outlier_prob_1(0) = predicted_norm_likelihood_1;
        
        // Predict outlier or not
        predicted_outlier_1 = predictCluster(curr_outlier_prob_1);
        
      }
      
      if(outlier_2){
        
        // Calculate the likelihood associated with the outlier component
        outlier_likelihood_2 =  calcTdistnLikelihood(data_2.row(j).t(), 
                                                     global_mean_2, 
                                                     global_variance_2, 
                                                     n_cols_2,
                                                     t_df_2);
        
        
        outlier_likelihood_2 += log(outlier_weight_2);
        
        curr_outlier_prob_2(1) = outlier_likelihood_2;
        
        // The normal likelihood for the current class allocation
        predicted_norm_likelihood_2 = calcNormalLikelihood(data_2.row(j).t(),
                                                           mu_n_2.col(predicted_class_2 - 1),
                                                           variance_n_2.slice(predicted_class_2 - 1),
                                                           n_cols_2);
        
        predicted_norm_likelihood_2 += log(1 - outlier_weight_2);
        curr_outlier_prob_2(0) = predicted_norm_likelihood_2;
        
        predicted_outlier_2 = predictCluster(curr_outlier_prob_2);
      }
      
      
      if (i >= burn && (i - burn) % thinning == 0) { // include  && record_posteriors in if?
        alloc_prob_curr_1.row(j) = arma::trans(curr_prob_vec_1);
        alloc_prob_curr_2.row(j) = arma::trans(curr_prob_vec_2);
      }
      
      // update labels - in gaussian data this is only if the current point is 
      // not fixed
      if(fix_vec_1[j] == 0){
        clust_labels_1(j) = predicted_class_1;
        if(outlier_1){
          outlier_vec_1(j) = predicted_outlier_1;
        }
      }
      
      if(fix_vec_2[j] == 0){
        clust_labels_2(j) = predicted_class_2;
        if(outlier_2){
          outlier_vec_2(j) = predicted_outlier_2;
        }
      }
    }
    
    // Update cluster labels in second dataset
    // Return the new labels, weights and similarity in a single vector
    // if(data_2_unsupervised){
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
    // }
    
    // Separate the output into the relevant components
    clust_labels_2 = arma::conv_to<arma::uvec>::from(labels_weights_phi.subvec(0, n - 1));
    clust_weights_2 = labels_weights_phi.subvec(n, n + n_clust_2 - 1);
    phi = arma::as_scalar(labels_weights_phi(n + n_clust_2));
    
    // if current iteration is a recorded iteration, save the labels
    if (i >= burn && (i - burn) % thinning == 0) {
      
      // record_iter = ((i - burn) / thinning) + num_load ;
      
      phi_record(phi_count) = phi;
      phi_count++;
      
      // std::cout << "Here.\n";
      record_1.col(record_iter) = clust_labels_1;
      record_2.col(record_iter) = clust_labels_2;
      outlier_probs_saved_1.col(record_iter) = outlier_vec_1;
      outlier_probs_saved_2.col(record_iter) = outlier_vec_2;
      
      alloc_prob_1 += alloc_prob_curr_1;
      alloc_prob_2 += alloc_prob_curr_2;
      
      record_iter++;
    }
  }
  
  // Normalise the allocation probabilities
  alloc_prob_1 = alloc_prob_1 / eff_count;
  alloc_prob_2 = alloc_prob_2 / eff_count;
  
  // construct similarity matrices
  sim_1 = createSimilarityMat(record_1);
  sim_2 = createSimilarityMat(record_2);
  
  return List::create(Named("similarity_1") = sim_1,
                      Named("similarity_2") = sim_2,
                      Named("allocation_mat_1") = alloc_prob_1,
                      Named("allocation_mat_2") = alloc_prob_2,
                      Named("class_record_1") = record_1,
                      Named("class_record_2") = record_2,
                      Named("context_similarity") = phi_record,
                      Named("entropy") = entropy_cw,
                      Named("outlier_1") = outlier_probs_saved_1,
                      Named("outlier_2") = outlier_probs_saved_2); // ,
  // Named("likelihood") = mdi_recorded_likelihood);
}