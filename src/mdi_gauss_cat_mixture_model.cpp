# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// Refers to categorical_parameters.cpp, mdi_parameters.cpp, entropy.cpp, t_dstn_clusters.cpp

// MDI clustering for a gaussian and cateogrical dataset
// [[Rcpp::export]]
Rcpp::List mdi_gauss_cat(arma::mat cont_data,
                         arma::umat cat_data,
                         arma::vec mu_0,
                         double lambda_0,
                         arma::mat scale_0,
                         int nu_0,
                         double a_0,
                         double b_0,
                         arma::vec clust_wgt_priors_gauss,
                         arma::vec clust_wgt_priors_cat,
                         arma::field<arma::vec> phi_prior,
                         arma::uvec clust_labels_gauss,
                         arma::uvec clust_labels_cat,
                         arma::uword n_clust_gauss,
                         arma::uword n_clust_cat,
                         arma::uvec fix_vec_1,
                         arma::uvec fix_vec_2,
                         arma::uword n_iter,
                         arma::uword burn,
                         arma::uword thin,
                         bool allow_outliers = false,
                         double t_df = 4.0,
                         bool normalise = false,
                         double u_1 = 2,
                         double v_1 = 10,
                         arma::uword rate_gauss_0 = 1,
                         arma::uword rate_cat_0 = 1
){
  
  // As we tend to iterate to less than n_iter, add 1 to itself
  // Done here as many objects are sized based on this.
  n_iter++;
  
  // Declare the sample size and dimensionality of the continuous and 
  // categorical data
  arma::uword eff_count = ceil((double)(n_iter + num_load - burn) / (double)thin);
  arma::uword n = cont_data.n_rows;
  arma::uword n_cols_cont = cont_data.n_cols;
  arma::uword n_cols_cat = cat_data.n_cols;
  arma::uword min_n_clust = std::min(n_clust_gauss, n_clust_cat);
  arma::uword record_iter = 0;
  arma::uword predicted_class = 0; // the predicted class (assuming the point is not an outlier)
  arma::uword predicted_outlier = 0;
  double v = 1.0; // strategic latent variable
  double b = 0.0; // a will be N - b, no need to declare
  double context_similarity = 0.0;
  double Z = 1.0; // Declare the normalising constant
  double curr_outlier_likelihood = 0.0;
  double outlier_weight = 1.0 - sample_beta(u_1, v_1);
  double predicted_norm_likelihood = 0.0; // likelihood of being non-outlier
  double gaussian_score = 0.0; // Likelihood of the current gaussian model
  double cat_score = 0.0; // Likelihood of the current categorical model
  double log_gamma_n = log_factorial(n - 1);
  arma::uvec outlier_vec(n); // binary vector denoting outlier assignment
  arma::uvec cat_count(n_cols_cat);
  arma::uvec relevant_labels(n); // Class labels of non-outlier points
  arma::vec global_mean(n_cols_cont);
  arma::vec clust_wgt_gauss(n_clust_gauss);
  arma::vec clust_wgt_cat(n_clust_cat);
  arma::vec curr_class_probs(n_clust_gauss);
  arma::vec curr_norm_likelihoods(n_clust_gauss);
  arma::vec curr_gauss_prob_vec(n_clust_gauss);
  arma::vec curr_cat_prob_vec(n_clust_cat);
  arma::vec context_similarity_record(eff_count);
  arma::vec labels_weights_phi(n + n_clust_cat + 1); // hold output of label flipping function
  arma::vec entropy_cw(n_iter - num_load); // record of entropy at each iteration
  arma::vec rate_0_gauss(n_clust_gauss); // Rate priors for weight sampling from Gamma distn
  arma::vec rate_0_cat(n_clust_cat);// Rate priors for weight sampling from Gamma distn
  arma::vec curr_outlier_prob(2);
  arma::umat gaussian_record(n, eff_count);
  arma::umat categorical_record(n, eff_count);
  arma::umat outlier_probs_saved(n, eff_count); // to save the outlier labels
  arma::mat mu_n(n_cols_cont, n_clust_gauss);
  arma::mat global_variance(n_cols_cont, n_cols_cont);
  arma::mat alloc_prob_gauss(n, n_clust_gauss); // hold the allocation probabilities
  arma::mat alloc_prob_cat(n, n_clust_cat); // hold the allocation probabilities
  arma::mat alloc_prob_gauss_curr(n, n_clust_gauss); // current allocation probs
  arma::mat alloc_prob_cat_curr(n, n_clust_cat); // current allocation probs
  arma::cube variance_n(n_cols_cont, n_cols_cont, n_clust_gauss);
  arma::field<arma::mat> class_probabilities(n_cols_cat);
  // arma::field<arma::cube> loc_mu_variance(2); // cubes of posterior mean and variance 
  

  
  // arma::cube gaussian_class_probs(eff_count, n_clust_gauss, n); // for recording the probabilities for each class
  // arma::cube cat_class_probs(eff_count, n_clust_cat, n); // for recording the probabilities for each class
  
  
  // Frequently will compare clusters across contexts so this is a useful bound
  // to iterate to
  
  
  //  Normalise the continuous data
  if(normalise){
    cont_data = arma::normalise(cont_data);
  }
  
  // Declare global variance and mean - used in outlier t-distribution 
  // (if included)
  global_variance = 0.5 * arma::cov(cont_data); // Olly's rec
  global_mean = arma::trans(arma::mean(cont_data, 0));
  
  // Cluster weights for each dataset
  clust_wgt_gauss.zeros();
  clust_wgt_cat.zeros();
  
  // Initialise the matrices within this field, defining their dimension 
  // Required for loading
  // for(arma::uword i = 0; i < n_clust_gauss; i++){
  //   loc_mu_variance(0) = arma::cube(n_cols_cont, n_cols_cont, n_clust_gauss);
  //   loc_mu_variance(0).zeros();
  //   loc_mu_variance(1) = arma::cube(n_cols_cont, 1, n_clust_gauss);
  //   loc_mu_variance(1).zeros();
  // }
  
  // Posterior mean and variance
  mu_n.zeros();
  variance_n.zeros();

  // Declare the field for the phi variable for the categorical data
  cat_count = cat_counter(cat_data);
  class_probabilities = DeclareClassProbsField(cat_count,
                                               n_cols_cat,
                                               n_clust_cat);
  
  // Initialise the context similarity parameter based on priors
  context_similarity = arma::randg(arma::distr_param(a_0, 1/b_0) );

  // These are the lists recording the posterior mean and 
  // variance for each class for each recorded iteration

  
  arma::cube cluster_var_record(n_cols_cont, n_cols_cont, eff_count);
  for(arma::uword i = 0; i < n_clust_gauss; i++){
    cluster_var_record.zeros();
    variance(i) = cluster_var_record;
  }

  // Records of allocation
  gaussian_record.zeros();
  categorical_record.zeros();
  
  rate_0_gauss.fill(rate_gauss_0);
  rate_0_cat.fill(rate_cat_0);

  // A positive integer to hold the current count of recorded iterations
  
  
  // OUTLIER COMPONENT
  // Variables to handle outlier component from tagm
  // Vector of 0 and 1 for points assigned to outlier group or not
  
  
  // Begin with all non-fixed points (i.e. unlabelled) to outlier component
  outlier_vec.zeros();
  if(allow_outliers && arma::any(fix_vec_1)){
    outlier_vec = 1 - fix_vec_1;
  }
  
  // the current iterations probabilities (overwritten - saved to the above cube
  // after burn in and as defined by thin)

  
  
  
  
  // the current iterations probabilities of being an outlier (non-outlier prob
  // is 1 - curr_outlier_prob)
  
  // For holding the log-likelihood of each recorded iteration
  // arma::vec mdi_recorded_likelihood(eff_count);
  // mdi_recorded_likelihood.zeros();
  
  // // The total likelihood of the gaussian model in the current iteration
  // double gaussian_score = 0.0;
  // 
  // // The total likelihood of the categorical model in the current iteration
  // double cat_score = 0.0;
  // 
  // double mdi_likelihood = 0.0;
  // double log_gamma_n = log_factorial(n - 1);
  

  
  // Set all entries to zero as we will be using += rather than = (i.e will 
  // never over-write the initial value held, only add to it)
  alloc_prob_gauss.zeros();
  alloc_prob_cat.zeros();
  
  // Local allocation probability matrix

  arma::uword start_iter = 0;
  arma::uword phi_count = 0;
  
  for(arma::uword i = 0; i < n_iter; i++){
    
    // Consider only the labels of points not considered to be outliers
    relevant_labels = clust_labels_gauss % (1 - outlier_vec);
    
    // Entropy for graphing convergence
    entropy_cw(i) = CalcEntropy(clust_wgt_gauss);
    
    // ## Sample the parameters for the two datasets ##
    
    // Sample the posterior mean and variance for the gaussian data
    // loc_mu_variance = mean_variance_sampling(cont_data,
    //                                          relevant_labels,
    //                                          n_clust_gauss,
    //                                          nu_0,
    //                                          n_cols_cont,
    //                                          scale_0,
    //                                          lambda_0,
    //                                          mu_0);
    
    variance_n = SampleVariancePosterior(cont_data,
                                   relevant_labels,
                                   n_clust_gauss,
                                   nu_0,
                                   n_cols_cont,
                                   scale_0,
                                   lambda_0,
                                   mu_0);
    
    mu_n = SampleMeanPosterior(cont_data,
                         relevant_labels,
                         n_clust_gauss,
                         n_cols_cont,
                         variance_n,
                         lambda_0,
                         mu_0);
    
    // For the categorical data, sample the probabilities for each class
    class_probabilities = SampleCategoryProbabilities(cat_data,
                                                     class_probabilities,
                                                     phi_prior,
                                                     clust_labels_cat,
                                                     cat_count,
                                                     n_clust_cat,
                                                     n_cols_cat);
    
    // ## Sample cluster weights for the two datasets ##
    
    // Gaussian weights
    clust_wgt_gauss = SampleMDIClusterWeights(clust_wgt_priors_gauss,
                                                   rate_0_gauss,
                                                   v,
                                                   n_clust_gauss,
                                                   n_clust_cat,
                                                   clust_wgt_cat,
                                                   clust_labels_gauss,
                                                   clust_labels_cat,
                                                   context_similarity);
    
    // Categorical weights
    clust_wgt_cat = SampleMDIClusterWeights(clust_wgt_priors_cat,
                                                      rate_0_cat,
                                                      v,
                                                      n_clust_cat,
                                                      n_clust_gauss,
                                                      clust_wgt_gauss,
                                                      clust_labels_cat,
                                                      clust_labels_gauss,
                                                      context_similarity);
    
    // Calculate the current normalising constant (consider being more clever
    // about this)
    Z = CalcNormalisingConst(clust_wgt_gauss,
                                       clust_wgt_cat,
                                       context_similarity,
                                       n_clust_gauss,
                                       n_clust_cat);
    
    // sample the strategic latent variable, v
    v = arma::randg( arma::distr_param(n, 1.0/Z ) );
    
    // sample the context similarity parameter (as only two contexts this is a
    // single double - number not a drink)
    context_similarity = SampleMDIPhi(clust_labels_gauss,
                                    clust_labels_cat,
                                    clust_wgt_gauss,
                                    clust_wgt_cat,
                                    v,
                                    n,
                                    min_n_clust,
                                    a_0,
                                    b_0);
    
    // Initialise the count of outputs
    // record_iter = num_load - 1;
    
    if(allow_outliers){
      
      // Components of outlier weight
      b = sum(outlier_vec);
      outlier_weight = SampleBetaDistn(b + u_1, n + v_1 - b);
    }
    
    // sample class for each observation
    for(arma::uword j = 0; j < n; j++){
      
      // for each point create the vector of probabilities associated with 
      // assignment to each cluster within the gaussian dataset
      curr_norm_likelihoods = SampleMDIGaussClustProbs(j,
                                                    cont_data,
                                                    n_clust_gauss,
                                                    mu_n,
                                                    variance_n,
                                                    // loc_mu_variance(1),
                                                    // loc_mu_variance(0),
                                                    context_similarity,
                                                    clust_wgt_gauss,
                                                    relevant_labels,
                                                    clust_labels_cat);
      
      // Normalise this vector
      // Predict membership
      predicted_class =  PredictIndex(curr_norm_likelihoods) + 1;
     
     
     
      // curr_gauss_prob_vec = over_flow_handling(curr_norm_likelihoods);
      
      // Using the rejection method sample a single class from this vector
      // predicted_class = cluster_predictor(curr_gauss_prob_vec);
      
      
      // The various probabilities to determine if the observation is considered 
      // an outlier or not (if using TAGM rather than traditional mixtures)
      if(allow_outliers){
        
        // The likelihood associated with the outlier t-distribution
        curr_outlier_likelihood = CalcTdistnLikelihood((cont_data.row(j)).t(), 
                                               global_mean, 
                                               global_variance, 
                                               n_cols_cont,
                                               t_df);
        
        curr_outlier_likelihood += log(outlier_weight);
        
        curr_outlier_prob(1) = curr_outlier_likelihood;
        
        // The normal likelihood for the current class allocation
        predicted_norm_likelihood = normal_likelihood(arma::trans(cont_data.row(j)),
                                                      mu_n.col(predicted_class - 1),
                                                      variance_n.slice(predicted_class - 1),
                                                      n_cols_cont);
        
        predicted_norm_likelihood += log(1 - outlier_weight);
        curr_outlier_prob(0) = predicted_norm_likelihood;
        
        // Predict if point is considered an outlier in current cluster
        predicted_outlier =  PredictIndex(curr_outlier_prob);
        
        // // Overflow handling
        // curr_outlier_prob = exp(curr_outlier_prob - max(curr_outlier_prob));
        // 
        // // Normalise the vector
        // curr_outlier_prob = curr_outlier_prob / sum(curr_outlier_prob);
        // 
        // // Check if we predict the current individual to be assigned as an 
        // // outlier using the rejection method
        // predicted_outlier = cluster_predictor(curr_outlier_prob) - 1;
        // 
        
        gaussian_score += predicted_norm_likelihood * predicted_outlier 
          + curr_outlier_likelihood * (1 - predicted_outlier);
        // gaussian_score = gaussian_score / sum(clust_wgt_gauss)
      }
      
      curr_cat_prob_vec = SampleMDICatClustProb(j,
                                                     cat_data,
                                                     class_probabilities,
                                                     n_clust_cat,
                                                     n_cols_cat,
                                                     context_similarity,
                                                     clust_wgt_cat,
                                                     clust_labels_gauss,
                                                     clust_labels_cat);
      
      
      
      if (i >= burn && (i - burn) % thin == 0) { // include  && record_posteriors in if?
        // record_iter = (i - burn) / thin;
        // gaussian_class_probs.slice(j).row(record_iter) = arma::trans(curr_gauss_prob_vec);
        // cat_class_probs.slice(j).row(record_iter) = arma::trans(curr_cat_prob_vec);
        
        // alloc_prob_gauss.row(j) += arma::trans(curr_gauss_prob_vec);
        // alloc_prob_cat.row(j) += arma::trans(curr_cat_prob_vec);
        
        alloc_prob_gauss_curr.row(j) = arma::trans(curr_gauss_prob_vec);
        alloc_prob_cat_curr.row(j) = arma::trans(curr_cat_prob_vec);
      }
      
      // update labels - in gaussian data this is only if the current point is 
      // not fixed
      if(fix_vec_1[j] == 0){
        clust_labels_gauss(j) = predicted_class; // cluster_predictor(curr_gauss_prob_vec);
        if(allow_outliers){
          outlier_vec(j) = predicted_outlier;
        }
      }
      
      if (fix_vec_2[j] == 0){
        clust_labels_cat(j) = PredictIndex(curr_cat_prob_vec) + 1;
        // clust_labels_cat(j) = cluster_predictor(curr_cat_prob_vec);
      }
    }
    
    
    
    // Update cluster labels in second dataset
    // Return the new labels, weights and similarity in a single vector
    labels_weights_phi = UpdateClusterLabels(clust_labels_gauss,
                                              clust_labels_cat,
                                              clust_wgt_gauss,
                                              clust_wgt_cat,
                                              n_clust_gauss,
                                              n_clust_cat,
                                              context_similarity,
                                              min_n_clust,
                                              v,
                                              n,
                                              a_0,
                                              b_0,
                                              Z);
    
    // Separate the output into the relevant components
    clust_labels_cat = arma::conv_to<arma::uvec>::from(labels_weights_phi.subvec(0, n - 1));
    clust_wgt_cat = labels_weights_phi.subvec(n, n + n_clust_cat - 1);
    context_similarity = arma::as_scalar(labels_weights_phi(n + n_clust_cat));
    
    
    
    // if current iteration is a recorded iteration, save the labels
    if (i >= burn && (i - burn) % thin == 0) {
      
      // record_iter = ((i - burn) / thin) + num_load ;
      
      context_similarity_record(phi_count) = context_similarity;
      phi_count++;
      
      // Record the current iteration's score (note this does not account for
      // label flipping)
      // mdi_recorded_likelihood(record_iter) = mdi_likelihood;
      
      // if(save_results){
      //   // Save to file
      //   
      //   // Make sure to avoid overwriting the loaded entries
      //   // record_iter += num_load - 1;
      //   
      //   // Create appending label of the current iteration
      //   i_str = std::to_string(record_iter);
      //   
      //   // Create the file names for saving the current iteration's information
      //   // The discrete componenet allocations
      //   gauss_lab_file_loc = gauss_lab_file + i_str;
      //   cat_lab_file_loc = cat_lab_file + i_str;
      //   
      //   // The outlier allocation
      //   out_lab_file_loc = out_lab_file + i_str;
      //   
      //   // The component allocation probabilities
      //   alloc_1_file_loc = alloc_1_file + i_str;
      //   alloc_2_file_loc = alloc_2_file + i_str;
      //   
      //   // Save the results
      //   clust_labels_gauss.save(gauss_lab_file_loc);
      //   clust_labels_cat.save(cat_lab_file_loc);
      //   outlier_vec.save(out_lab_file_loc);
      //   
      //   alloc_prob_gauss_curr.save(alloc_1_file_loc);
      //   alloc_prob_cat_curr.save(alloc_2_file_loc);
      //   
      //   // If not saving to file record results in appropriate matrices
      // } else {
        gaussian_record.col(record_iter) = clust_labels_gauss;
        categorical_record.col(record_iter) = clust_labels_cat;
        outlier_probs_saved.col(record_iter) = outlier_vec;
        
        alloc_prob_gauss += alloc_prob_gauss_curr;
        alloc_prob_cat += alloc_prob_cat_curr;
        
      // }
      
      // Record posteriors of parameters for Gaussian and Categorical
      // distributions
      // if(record_posteriors){
      //   for(arma::uword j = 0; j < n_clust_gauss; j++){
      //     if(save_results){
      //       // Save to file
      //       
      //       // The current component
      //       comp_str = std::to_string(j + 1);
      //       mu_file_loc = mean_file + comp_str + "/iter_" + i_str;
      //       var_file_loc = var_file + comp_str + "/iter_" + i_str;
      //       
      //       loc_mu_variance(1).slice(j).save(mu_file_loc);
      //       loc_mu_variance(0).slice(j).save(var_file_loc);
      //     } else {
      //       // mu(record_iter, j) = loc_mu_variance(1).slice(j);
      //       // variance(record_iter, j) = loc_mu_variance(0).slice(j);
      //       
      //       mu.slice(j).col(record_iter) = loc_mu_variance(1).slice(j); 
      //       variance(j).slice(record_iter) = loc_mu_variance(0).slice(j);
      //     }
      //   }
      //   
      //   for(arma::uword j = 0; j < n_clust_cat; j++){
      //     
      //     // Create the component identifying part of the string
      //     comp_str = std::to_string(j + 1);
      //     
      //     for(arma::uword k = 0; k < n_cols_cat; k++){
      //       
      //       // If saving posteriors to file, save the component relevant parts
      //       if(save_results){
      //         
      //         // For each component save the component specific information in 
      //         // a single field to make for more concise information
      //         class_probs_comp(k) = arma::trans(class_probabilities(k).row(j));
      //         
      //         // If we have recorded all of the probabilities for each variable
      //         // (i.e. we are in the last iteraion over the number of columns)
      //         // save the file
      //         if(k == n_cols_cat - 1){
      //           // Create a file name for the current component with a label of 
      //           // the current iteration
      //           class_probs_file_loc = class_probs_file + comp_str + "/iter_" + i_str;
      //           class_probs_comp.save(class_probs_file_loc);
      //         }
      //         
      //         // If not saving the results record the posterior in the 
      //         // predesignated field
      //       } else{
      //         class_probabilities_saved(j)(k).row(record_iter) = class_probabilities(k).row(j);
      //       }
      //     }
      //   }
      // }
      record_iter++;
    }
  }
  
  // Loading saved posterior objects
  // if(save_results){
  //   for(arma::uword i = 0; i < eff_count; i++){
  //     
  //     
  //     i_str = std::to_string(i);
  //     
  //     gauss_lab_file_loc = gauss_lab_file + i_str;
  //     cat_lab_file_loc = cat_lab_file + i_str;
  //     out_lab_file_loc = out_lab_file + i_str;
  //     
  //     clust_labels_gauss.load(gauss_lab_file_loc);
  //     clust_labels_cat.load(cat_lab_file_loc);
  //     outlier_vec.load(out_lab_file_loc);
  //     
  //     gaussian_record.col(i) = clust_labels_gauss;
  //     categorical_record.col(i) = clust_labels_cat;
  //     outlier_probs_saved.col(i) = outlier_vec;
  //     
  //     // The component allocation probabilities
  //     alloc_1_file_loc = alloc_1_file + i_str;
  //     alloc_2_file_loc = alloc_2_file + i_str;
  //     
  //     alloc_prob_gauss_curr.load(alloc_1_file_loc);
  //     alloc_prob_cat_curr.load(alloc_2_file_loc);
  //     
  //     alloc_prob_gauss += alloc_prob_gauss_curr;
  //     alloc_prob_cat += alloc_prob_cat_curr;
  //     
  //     if(record_posteriors){
  //       for(arma::uword j = 0; j < n_clust_gauss; j++){
  //         
  //         // The current component
  //         comp_str = std::to_string(j + 1);
  //         mu_file_loc = mean_file + comp_str + "/iter_" + i_str;
  //         var_file_loc = var_file + comp_str + "/iter_" + i_str;
  //         
  //         loc_mu_variance(1).slice(j).load(mu_file_loc);
  //         loc_mu_variance(0).slice(j).load(var_file_loc);
  //         
  //         // mu(i, j) = loc_mu_variance(1).slice(j);
  //         // variance(i, j) = loc_mu_variance(0).slice(j);
  //         
  //         mu.slice(j).col(i) = loc_mu_variance(1).slice(j); 
  //         variance(j).slice(i) = loc_mu_variance(0).slice(j);
  //       }
  //       
  //       for(arma::uword j = 0; j < n_clust_cat; j++){
  //         // The current component
  //         comp_str = std::to_string(j + 1);
  //         class_probs_file_loc = class_probs_file + comp_str + "/iter_" + i_str;
  //         class_probs_comp.load(class_probs_file_loc);
  //         
  //         for(arma::uword k = 0; k < n_cols_cat; k++){
  //           // For each component save the component specific information in 
  //           // a single field to make for more concise information
  //           // class_probabilities(k).row(j) = arma::trans(class_probs_comp(k));
  //           class_probabilities_saved(j)(k).row(i) = arma::trans(class_probs_comp(k));
  //           
  //         }
  //       }
  //     }
  //   }
  // }
  
  // Normalise the allocation probabilities
  alloc_prob_gauss = alloc_prob_gauss / eff_count;
  alloc_prob_cat = alloc_prob_cat / eff_count;
  
  // construct similarity matrices
  arma::mat sim(n, n); 
  arma::mat cat_sim(n, n);
  sim = CreateSimilarityMat(gaussian_record);
  cat_sim = CreateSimilarityMat(categorical_record);
  
  return List::create(Named("similarity_1") = sim,
                      Named("similarity_2") = cat_sim,
                      Named("allocation_mat_1") = alloc_prob_gauss,
                      Named("allocation_mat_2") = alloc_prob_cat,
                      Named("class_record_1") = gaussian_record,
                      Named("class_record_2") = categorical_record,
                      // Named("mean_posterior") = mu,
                      // Named("variance_posterior") = variance,
                      // Named("class_prob_posterior") = class_probabilities_saved,
                      Named("context_similarity") = context_similarity_record,
                      Named("entropy") = entropy_cw,
                      Named("outlier") = outlier_probs_saved); // ,
  // Named("likelihood") = mdi_recorded_likelihood);
}
