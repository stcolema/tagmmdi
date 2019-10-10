# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// MDI clustering for two gaussian datasets
// [[Rcpp::export]]
Rcpp::List mdi_gauss_gauss(arma::mat data_1,
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
                           double a0,
                           double b0,
                           arma::uword num_iter,
                           arma::uword burn,
                           arma::uword thinning,
                           bool outlier_1 = false,
                           double t_df_1 = 4.0,
                           bool outlier_2 = false,
                           double t_df_2 = 4.0,
                           bool record_posteriors = false,
                           bool normalise_1 = false,
                           bool normalise_2 = false,
                           double u_1 = 2,
                           double v_1 = 10,
                           double u_2 = 2,
                           double v_2 = 10,
                           arma::uword rate_1_0 = 1,
                           arma::uword rate_2_0 = 1,
                           bool save_results = false,
                           bool load_results = false,
                           arma::uword num_load = 0
) {
  
  
  
  // As we tend to iterate to less than num_iter, add 1 to itself
  // Done here as many objects are sized based on this.
  num_iter++;
  
  // Declare the sample size and dimensionality of the continuous and 
  // categorical data
  arma::uword n = data_1.n_rows;
  arma::uword n_cols_1 = data_1.n_cols;
  
  // arma::uword n_cat = categorical_data.n_rows;
  arma::uword n_cols_2 = data_2.n_cols;
  
  // Frequently will compare clusters across contexts so this is a useful bound
  // to iterate to
  arma::uword min_n_clust = std::min(n_clust_1, n_clust_2);
  
  //  Normalise the continuous data
  if(normalise_1){
    data_1 = arma::normalise(data_1);
  }
  
  if(normalise_2){
    data_2 = arma::normalise(data_2);
  }
  
  // Declare global variance and mean - used in outlier t-distribution 
  // (if included)
  arma::mat global_variance_1(n_cols_1, n_cols_1);
  global_variance_1 = 0.5 * arma::cov(data_1); // Olly's rec
  
  arma::vec global_mean_1(n_cols_1);
  global_mean_1 = arma::trans(arma::mean(data_1, 0));
  
  arma::mat global_variance_2(n_cols_2, n_cols_2);
  global_variance_2 = 0.5 * arma::cov(data_2); // Olly's rec
  
  arma::vec global_mean_2(n_cols_2);
  global_mean_2 = arma::trans(arma::mean(data_2, 0));
  
  // strategic latent variable
  double v = 1.0;
  
  // Cluster weights for each dataset
  arma::vec clust_weights_1(n_clust_1);
  clust_weights_1.zeros();
  
  arma::vec clust_weights_2(n_clust_2);
  clust_weights_2.zeros();
  
  // These are the local cubes of posterior mean and variance overwritten each
  // iteration
  arma::field<arma::cube> loc_mu_variance_1(2);
  
  // Initialise the netroies within this field, defining their dimension 
  // Required for loading
  for(arma::uword i = 0; i < n_clust_1; i++){
    loc_mu_variance_1(0) = arma::cube(n_cols_1, n_cols_1, n_clust_1);
    loc_mu_variance_1(0).zeros();
    loc_mu_variance_1(1) = arma::cube(n_cols_1, 1, n_clust_1);
    loc_mu_variance_1(1).zeros();
  }
  
  arma::mat mu_n_1(n_cols_1, n_clust_1);
  mu_n_1.zeros();
  arma::cube variance_n_1(n_cols_1, n_cols_1, n_clust_1);
  variance_n_1.zeros();
  
  // And for the second dataset
  arma::field<arma::cube> loc_mu_variance_2(2);
  
  for(arma::uword i = 0; i < n_clust_2; i++){
    loc_mu_variance_2(0) = arma::cube(n_cols_2, n_cols_2, n_clust_2);
    loc_mu_variance_2(0).zeros();
    loc_mu_variance_2(1) = arma::cube(n_cols_2, 1, n_clust_2);
    loc_mu_variance_2(1).zeros();
  }
  
  arma::mat mu_n_2(n_cols_2, n_clust_2);
  mu_n_2.zeros();
  arma::cube variance_n_2(n_cols_2, n_cols_2, n_clust_2);
  variance_n_2.zeros();
  
  // Context similarity - sample prior
  double phi = arma::randg(arma::distr_param(a0, 1/b0) );
  
  // Declare normalising constant
  double Z = 1.0;
  
  // Used in each iteration
  arma::vec curr_prob_vec_1(n_clust_1);
  arma::vec curr_prob_vec_2(n_clust_2);
  
  // Various objects to record values for posterior distributions and clustering
  // the record for similarity in each clustering
  arma::uword eff_count = 0;
  
  // To handle loading old runs
  if(num_load > burn){
    eff_count = ceil((double)(num_iter) / (double)thinning) + num_load - 1;
  } else {
    eff_count = ceil((double)(num_iter + num_load - burn) / (double)thinning) + num_load;
  }
  
  // Increase these values by the number of pre-existing objects
  num_iter += num_load;
  
  // These are the lists recording the posterior mean and 
  // variance for each class for each recorded iteration
  arma::field<arma::cube> variance_1(n_clust_1);
  arma::cube cluster_var_record_1(n_cols_1, n_cols_1, eff_count);
  cluster_var_record_1.zeros();
  for(arma::uword i = 0; i < n_clust_1; i++){
    variance_1(i) = cluster_var_record_1;
  }
  arma::cube mu_1(n_cols_1, eff_count, n_clust_1);
  mu_1.zeros();
  
  // For the second dataset
  arma::field<arma::cube> variance_2(n_clust_2);
  arma::cube cluster_var_record_2(n_cols_2, n_cols_2, eff_count);
  cluster_var_record_2.zeros();
  for(arma::uword i = 0; i < n_clust_2; i++){
    variance_2(i) = cluster_var_record_2;
  }
  arma::cube mu_2(n_cols_2, eff_count, n_clust_2);
  mu_2.zeros();
  
  // Record the class allocations
  arma::umat record_1(n, eff_count);
  record_1.zeros();
  
  arma::umat record_2(n, eff_count);
  record_2.zeros();
  
  // The vector to record the context similarity parameter
  arma::vec phi_record(eff_count - num_load);
  
  // A positive integer to hold the current count of recorded iterations
  arma::uword record_iter = num_load;
  
  // To hold output of label flipping function
  arma::vec labels_weights_phi(n + n_clust_2 + 1);
  
  // A vector to record the entropy at each iteration
  arma::vec entropy_cw(num_iter - num_load);
  
  // Rate priors for weight smapling from Gamma distn
  arma::vec rate_0_1(n_clust_1);
  arma::vec rate_0_2(n_clust_2);
  
  // Placeholder prior
  rate_0_1.fill(rate_1_0);
  rate_0_2.fill(rate_2_0);
  
  // ## Filenames for saving posteriors ##
  // Create the base file names that we will append with local iteration / 
  // component information
  
  // Discrete allocation values
  std::string dataset_1_lab_file = "output/dataset_1/class/labels_iter_";
  std::string out_1_lab_file = "output/dataset_1/outlier_allocation/outlier_labels_iter_";
  
  std::string dataset_2_lab_file = "output/dataset_2/class/labels_iter_";
  std::string out_2_lab_file = "output/dataset_2/outlier_allocation/outlier_labels_iter_";
  
  // Parameter posteriors
  std::string mean_1_file = "output/dataset_1/mean/mean_";
  std::string var_1_file = "output/dataset_1/variance/var_";
  
  std::string mean_2_file = "output/dataset_2/mean/mean_";
  std::string var_2_file = "output/dataset_2/variance/var_";
  
  // Allocation probability matrices
  std::string alloc_1_file = "output/dataset_1/allocation/alloc_";
  std::string alloc_2_file = "output/dataset_2/allocation/alloc_";
  
  
  // The file name endings that change in each iteration / component
  std::string i_str;
  std::string comp_str;
  
  std::string lab_file_1_loc;
  std::string lab_file_2_loc;
  std::string out_lab_file_1_loc;
  std::string out_lab_file_2_loc;
  
  std::string mu_file_1_loc;
  std::string var_file_1_loc;
  
  std::string mu_file_2_loc;
  std::string var_file_2_loc;
  
  std::string alloc_1_file_loc;
  std::string alloc_2_file_loc;
  
  // OUTLIER COMPONENT
  // Variables to handle outlier component from tagm
  // Vector of 0 and 1 for points assigned to outlier group or not
  arma::uvec outlier_vec_1(n);
  arma::uvec outlier_vec_2(n);
  
  outlier_vec_1.zeros();
  outlier_vec_2.zeros();
  
  // Begin with all non-fixed points (i.e. unlabelled) to outlier component
  if(outlier_1 && arma::any(fix_vec_1)){
    outlier_vec_1 = 1 - fix_vec_1;
  }
  
  if(outlier_2 && arma::any(fix_vec_2)){
    outlier_vec_2 = 1 - fix_vec_2;
  }
  
  // Declare object for counting the number of items in the outlier component
  double b_1 = 0.0; 
  double b_2 = 0.0;
  
  // for recording the probabilities for each class
  arma::cube class_probs_1(eff_count, n_clust_1, n);
  arma::cube class_probs_2(eff_count, n_clust_2, n); 
  
  // the current iterations probabilities (overwritten - saved to the above cube
  // after burn in and as defined by thinning)
  arma::vec curr_class_probs_1(n_clust_1);
  arma::vec norm_likelihoods_1(n_clust_1);
  
  arma::vec curr_class_probs_2(n_clust_2);
  arma::vec norm_likelihoods_2(n_clust_2);
  
  // Class labels of points not currently assigned as outliers
  arma::uvec relevant_labels_1(n);
  arma::uvec relevant_labels_2(n);
  
  // the current iterations probabilities of being an outlier (non-outlier prob
  // is 1 - curr_outlier_prob)
  arma::vec curr_outlier_prob_1(2);
  double outlier_likelihood_1 = 0.0;
  arma::uword predicted_outlier_1 = 0;
  double outlier_weight_1 = 1 - sample_beta(u_1, v_1);
  
  arma::vec curr_outlier_prob_2(2);
  double outlier_likelihood_2 = 0.0;
  arma::uword predicted_outlier_2 = 0;
  double outlier_weight_2 = 1 - sample_beta(u_2, v_2);
  
  // the predicted class assuming the point is not an outlier for the two contexts
  arma::uword predicted_class_1 = 0;
  arma::uword predicted_class_2 = 0;
  
  // This is where we save the outlier labels
  arma::umat outlier_probs_saved_1(n, eff_count);
  arma::umat outlier_probs_saved_2(n, eff_count);
  
  // Declare the variable to hold the likelihood of bein non-outlier
  double predicted_norm_likelihood_1 = 0.0;
  double predicted_norm_likelihood_2 = 0.0;
  
  // The matrices to hold the allocation probabilities for each class
  arma::mat alloc_prob_1(n, n_clust_1);
  arma::mat alloc_prob_2(n, n_clust_2);
  
  // Set all entries to zero as we will be using += rather than = (i.e will 
  // never over-write the initial value held, only add to it)
  alloc_prob_1.zeros();
  alloc_prob_2.zeros();
  
  // Local allocation probability matrix
  arma::mat alloc_prob_curr_1(n, n_clust_1);
  arma::mat alloc_prob_curr_2(n, n_clust_2);
  
  arma::uword start_iter = 0;
  arma::uword phi_count = 0;
  
  // Loading posterior objects from previous runs
  if(load_results && num_load > 0){
    // Only need to load most recent objects
    start_iter = num_load - 1;
    
    std::string start_iter_str = std::to_string(start_iter);
    
    lab_file_1_loc = dataset_1_lab_file + start_iter_str;
    out_lab_file_1_loc = out_1_lab_file + start_iter_str;
    
    lab_file_2_loc = dataset_2_lab_file + start_iter_str;
    out_lab_file_2_loc = out_2_lab_file + start_iter_str;
    
    clust_labels_1.load(lab_file_1_loc);
    outlier_vec_1.load(out_lab_file_1_loc);
    
    clust_labels_2.load(lab_file_2_loc);
    outlier_vec_2.load(out_lab_file_2_loc);
    
    for(arma::uword j = 0; j < n_clust_1; j++){
      
      // The current component
      comp_str = std::to_string(j + 1);
      
      // The relevant files
      mu_file_1_loc = mean_1_file + comp_str + "/iter_" + start_iter_str;
      var_file_1_loc = var_1_file + comp_str + "/iter_" + start_iter_str;
      
      // Update the current value of the means and covariance matrices
      loc_mu_variance_1(1).slice(j).load(mu_file_1_loc);
      loc_mu_variance_1(0).slice(j).load(var_file_1_loc);
      
    }
    
    for(arma::uword j = 0; j < n_clust_2; j++){
      
      // The current component
      comp_str = std::to_string(j + 1);
      
      // The relevant files
      mu_file_2_loc = mean_2_file + comp_str + "/iter_" + start_iter_str;
      var_file_2_loc = var_2_file + comp_str + "/iter_" + start_iter_str;
      
      // Update the current value of the means and covariance matrices
      loc_mu_variance_2(1).slice(j).load(mu_file_2_loc);
      loc_mu_variance_2(0).slice(j).load(var_file_2_loc);
      
    }
    
    
  }
  
  // Initialise v 
  v = 1.0;
  
  for(arma::uword i = num_load; i < num_iter; i++){
    
    // Consider only the labels of points not considered to be outliers
    relevant_labels_1 = clust_labels_1 % (1 - outlier_vec_1);
    relevant_labels_2 = clust_labels_2 % (1 - outlier_vec_2);
    
    // Entropy for graphing convergence
    entropy_cw(i - num_load) = entropy(clust_weights_1);
    
    // ## Sample the parameters for the two datasets ##
    
    // Sample the posterior mean and variance for the first dataset
    loc_mu_variance_1 = mean_variance_sampling(data_1,
                                               relevant_labels_1,
                                               n_clust_1,
                                               df_0_1,
                                               n_cols_1,
                                               scale_0_1,
                                               lambda_0_1,
                                               mu_0_1);
    
    variance_n_1 = variance_sampling(data_1,
                                     relevant_labels_1,
                                     n_clust_1,
                                     df_0_1,
                                     n_cols_1,
                                     scale_0_1,
                                     lambda_0_1,
                                     mu_0_1);
    
    mu_n_1 = mean_sampling(data_1,
                           relevant_labels_1,
                           n_clust_1,
                           n_cols_1,
                           variance_n_1,
                           lambda_0_1,
                           mu_0_1);
    
    // // Sample the posterior mean and variance for the second dataset
    loc_mu_variance_2 = mean_variance_sampling(data_2,
                                               relevant_labels_2,
                                               n_clust_2,
                                               df_0_2,
                                               n_cols_2,
                                               scale_0_2,
                                               lambda_0_2,
                                               mu_0_2);
    
    
    variance_n_2 = variance_sampling(data_2,
                                     clust_labels_2,
                                     n_clust_2,
                                     df_0_2,
                                     n_cols_2,
                                     scale_0_2,
                                     lambda_0_2,
                                     mu_0_2);
    
    mu_n_2 = mean_sampling(data_2,
                           clust_labels_2,
                           n_clust_2,
                           n_cols_2,
                           variance_n_2,
                           lambda_0_2,
                           mu_0_2);
    
    // ## Sample cluster weights for the two datasets ##
    clust_weights_1 = mdi_cluster_weights(clust_weight_priors_1,
                                          rate_0_1,
                                          v,
                                          n_clust_1,
                                          n_clust_2,
                                          clust_weights_2,
                                          clust_labels_1,
                                          clust_labels_2,
                                          phi);
    
    clust_weights_2 = mdi_cluster_weights(clust_weight_priors_2,
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
    Z = calculate_normalising_constant(clust_weights_1,
                                       clust_weights_2,
                                       phi,
                                       n_clust_1,
                                       n_clust_2);
    
    // sample the strategic latent variable, v
    v = arma::randg( arma::distr_param(n, 1.0/Z ) );
    
    // sample the context similarity parameter (as only two contexts this is a
    // single double - number not a drink)
    phi = sample_phi(clust_labels_1,
                     clust_labels_2,
                     clust_weights_1,
                     clust_weights_2,
                     v,
                     n,
                     min_n_clust,
                     a0,
                     b0);
    
    if(outlier_1){
      b_1 = sum(outlier_vec_1);
      outlier_weight_1 = sample_beta(b_1 + u_1, n + v_1 - b_1);
    }
    
    if(outlier_2){
      b_2 = sum(outlier_vec_2);
      outlier_weight_2 = sample_beta(b_2 + u_2, n + v_2 - b_2);
    }
    
    // sample class for each observation
    for(arma::uword j = 0; j < n; j++){
      
      // for each point create the vector of probabilities associated with 
      // assignment to each cluster within the gaussian dataset
      // Calculate the log-likelihoods for each cluster in each context
      norm_likelihoods_1 = mdi_gauss_clust_probs(j,
                                                 data_1,
                                                 n_clust_1,
                                                 mu_n_1,
                                                 variance_n_1,
                                                 // loc_mu_variance_1(1),
                                                 // loc_mu_variance_1(0),
                                                 phi,
                                                 clust_weights_1,
                                                 relevant_labels_1,
                                                 clust_labels_2);
      
      norm_likelihoods_2 = mdi_gauss_clust_probs(j,
                                                 data_2,
                                                 n_clust_2,
                                                 mu_n_2,
                                                 variance_n_2,
                                                 // loc_mu_variance_2(1),
                                                 // loc_mu_variance_2(0),
                                                 phi,
                                                 clust_weights_2,
                                                 relevant_labels_2,
                                                 clust_labels_1);
      
      // Convert to likelihoods, handle overflow and normalise
      curr_prob_vec_1 = over_flow_handling(norm_likelihoods_1);
      curr_prob_vec_2 = over_flow_handling(norm_likelihoods_2);
      
      // Predict the component to which the obersation belongs
      predicted_class_1 = cluster_predictor(curr_prob_vec_1);
      predicted_class_2 = cluster_predictor(curr_prob_vec_2);
      
      // The various probabilities to determine if the observation is considered 
      // an outlier or not (if using TAGM rather than traditional mixtures)
      if(outlier_1){
        
        // Calculate the likelihood associated with the outlier component
        outlier_likelihood_1 =  t_likelihood(arma::trans(data_1.row(j)), 
                                             global_mean_1, 
                                             global_variance_1, 
                                             n_cols_1,
                                             t_df_1);
        
        
        outlier_likelihood_1 += log(outlier_weight_1);
        
        curr_outlier_prob_1(1) = outlier_likelihood_1;
        
        // The normal likelihood for the current class allocation
        predicted_norm_likelihood_1 = normal_likelihood(arma::trans(data_1.row(j)),
                                                        mu_n_1.col(predicted_class_1 - 1),
                                                        variance_n_1.slice(predicted_class_1 - 1),
                                                        // loc_mu_variance_1(1).slice(predicted_class_1 - 1),
                                                        // loc_mu_variance_1(0).slice(predicted_class_1 - 1),
                                                        n_cols_1);
        
        predicted_norm_likelihood_1 += log(1 - outlier_weight_1);
        curr_outlier_prob_1(0) = predicted_norm_likelihood_1;
        
        // Overflow handling
        curr_outlier_prob_1 = exp(curr_outlier_prob_1 - max(curr_outlier_prob_1));
        
        // Normalise the vector
        curr_outlier_prob_1 = curr_outlier_prob_1 / sum(curr_outlier_prob_1);
        
        // Predict if the current observation is an outlier or not  
        predicted_outlier_1 = cluster_predictor(curr_outlier_prob_1) - 1; 
      }
      
      // std::cout << "out of outlier 1.\n";
      
      if(outlier_2){
        
        // Calculate the likelihood associated with the outlier component
        outlier_likelihood_2 =  t_likelihood(arma::trans(data_2.row(j)), 
                                             global_mean_2, 
                                             global_variance_2, 
                                             n_cols_2,
                                             t_df_2);
        
        
        outlier_likelihood_2 += log(outlier_weight_2);
        
        curr_outlier_prob_2(1) = outlier_likelihood_2;
        
        // The normal likelihood for the current class allocation
        predicted_norm_likelihood_2 = normal_likelihood(arma::trans(data_2.row(j)),
                                                        mu_n_2.col(predicted_class_2 - 1),
                                                        variance_n_2.slice(predicted_class_2 - 1),
                                                        // loc_mu_variance_2(1).slice(predicted_class_2 - 1),
                                                        // loc_mu_variance_2(0).slice(predicted_class_2 - 1),
                                                        n_cols_2);
        
        predicted_norm_likelihood_2 += log(1 - outlier_weight_2);
        curr_outlier_prob_2(0) = predicted_norm_likelihood_2;
        
        // Overflow handling
        curr_outlier_prob_2 = exp(curr_outlier_prob_2 - max(curr_outlier_prob_2));
        
        // Normalise the vector
        curr_outlier_prob_2 = curr_outlier_prob_2 / sum(curr_outlier_prob_2);
        
        // Predict if the current observation is an outlier or not  
        predicted_outlier_2 = cluster_predictor(curr_outlier_prob_2) - 1; 
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
    labels_weights_phi = cluster_label_update(clust_labels_1,
                                              clust_labels_2,
                                              clust_weights_1,
                                              clust_weights_2,
                                              n_clust_1,
                                              n_clust_2,
                                              phi,
                                              min_n_clust,
                                              v,
                                              n,
                                              a0,
                                              b0,
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
      
      if(save_results){
        // Save to file
        
        // Make sure to avoid overwriting the loaded entries
        // record_iter += num_load - 1;
        
        // Create appending label of the current iteration
        i_str = std::to_string(record_iter);
        
        // Create the file names for saving the current iteration's information
        // The discrete componenet allocations
        lab_file_1_loc = dataset_1_lab_file + i_str;
        lab_file_2_loc = dataset_2_lab_file + i_str;
        
        // The outlier allocation
        out_lab_file_1_loc = out_1_lab_file + i_str;
        out_lab_file_2_loc = out_2_lab_file + i_str;
        
        // The component allocation probabilities
        alloc_1_file_loc = alloc_1_file + i_str;
        alloc_2_file_loc = alloc_2_file + i_str;
        
        // Save the results
        clust_labels_1.save(lab_file_1_loc);
        outlier_vec_1.save(out_lab_file_1_loc);
        
        clust_labels_2.save(lab_file_2_loc);
        outlier_vec_2.save(out_lab_file_2_loc);
        
        alloc_prob_curr_1.save(alloc_1_file_loc);
        alloc_prob_curr_2.save(alloc_2_file_loc);
        
        // If not saving to file record results in appropriate matrices
      } else {
        // std::cout << "Here.\n";
        record_1.col(record_iter) = clust_labels_1;
        record_2.col(record_iter) = clust_labels_2;
        outlier_probs_saved_1.col(record_iter) = outlier_vec_1;
        outlier_probs_saved_2.col(record_iter) = outlier_vec_2;
        
        alloc_prob_1 += alloc_prob_curr_1;
        alloc_prob_2 += alloc_prob_curr_2;
      }
      
      // Record posteriors of parameters for Gaussian and Categorical
      // distributions
      if(record_posteriors){
        std::cout << "In here.\n";
        for(arma::uword j = 0; j < n_clust_1; j++){
          if(save_results){
            // Save to file
            
            // The current component
            comp_str = std::to_string(j + 1);
            mu_file_1_loc = mean_1_file + comp_str + "/iter_" + i_str;
            var_file_1_loc = var_1_file + comp_str + "/iter_" + i_str;
            
            loc_mu_variance_1(1).slice(j).save(mu_file_1_loc);
            loc_mu_variance_1(0).slice(j).save(var_file_1_loc);
          } else {
            mu_1.slice(j).col(record_iter) = loc_mu_variance_1(1).slice(j); 
            variance_1(j).slice(record_iter) = loc_mu_variance_1(0).slice(j);
          }
        }
        
        for(arma::uword j = 0; j < n_clust_2; j++){
          if(save_results){
            // Save to file
            
            // The current component
            comp_str = std::to_string(j + 1);
            mu_file_2_loc = mean_2_file + comp_str + "/iter_" + i_str;
            var_file_2_loc = var_2_file + comp_str + "/iter_" + i_str;
            
            loc_mu_variance_2(1).slice(j).save(mu_file_2_loc);
            loc_mu_variance_2(0).slice(j).save(var_file_2_loc);
          } else {
            mu_2.slice(j).col(record_iter) = loc_mu_variance_2(1).slice(j); 
            variance_2(j).slice(record_iter) = loc_mu_variance_2(0).slice(j);
          }
        }
      }
      record_iter++;
    }
  }
  
  // Loading saved posterior objects
  if(save_results){
    for(arma::uword i = 0; i < eff_count; i++){
      
      
      i_str = std::to_string(i);
      
      lab_file_1_loc = dataset_1_lab_file + i_str;
      out_lab_file_1_loc = out_1_lab_file + i_str;
      
      lab_file_2_loc = dataset_2_lab_file + i_str;
      out_lab_file_2_loc = out_2_lab_file + i_str;
      
      clust_labels_1.load(lab_file_1_loc);
      outlier_vec_1.load(out_lab_file_1_loc);
      
      clust_labels_2.load(lab_file_2_loc);
      outlier_vec_2.load(out_lab_file_2_loc);
      
      record_1.col(i) = clust_labels_1;
      record_2.col(i) = clust_labels_2;
      
      outlier_probs_saved_1.col(i) = outlier_vec_1;
      outlier_probs_saved_2.col(i) = outlier_vec_2;
      
      // The component allocation probabilities
      alloc_1_file_loc = alloc_1_file + i_str;
      alloc_2_file_loc = alloc_2_file + i_str;
      
      alloc_prob_curr_1.load(alloc_1_file_loc);
      alloc_prob_curr_2.load(alloc_2_file_loc);
      
      alloc_prob_1 += alloc_prob_curr_1;
      alloc_prob_2 += alloc_prob_curr_2;
      
      if(record_posteriors){
        for(arma::uword j = 0; j < n_clust_1; j++){
          
          // The current component
          comp_str = std::to_string(j + 1);
          
          // The relevant files
          mu_file_1_loc = mean_1_file + comp_str + "/iter_" + i_str;
          var_file_1_loc = var_1_file + comp_str + "/iter_" + i_str;
          
          // Update the current value of the means and covariance matrices
          loc_mu_variance_1(1).slice(j).load(mu_file_1_loc);
          loc_mu_variance_1(0).slice(j).load(var_file_1_loc);
          
          mu_1.slice(j).col(i) = loc_mu_variance_1(1).slice(j); 
          variance_1(j).slice(i) = loc_mu_variance_1(0).slice(j);
        }
        
        for(arma::uword j = 0; j < n_clust_2; j++){
          
          // The current component
          comp_str = std::to_string(j + 1);
          
          // The relevant files
          mu_file_2_loc = mean_2_file + comp_str + "/iter_" + i_str;
          var_file_2_loc = var_2_file + comp_str + "/iter_" + i_str;
          
          // Update the current value of the means and covariance matrices
          loc_mu_variance_2(1).slice(j).load(mu_file_2_loc);
          loc_mu_variance_2(0).slice(j).load(var_file_2_loc);
          
          mu_2.slice(j).col(i) = loc_mu_variance_2(1).slice(j); 
          variance_2(j).slice(i) = loc_mu_variance_2(0).slice(j);
        }
      }
    }
  }
  
  // Normalise the allocation probabilities
  alloc_prob_1 = alloc_prob_1 / eff_count;
  alloc_prob_2 = alloc_prob_2 / eff_count;
  
  // construct similarity matrices
  arma::mat sim_1(n, n); 
  arma::mat sim_2(n, n);
  sim_1 = similarity_mat(record_1);
  sim_2 = similarity_mat(record_2);
  
  return List::create(Named("similarity_1") = sim_1,
                      Named("similarity_2") = sim_2,
                      Named("allocation_mat_1") = alloc_prob_1,
                      Named("allocation_mat_2") = alloc_prob_2,
                      Named("class_record_1") = record_1,
                      Named("class_record_2") = record_2,
                      Named("mean_posterior_1") = mu_1,
                      Named("variance_posterior_1") = variance_1,
                      Named("mean_posterior_2") = mu_2,
                      Named("variance_posterior_2") = variance_2,
                      Named("context_similarity") = phi_record,
                      Named("entropy") = entropy_cw,
                      Named("outlier_1") = outlier_probs_saved_1,
                      Named("outlier_2") = outlier_probs_saved_2); // ,
  // Named("likelihood") = mdi_recorded_likelihood);
}