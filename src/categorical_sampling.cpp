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
  double out = 0.0;
  
  for (arma::uword i = 0; i < num_iter; i++){
    if(cluster_record(point, i) == cluster_record(comparison_point, i)){
      out++;
    }
    
  }
  out = out / num_iter;
  return out;
}

// Constructs a similarity matrix comparing all points clustering across the 
// iterations
// [[Rcpp::export]]
arma::mat similarity_mat(arma::umat cluster_record){
  arma::uword sample_size = cluster_record.n_rows;
  arma::uword num_iter = cluster_record.n_cols;
  arma::mat out(sample_size, sample_size);
  
  for (arma::uword point = 0; point < sample_size; point++){ // if not doing diagonal, restrict to sample size - 1
    for (arma::uword comparison_point = point; // + 1; 
         comparison_point < sample_size;
         comparison_point++){
      out(point, comparison_point) = point_similarity(point, 
          comparison_point,
          cluster_record,
          num_iter);
      out(comparison_point, point) = out(point, comparison_point);
    }
  }
  return out;
}

// Calculate the entropy for the current cluster weights
// [[Rcpp::export]]
double entropy(arma::vec class_weights){
  arma::uword n = class_weights.n_elem;
  arma::vec entropy_components(n);
  // std::cout << "\nDeclared\n";
  
  
  for(arma::uword i = 0; i < n; i++){
    entropy_components(i) = - class_weights(i) * log(class_weights(i));
    if (entropy_components.has_nan()){
      entropy_components(i) = 0;
    }
  }
  // std::cout << "Inter";
  double entropy_out = sum(entropy_components);
  return entropy_out;
}

// === Gaussian ================================================================

// Returns a variable involved in updating the scale, it is similar to sample 
// covariance
arma::mat S_n(arma::mat data,
              arma::vec sample_mean,
              arma::uword sample_size,
              arma::uword num_cols){
  arma::mat sample_covariance(num_cols, num_cols);
  sample_covariance.zeros();
  if(sample_size > 0){
    for(arma::uword i = 0; i < sample_size; i++){
      
      arma::vec row_i = trans(data.row(i));
      
      sample_covariance = (sample_covariance 
                             + ((row_i - sample_mean) 
                                  * arma::trans(row_i - sample_mean)
                             )
      );
      
    }
  }
  return sample_covariance;
}

// Returns a vector of the mean of a cluster of size n
arma::vec mean_n(double lambda_0,
                 arma::vec mu_0,
                 arma::uword sample_size,
                 arma::uword num_cols,
                 arma::vec sample_mean){
  arma::vec mu_n(num_cols);
  mu_n = ((lambda_0 * mu_0 + sample_size * sample_mean)
            / (lambda_0 + sample_size));
  return mu_n;
}

// Returns the matrix for the scale hyperparameter for n observations
arma::mat scale_n(arma::mat scale_0,
                  arma::vec mu_0,
                  double lambda_0,
                  arma::mat sample_covariance,
                  arma::uword sample_size,
                  arma::vec sample_mean){
  arma::mat scale_out;
  
  
  if(sample_size > 0){
    scale_out = (scale_0
                   + sample_covariance
                   + ((lambda_0 * sample_size) / (lambda_0 + sample_size))
                   * (sample_mean - mu_0) * arma::trans(sample_mean - mu_0)
    );
    return scale_out;
  }
  
  scale_out = scale_0; 
  return scale_out;
}

// sample a multivariate normal
arma::mat mvrnormArma(arma::uword n,
                      arma::vec mu,
                      arma::mat sigma) {
  arma::uword ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() +
    Y * arma::chol(sigma);
}

// sample from the mean posterior
arma::vec mean_posterior(arma::vec mu_0,
                         arma::mat variance,
                         double lambda_0,
                         arma::mat data){
  arma::uword ncols = data.n_cols;
  arma::vec sample_mean(ncols); sample_mean.zeros();
  arma::vec mu_out(ncols);
  arma::vec mu_n(ncols);
  arma::mat variance_n(ncols, ncols);
  
  arma::uword sample_size = data.n_rows;
  if (sample_size > 0){
    sample_mean = trans(arma::mean(data, 0));
  }
  
  // std::cout << "\nPast initial if\n";
  
  double lambda_n = lambda_0 + sample_size;
  
  mu_n = mean_n(lambda_0, mu_0, sample_size, ncols, sample_mean);
  
  // std::cout << "\nPast initial mu_n\n";
  // std::cout << "\n" << mu_n << "\n";
  
  variance_n = variance / lambda_n;
  
  // std::cout << "\n" << variance_n << "\n";
  
  arma::mat x = mvrnormArma(1, mu_n, variance_n);
  
  mu_out = arma::conv_to<arma::vec>::from(x.row(0));
  
  return mu_out;
  
}


// Sample the variance for after n observations
arma::mat variance_posterior(int df_0,
                             arma::mat scale_0,
                             double lambda_0,
                             arma::vec mu_0,
                             arma::mat data){
  
  arma::uword sample_size = data.n_rows, num_cols = data.n_cols;
  
  // std::cout << "\nCluster data:\n" << data;
  
  arma::vec sample_mean(num_cols); sample_mean.zeros();
  
  int df_n = df_0 + sample_size;
  
  // std::cout << "Reached another dec\n";
  
  arma::mat scale_n_value(num_cols, num_cols);
  arma::mat sample_covariance(num_cols, num_cols);
  arma::mat variance_out(num_cols, num_cols);
  
  // std::cout << sample_size << "\n";
  
  if (sample_size > 0){
    
    sample_mean = arma::trans(arma::mean(data, 0));
    
  } else{
    sample_mean.fill(0.0);
  }
  
  // std::cout << "\nSample covariance reached\n";
  
  sample_covariance = S_n(data, sample_mean, sample_size, num_cols);
  
  // std::cout << "Scale_n reached\n";
  
  arma::mat samp_cov(num_cols, num_cols);
  samp_cov = (sample_size - 1) * arma::cov(data);
  
  // std::cout  << samp_cov << "\n\n";
  
  scale_n_value = scale_n(
    scale_0,
    mu_0,
    lambda_0,
    sample_covariance,
    sample_size,
    sample_mean
  );
  
  
  // std::cout << scale_n_value << "\n";
  
  variance_out = arma::iwishrnd(scale_n_value, df_n);
  
  return variance_out;
  
} 

// Returns a 2-D field of cubes (i.e. a 5D object) for the current iterations 
// means and variances across each cluster (hence a field)
arma::field<arma::cube> mean_variance_sampling(arma::mat data,
                                               arma::uvec cluster_labels,
                                               arma::uword k,
                                               int df_0,
                                               arma::uword num_cols,
                                               arma::mat scale_0,
                                               double lambda_0,
                                               arma::vec mu_0){
  arma::mat cluster_data;
  arma::mat variance(num_cols, num_cols);
  arma::vec mu(num_cols);
  
  arma::field<arma::cube> mean_variance_field(2);
  
  arma::cube var_entry = arma::zeros<arma::cube>(num_cols, num_cols, k);
  arma::cube mu_entry = arma::zeros<arma::cube>(num_cols, 1, k);
  
  mean_variance_field(0) = var_entry;
  mean_variance_field(1) = mu_entry;
  
  for (arma::uword j = 1; j < k + 1; j++) {
    // std::cout << "\nj for loop";
    cluster_data = data.rows(find(cluster_labels == j ));
    
    // std::cout << mean_variance_field(1).slice(j - 1) << "\n";
    
    // std::cout << "\n" << df_0 << "\n";
    
    // if(cluster_data.n_rows < 10){
      // std::cout << cluster_data << "\n\n";
    // }
    
    mean_variance_field(0).slice(j - 1) = variance_posterior(
      df_0,
      scale_0,
      lambda_0,
      mu_0,
      cluster_data
    );
    
    // std::cout << "\nVariance sampled\n";
    
    // std::cout << mean_variance_field(1).slice(j - 1) << "\n";
    
    mean_variance_field(1).slice(j - 1) = mean_posterior(mu_0, 
                        mean_variance_field(0).slice(j - 1), 
                        lambda_0,
                        cluster_data);
    
    // std::cout << "\nAccessed cubes";
  }
  return mean_variance_field;
}

// === Dirichlet ===============================================================

// update the concentration parameter in the Dirichlet distribution
arma::vec concentration_n(arma::vec concentration_0,
                          arma::uvec cluster_labels,
                          arma::uword num_cat){
  
  // arma::uword n = cluster_labels.n_elem;
  arma::uvec class_members;
  arma::uword class_count;
  arma::vec concentration(num_cat);
  
  for (arma::uword i = 1; i < num_cat + 1; i++) {
    class_count = 0;
    class_members = find(cluster_labels == i);
    class_count = class_members.n_elem;

    
    concentration(i - 1) = arma::as_scalar(concentration_0(i - 1)) + class_count;
  }
  return concentration;
}

arma::vec concentration_n_class(arma::vec concentration_0,
                                arma::uvec cluster_labels,
                                arma::uword num_cat){
  
  // arma::uword n = cluster_labels.n_elem;
  
  arma::uword class_count = 0;
  arma::vec concentration(num_cat);
  
  arma::uvec class_members;
  
  for (arma::uword i = 0; i < num_cat; i++) {
    // class_count = 0;
    class_members = find(cluster_labels == i);
    class_count = class_members.n_elem;
    // for (arma::uword j = 0; j < n; j++ ) {
    //   if (cluster_labels(j) == i) {
    //     class_count++;
    //   }
    // }
    
    concentration(i) = arma::as_scalar(concentration_0(i)) + class_count;
  }
  return concentration;
}

// sample parameters for a dirichlet distribution (normally for the clusters)
// [[Rcpp::export]]
arma::vec dirichlet_posterior(arma::vec concentration_0,
                              arma::uvec cluster_labels,
                              arma::uword num_clusters){
  arma::vec cluster_weight = arma::zeros<arma::vec>(num_clusters);
  
  arma::vec concentration(num_clusters);
  concentration = concentration_n(concentration_0,
                                  cluster_labels,
                                  num_clusters);
  
  
  for (arma::uword i = 0; i < num_clusters; i++) {
    
    // cluster_weight(i) = Rf_rgamma(arma::as_scalar(concentration(i)), 1);
    cluster_weight(i) = arma::randg( arma::distr_param(arma::as_scalar(concentration(i)), 1.0) );
  }
  
  double total_cluster_weight = sum(cluster_weight);
  cluster_weight = cluster_weight / total_cluster_weight;
  return cluster_weight;
}

// Dirichlet posterior for class weights (difference of base 0 compared to vanilla
// dirichlet_posterior function
arma::vec dirichlet_posterior_class(arma::vec concentration_0,
                                    arma::uvec cluster_labels,
                                    arma::uword num_clusters){
  arma::vec cluster_weight = arma::zeros<arma::vec>(num_clusters);
  
  arma::vec concentration(num_clusters);
  concentration = concentration_n_class(concentration_0,
                                        cluster_labels,
                                        num_clusters);
  
  
  for (arma::uword i = 0; i < num_clusters; i++) {
    
    // cluster_weight(i) = Rf_rgamma(arma::as_scalar(concentration(i)), 1);
    // double a = arma::as_scalar(concentration(i));
    cluster_weight(i) = arma::randg(arma::distr_param(arma::as_scalar(concentration(i)), 1.0) );
  }
  
  double total_cluster_weight = sum(cluster_weight);
  cluster_weight = cluster_weight / total_cluster_weight;
  return cluster_weight;
}


//  TEST FUNCTIONS FOR CLASSES WITH NULL OF 0 //////////////////////////////////

// MDI class weights (Gamma rather than dirichlet)
arma::vec gamma_posterior(arma::vec concentration_0,
                          arma::uvec cluster_labels,
                          arma::uword num_clusters,
                          double rate
                          ){
  arma::vec cluster_weight = arma::zeros<arma::vec>(num_clusters);
  
  arma::vec concentration(num_clusters);
  concentration = concentration_n(concentration_0,
                                  cluster_labels,
                                  num_clusters);
  
  for (arma::uword i = 0; i < num_clusters; i++) {
    cluster_weight(i) = arma::randg( arma::distr_param(arma::as_scalar(concentration(i)), 1.0) );
  }
  return cluster_weight;
}

double mdi_cluster_rate(double v,
                        arma::uword n_clust_comp,
                        arma::uword cluster_index,
                        arma::vec cluster_weights_comp,
                        double phi){
  double b = 0.0;
  
  // b = arma::sum(cluster_weights_comp) 
  //     + cluster_weights_comp(cluster_index) * phi;
  
  for(arma::uword i = 0; i < n_clust_comp; i++){
    b += cluster_weights_comp(i) * (1 + phi * (cluster_index == i));
  }

  b = b * v;
  
  return b;
}

arma::vec mdi_cluster_weights(arma::vec shape_0,
                              arma::vec rate_0,
                              double v,
                              arma::uword n_clust,
                              arma::uword n_clust_comp,
                              arma::vec cluster_weights_comp,
                              arma::uvec cluster_labels,
                              arma::uvec cluster_labels_comp,
                              double phi){
  
  arma::vec cluster_weight = arma::zeros<arma::vec>(n_clust);
  
  double b = 0.0;
  // double b1 = 0.0; 

  arma::vec shape_n(n_clust);
  
  shape_n = concentration_n(shape_0,
                            cluster_labels,
                            n_clust) + 1;
  
  b = v * arma::sum(cluster_weights_comp);
  
  // std::cout << "\nShape: " << shape_n << "\nRate: " << b << "\n";
  
  for (arma::uword i = 0; i < n_clust; i++) {
    
    // if(i < n_clust_comp){
    //   b += v * cluster_weights_comp(i) * phi;
    // }
    
    b = mdi_cluster_rate(v,
                         n_clust_comp,
                         i,
                         cluster_weights_comp,
                         phi);
    // 
    // if (abs(b - b1) > 0.001){
    //   std::cout << "\nRates differ:\n" << b << "\n" << b1 << "\n\n";
    // }
    // 
    // std::cout << "\nShape: " << shape_n(i) << "\nRate: " << b + rate_0(i) 
    //           << "\nStrategic latent variable: " << v << "\n\n";
    
    cluster_weight(i) = arma::randg(arma::distr_param(shape_n(i), 1 / (b + rate_0(i))));
    
    // if(i < n_clust_comp){
    //   b -= v * cluster_weights_comp(i) * phi;
    // }
    
  }
  return cluster_weight;
}

// Several functions to initialise the 3D array required for the different
// classes for each variable within each cluster

// count unique entries in a vector
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
                                                  arma::field<arma::mat> class_probabilities,
                                                  arma::field<arma::vec> phi_prior,
                                                  arma::uvec cluster_labels,
                                                  arma::uvec cat_count,
                                                  arma::uword num_clusters,
                                                  arma::uword num_cols
){

  arma::umat cluster_data;
  // arma::umat indices;
  for(arma::uword k = 1; k < num_clusters + 1; k++){
    
    // std::cout << "In for loop, k = " << k << "\n";
    
    // indices = find(cluster_labels == k);
    
    // std::cout << "\n" << indices << "\n";
    
    cluster_data = data.rows(find(cluster_labels == k));
    
    // std::cout << "Generic message\n";
    
    for(arma::uword j = 0; j < num_cols; j++){
      
      
      class_probabilities(j).row(k - 1) = arma::trans(dirichlet_posterior_class(phi_prior(j),
                                                                            cluster_data.col(j),
                                                                            cat_count(j)
                                                                            )
      );
      
      // std::cout << "Another message\n";
    }
  }
  return class_probabilities;
}

// Sample the cluster membership of point
// [[Rcpp::export]]
arma::vec categorical_cluster_probabilities(arma::urowvec point,
                                            arma::umat data,
                                            arma::field<arma::mat> class_probabilities,
                                            arma::vec cluster_weights,
                                            arma::uword num_clusters,
                                            arma::uword num_cols){
  
  // std::cout << "In function cluster_probabilities\n";
  arma::vec probabilities = arma::zeros<arma::vec>(num_clusters);

  double curr_weight = 0.0;
  
  // std::cout << "\n\n" << class_probabilities << "\n\n";
  for(arma::uword i = 0; i < num_clusters; i++){
    curr_weight = log(cluster_weights(i));
    for(arma::uword j = 0; j < num_cols; j++){

      probabilities(i) = probabilities(i) + std::log(class_probabilities(j)(i, point(j)));
    }
    probabilities(i) = probabilities(i) + curr_weight;
  }
  probabilities = exp(probabilities - max(probabilities));
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
  pred = 1 + sum(u > cumsum(probabilities)); // 1 + 
  return pred;
}

// The actual categorical clustering all wrapped up in one function
// [[Rcpp::export]]
Rcpp::List categorical_clustering(arma::umat data,
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
  
  arma::uword n = data.n_rows;
  arma::uword num_cols = data.n_cols;
  
  arma::uvec cat_count(num_cols);
  cat_count = cat_counter(data);
  
  arma::field<arma::mat> class_probabilities(num_cols);
  
  // std::cout << "Declaring field of matrices for class probs\n";
  
  class_probabilities = declare_class_probs_field(cat_count,
                            num_cols,
                            num_clusters);
  
  // std::cout << "Class probabilities declared\n";
  
  // std::cout << "Num clusters: " << num_clusters << "\n";
  
  arma::vec cluster_weights(num_clusters);
  
  arma::vec curr_cluster_probs(num_clusters);
  
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  
  // std::cout << "N: " << n << "\nEff count: " << eff_count << "\n";
  arma::umat record(n, eff_count);
  record.zeros();
  
  // arma::uword prediction = 0;
  
  // std::cout << "Reached for loop\n";
  
  for(arma::uword i = 0; i < num_iter; i++){
    
    cluster_weights = dirichlet_posterior(cluster_weight_priors,
                                          cluster_labels,
                                          num_clusters);
    
    // std::cout << "Cluster weights calculated\n";
    
    class_probabilities = sample_class_probabilities(data,
                                                     class_probabilities,
                                                     phi_prior,
                                                     cluster_labels,
                                                     cat_count,
                                                     num_clusters,
                                                     num_cols
                                 
    );
    
    // std::cout << "Class probs calculated\n";
    // std::cout << class_probabilities << "\n\n";
    
    for(arma::uword j = 0; j < n; j++){
      
      // sample cluster for each point here
      curr_cluster_probs = categorical_cluster_probabilities(data.row(j),
                                                             data,
                                                             class_probabilities,
                                                             cluster_weights,
                                                             num_clusters,
                                                             num_cols);
      
      // std::cout << "Cluster sampled\n";
      
      if(fix_vec(j) == 0){
        // std::cout << "Prediction: " << prediction << "\n";
        cluster_labels(j) = cluster_predictor(curr_cluster_probs);
      }
    }
    if (i >= burn && (i - burn) % thinning == 0) {
      // std::cout << "\n" << i << "\n";
      record.col((i - burn) / thinning) = cluster_labels;
    }
  }
  arma::mat sim(n, n);
  sim = similarity_mat(record);
  return List::create(Named("similarity") = sim,
                      Named("class_record") = record);
  // return sim;
}

// === Gaussian clustering =====================================================

// arma::vec sample_gaussian_cluster(arma::vec point,
//                                   arma::mat data,
//                                   arma::uword k,
//                                   arma::vec class_weights,
//                                   arma::uvec class_labels,
//                                   arma::cube mu,
//                                   arma::cube variance,
//                                   bool outlier = false,
//                                   arma::vec global_mean = arma::zeros<arma::vec>(1),
//                                   arma::mat global_variance = arma::zeros<arma::mat>(1, 1),
//                                   double t_df = 4.0){
//   
//   double curr_weight;
//   double exponent;
//   double log_likelihood;
//   
//   
//   double log_det;
//   arma::vec prob_vec(k);
//   
//   arma::uvec count_probs;
//   
//   arma::uword d = data.n_cols;
//   
//   for(arma::uword i = 1; i < k + 1; i++){
//     curr_weight = log(class_weights(i - 1));
//     
//     if(outlier && i == k){
//       
//       exponent = arma::as_scalar(arma::trans(point - global_mean) 
//                                  * arma::inv(global_variance)
//                                  * (point - global_mean));
//                                    
//       log_det = arma::log_det(global_variance).real();
//                                    
//       log_likelihood = lgamma((t_df + d)/2.0)
//                        - lgamma(t_df/2.0)
//                        + d/2.0 * log(t_df * M_PI) 
//                        + log_det 
//                        - ((t_df + d)/2.0) * log(1 + (1/t_df) * exponent);
//                                    
//     }
//     else {
//       exponent = arma::as_scalar(arma::trans(point - mu.slice(i - 1)) 
//                                  * arma::inv(variance.slice(i - 1))
//                                  * (point - mu.slice(i - 1)));
//                                    
//       log_det = arma::log_det(variance.slice(i - 1)).real();
//       log_likelihood = -0.5 *(log_det + exponent + d * log(2 * M_PI));
//     }
//     prob_vec(i - 1) = curr_weight + log_likelihood;
//     // std::cout <<  "\nDIRICHLET:\nCluster " << i << "\nProbability: " 
//     //           << prob_vec(i - 1) << "\nWeight: " << exp(curr_weight) 
//     //           << "\nLog likelihood: " << log_likelihood << "\n\n";
//   } 
//   prob_vec = exp(prob_vec - max(prob_vec));
//   prob_vec = prob_vec / sum(prob_vec);
//   
//   return prob_vec;
// }
// 
// The actual clustering/sampling
// // [[Rcpp::export]]
// Rcpp::List gaussian_clustering(arma::uword num_iter,
//                                arma::vec concentration_0,
//                                arma::mat scale_0,
//                                arma::uvec class_labels,
//                                std::vector<bool> fix_vec,
//                                arma::vec mu_0,
//                                double lambda_0,
//                                arma::mat data,
//                                int df_0,
//                                arma::uword k,
//                                arma::uword burn,
//                                arma::uword thinning,
//                                bool outlier = false,
//                                double t_df = 4.0,
//                                bool record_posteriors = false,
//                                bool normalise = false){
//   
//   // std::cout << "In function";
//   arma::uword N;
//   arma::uword num_cols;
//   N = data.n_rows;
//   num_cols = data.n_cols;
//   
//   //  Normalise the continuous data
//   if(normalise){
//     data = arma::normalise(data);
//   }
//   
//   // for use in the outlier distribution
//   arma::mat global_variance(num_cols, num_cols);
//   
//   // std::cout << arma::cov(data) << "\n";
//   
//   global_variance = 0.5 * arma::cov(data); // Olly's rec
//   
//   arma::vec global_mean(num_cols);
//   
//   global_mean = arma::trans(arma::mean(data, 0));
//   
//   arma::vec entropy_cw(num_iter);
//   
//   arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
//   arma::uword record_ind;
//   
//   arma::umat record(N, eff_count);
//   record.zeros();
//   
//   // std::cout << "Record out \n";
//   
//   arma::mat sim(N, N);
//   arma::mat cluster_data;
//   
//   // Add +1 to k to allow outlier class
//   if(outlier){
//     k++;
//   }
//   
//   arma::vec class_weights(k);
//   
//   // These are the lists recording the posterior mean and 
//   // variance for each class for each recorded iteration
//   ListMatrix variance(eff_count, k);
//   ListMatrix mu(eff_count, k);
//   
//   arma::vec point;
//   
//   // std::cout << "Faux output sentence\n";
//   
//   arma::cube class_probs(eff_count, k, N);
//   arma::vec curr_class_probs(k);
//   
//   // std::cout << "Output sentence\n";
//   
//   // These are the local cubes of posterior mean and variance overwritten each
//   // iteration
//   arma::field<arma::cube> loc_mu_variance(2);
//   
//   for(arma::uword i = 0; i < num_iter; i++){
//     
//     class_weights = dirichlet_posterior(concentration_0, class_labels, k);
//     
//     // std::cout << class_weights << "\n\n";
//     // std::cout << "\nENTROPY";
//     
//     entropy_cw(i) = entropy(class_weights);
//     
//     // std::cout << "\nBegin sampling parameters\n";
//     
//     loc_mu_variance = mean_variance_sampling(data,
//                                              class_labels,
//                                              k,
//                                              df_0,
//                                              num_cols,
//                                              scale_0,
//                                              lambda_0,
//                                              mu_0);
//     
//     // std::cout << "\nAccessed cubes\n";
//     
//     for (arma::uword jj = 0; jj < N; jj++){
//       // if the current point is not fixed, sample its class
//       point = arma::trans(data.row(jj));
//       
//       curr_class_probs = sample_gaussian_cluster(point, 
//                                                  data,
//                                                  k, 
//                                                  class_weights, 
//                                                  class_labels,
//                                                  loc_mu_variance(1),
//                                                  loc_mu_variance(0),
//                                                  outlier,
//                                                  global_mean,
//                                                  global_variance,
//                                                  t_df
//       );
//       
//       // std::cout << curr_class_probs << "\n\n";
//       
//       if (i >= burn && (i - burn) % thinning == 0) {
//         // std::cout << "record accessed" << "\n";
//         record_ind = (i - burn) / thinning;
//         class_probs.slice(jj).row(record_ind) = arma::trans(curr_class_probs);
//         
//       }
//       if(! fix_vec[jj]){
//         class_labels(jj) = cluster_predictor(curr_class_probs);
//         // std::cout << "New label\n" << class_labels(jj) << "\n";
//       }
//     }
//     // std::cout << "Labels\n" << class_labels << "\n";
//     // std::cout << "Generic message\n" << "\n";
//     
//     if (i >= burn && (i - burn) % thinning == 0) {
//       // std::cout << i << "\n";
//       
//       record_ind = (i - burn) / thinning;
//       record.col(record_ind) = class_labels;
//       // std::cout << "record accessed" << "\n";
//       
//       if(record_posteriors){
//         for(arma::uword j = 0; j < k; j ++){
//           
//           // std::cout << "Recording params" << j << "\n";
//           
//           mu(record_ind, j) = loc_mu_variance(1).slice(j);
//           variance(record_ind, j) = loc_mu_variance(0).slice(j);
//         }
//       }
//       
//     }
//   }
//   
//   // std::cout << "Issue is here";
//   sim = similarity_mat(record);
//   // std::cout << "Y is here";
//   // return sim;
//   
//   if(record_posteriors){
//     
//     return List::create(Named("similarity") = sim,
//                         Named("class_record") = record,
//                         Named("mean_posterior") = mu,
//                         Named("variance_posterior") = variance,
//                         Named("entropy") = entropy_cw,
//                         Named("class_prob") = class_probs);
//   }
//   
//   return List::create(Named("similarity") = sim,
//                       Named("class_record") = record,
//                       Named("entropy") = entropy_cw,
//                       Named("class_prob") = class_probs);
// }

// Returns the normal distribution log-likelihood
double normal_likelihood(arma::vec point,
                         arma::vec mu,
                         arma::mat variance,
                         arma::uword d){
  double log_det = 0.0;
  double exponent = 0.0;
  double log_likelihood = 0.0;
  
  
  exponent = arma::as_scalar(arma::trans(point - mu) 
                               * arma::inv(variance)
                               * (point - mu));
                               
  log_det = arma::log_det(variance).real();
  log_likelihood = -0.5 *(log_det + exponent + d * log(2 * M_PI));
                               
  return log_likelihood;
}

// For k classes returns a k-vector of probabilities for point to belong to said
// classes
arma::vec sample_gaussian_cluster(arma::vec point,
                                  arma::mat data,
                                  arma::uword k,
                                  arma::vec class_weights,
                                  arma::cube mu,
                                  arma::cube variance){
  
  double curr_weight;
  double log_likelihood;
  
  arma::vec prob_vec(k);
  
  arma::uword d = data.n_cols;
  
  // Calculate the likelihood for each class
  for(arma::uword i = 0; i < k; i++){
    curr_weight = log(class_weights(i));
    
    log_likelihood = normal_likelihood(point, 
                                       mu.slice(i),
                                       variance.slice(i),
                                       d);
    
    prob_vec(i) = curr_weight + log_likelihood;
  } 
  return prob_vec;
}

arma::vec over_flow_handling(arma::vec prob_vec_in){
  
  arma::vec prob_vec;
  
  // Overflow handling
  prob_vec = exp(prob_vec_in - max(prob_vec_in));
  
  // Normalise the vector
  prob_vec = prob_vec / sum(prob_vec);
  
  return prob_vec;
}

// Returns the t-distribution log likelihood for a given point
double t_likelihood(arma::vec point,
                    arma::vec mu,
                    arma::mat variance,
                    arma::uword d,
                    double df
){
  
  double log_det = 0.0;
  double exponent = 0.0;
  double log_likelihood = 0.0;
  
  exponent = arma::as_scalar(arma::trans(point - mu) 
                             * arma::inv(variance)
                             * (point - mu));
                               
  log_det = arma::log_det(variance).real();
  
  log_likelihood = lgamma((df + d)/2.0) 
    - lgamma(df/2.0) 
    - d/2.0 * log(df * M_PI) 
    - 0.5 * log_det 
    - ((df + d)/2.0) * log(1 + (1/df) * exponent);
    
  return log_likelihood;
}

// Samples if the point is an outlier or not (comparing to assigned class)
double sample_outlier(arma::vec point,
                      arma::mat data,
                      double outlier_weight,
                      arma::vec global_mean,
                      arma::mat global_variance,
                      double t_df = 4.0,
                      arma::uword u = 2,
                      arma::uword v = 10){
  
  double log_likelihood = 0.0;
  double class_weight = 0.0;
  double prob = 0.0;
  
  arma::uword n = data.n_rows;
  arma::uword d = data.n_cols;
  
  // The probability of belonging to the outlier class
  class_weight = (outlier_weight + u)/(n + u + v - 1);
  log_likelihood = t_likelihood(point, global_mean, global_variance, d, t_df);
  
  // std::cout << "\nLog outlier weight: " << class_weight << "\nLog-likelihood: "
  //           << log_likelihood << "\n";
  prob = log(class_weight) + log_likelihood;
  
  return prob;
}
//   // Overflow handling
//   prob_vec = exp(prob_vec - max(prob_vec));
//   
//   // Normalise
//   prob_vec = prob_vec / sum(prob_vec);
//   
//   // Handle overflow
//   prob_vec = over_flow_handling(prob_vec);
//   
//   return prob_vec;
// }

// Calculates a, the secon parameter of the Beta posterior for outlier_weight
double calculate_a(arma::vec norm_likelihoods,
                   double outlier_weight,
                   double outlier_likelihood,
                   arma::uword k){
  arma::vec weighted_likelihood(k);
  arma::vec denom(k);
  
  double a = 0.0;
  
  weighted_likelihood = norm_likelihoods * (1 - outlier_weight);
  denom = weighted_likelihood + outlier_likelihood;
  
  a = 1 / k * arma::sum(weighted_likelihood / denom);
  
  return a;
}

// Samples from a Beta distribution based on the idea that two independent gamma
// functions can be used. See wikipedia.
double sample_beta(double a, double b, double theta = 1.0){
  double X = arma::randg( arma::distr_param(a, 1/theta) );
  double Y = arma::randg( arma::distr_param(b, 1/theta) );
  
  double beta = X / (X + Y);
  
  return beta;
}

// [[Rcpp::export]]
Rcpp::List gaussian_clustering(arma::uword num_iter,
                               arma::vec concentration_0,
                               arma::mat scale_0,
                               arma::uvec class_labels,
                               arma::uvec fix_vec,
                               arma::vec mu_0,
                               double lambda_0,
                               arma::mat data,
                               int df_0,
                               arma::uword k,
                               arma::uword burn,
                               arma::uword thinning,
                               bool outlier = false,
                               double t_df = 4.0,
                               bool record_posteriors = false,
                               bool normalise = false,
                               double u = 2,
                               double v = 10){
  
  // std::cout << "In function";
  arma::uword N;
  arma::uword num_cols;
  N = data.n_rows;
  num_cols = data.n_cols;
  
  // To allow using < and keeping in line with object sizes
  num_iter++;
  
  //  Normalise the continuous data
  if(normalise){
    data = arma::normalise(data);
  }
  
  // for use in the outlier distribution
  arma::mat global_variance(num_cols, num_cols);
  
  // std::cout << arma::cov(data) << "\n";
  
  global_variance = 0.5 * arma::cov(data); // Olly's rec
  
  arma::vec global_mean(num_cols);
  
  global_mean = arma::trans(arma::mean(data, 0));
  
  arma::vec entropy_cw(num_iter);
  
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  arma::uword record_ind;
  
  arma::umat record(N, eff_count);
  
  // for recording the probabilities for being declared an outlier
  arma::umat outlier_probs_saved(N, eff_count);
  
  record.zeros();
  
  // std::cout << "Record out \n";
  
  arma::mat sim(N, N);
  arma::mat cluster_data;
  
  
  // These are the lists recording the posterior mean and 
  // variance for each class for each recorded iteration
  ListMatrix variance(eff_count, k);
  ListMatrix mu(eff_count, k);
  
  arma::vec point;
  
  // std::cout << "Faux output sentence\n";
  
  // arma::uvec out_class = 0;
  // 
  // // Add +1 to k to allow outlier class
  // if(outlier){
  //   out_class++;
  //   // k++;
  // }
  // 
  
  arma::vec class_weights(k);
  arma::vec outlier_weights(2);
  
  // for the calculation of the a, b parameters in the posterior of the outlier 
  // class weight (epsilon in Crook et al 2018)
  arma::vec outlier_probs(N);
  
  arma::uvec b_k(N); // the vector of b_k's the sum of which gives b
  double b = 0.0; // a will be N - b, no need to declare
  
  // for recording the probabilities for each class
  arma::cube class_probs(eff_count, k, N); 
  
  // Vector of 0 and 1 for points assigned to outlier group or not
  arma::uvec outlier_vec(N);
  
  // Begin with all non-fixed points (i.e. unlabelled) to outlier component
  outlier_vec = 1 - fix_vec;
  
  // the current iterations probabilities (overwritten - saved to the above cube
  // after burn in and as defined by thinning)
  arma::vec curr_class_probs(k);
  arma::vec curr_norm_likelihoods(k);
  
  // Class labels of points not currently assigned as outliers
  arma::uvec relevant_labels(N);
  
  // the current iterations probabilities of being an outlier (non-outlier prob
  // is 1 - curr_outlier_prob)
  arma::vec curr_outlier_prob(2);
  double curr_outlier_likelihood = 0.0;
  arma::uword predicted_outlier = 0;
  double outlier_weight = 1 - sample_beta(u, v);

  b_k = arma::find(outlier_vec == 0);
  b = b_k.n_elem;
  
  outlier_weight = 1 - sample_beta(u + b, v + N - b);
  
  double predicted_norm_likelihood = 0.0;
  
  // the predicted class assuming the point is not an outlier
  arma::uword predicted_class = 0;
  
  // std::cout << "Output sentence\n";
  
  // These are the local cubes of posterior mean and variance overwritten each
  // iteration
  arma::field<arma::cube> loc_mu_variance(2);
  
  // The matrices to hold the allocation probabilities for each class
  arma::mat alloc_prob_gauss(N, k);

  // Set all entries to zero as we will be using += rather than = (i.e will 
  // never over-write the initial value held, only add to it)
  alloc_prob_gauss.zeros();

  for(arma::uword i = 0; i < num_iter; i++){
    
    // To see which points are relevant for defining component parameters
    // use pairwise multiplication between the current label and the outlier
    relevant_labels = class_labels % (1 - outlier_vec);
    
    // std::cout << "\nRelevant labels:\n" << relevant_labels << "\n";
    
    class_weights = dirichlet_posterior(concentration_0, relevant_labels, k);
    
    // std::cout << class_weights << "\n\n";
    // std::cout << "\nENTROPY";
    
    entropy_cw(i) = entropy(class_weights);

    // std::cout << "\nBegin sampling parameters\n";
    
    loc_mu_variance = mean_variance_sampling(data,
                                             relevant_labels,
                                             k,
                                             df_0,
                                             num_cols,
                                             scale_0,
                                             lambda_0,
                                             mu_0);
    
    // std::cout << "\nAccessed cubes\n";
    
    for (arma::uword jj = 0; jj < N; jj++){
      // if the current point is not fixed, sample its class
      point = arma::trans(data.row(jj));
      
      curr_norm_likelihoods = sample_gaussian_cluster(point, 
                                                      data,
                                                      k, 
                                                      class_weights, 
                                                      loc_mu_variance(1),
                                                      loc_mu_variance(0)
      );
      
      curr_class_probs = over_flow_handling(curr_norm_likelihoods);


      
      if (i >= burn && (i - burn) % thinning == 0) {
        
        // std::cout << "\nNormal probabilities:\n" << curr_class_probs << "\n";
        // std::cout << "record accessed" << "\n";
        record_ind = (i - burn) / thinning;
        class_probs.slice(jj).row(record_ind) = arma::trans(curr_class_probs);
        alloc_prob_gauss.row(jj) += arma::trans(curr_class_probs);
        
      }
      
      predicted_class = cluster_predictor(curr_class_probs);
      
      if(outlier){
        
        // std::cout << "Outlier stuff\n";
        
        curr_outlier_likelihood = sample_outlier(point,
                                                 data,
                                                 outlier_weight,
                                                 global_mean,
                                                 global_variance,
                                                 t_df,
                                                 u = 2 + b,
                                                 v = 10 + N - b);
        
        // curr_outlier_prob(1) = curr_outlier_likelihood;
        // 
        // 
        // 
        // // std::cout << "t likelihood OK!\n";
        // 
        // curr_outlier_prob(0)= curr_norm_likelihoods(predicted_class - 1) 
        //   + log(1 - outlier_weight)
        //   - log(class_weights(predicted_class - 1));
        // 
        // // std::cout << "normal likelihood OK!\n";
        // curr_outlier_prob = over_flow_handling(curr_outlier_prob);
        //   
        curr_outlier_prob(1) = log(1 + exp(curr_outlier_likelihood));
        
        // std::cout << "t likelihood OK!\n";
        
        predicted_norm_likelihood = normal_likelihood(point,
                                                      loc_mu_variance(1).slice(predicted_class - 1),
                                                      loc_mu_variance(0).slice(predicted_class - 1),
                                                      num_cols);
        
        predicted_norm_likelihood += log(1 - outlier_weight);
        curr_outlier_prob(0) = predicted_norm_likelihood;
  

        
        // std::cout << "normal likelihood OK!\n";
        // std::cout << "\nOutlier probs pre-normalising:\n" << curr_outlier_prob 
        // << "\n";
        
        // Overflow handling
        curr_outlier_prob = exp(curr_outlier_prob - max(curr_outlier_prob)) - 1;
        
        // Normalise the vector
        curr_outlier_prob = curr_outlier_prob / sum(curr_outlier_prob);
            
          
          
        // std::cout << "\nCurr outlier prob:\n" << curr_outlier_prob << "\n";
        predicted_outlier = cluster_predictor(curr_outlier_prob) - 1; // as +1 to handle R
        
        // std::cout << "Outlier prediction done.\n";
      }
      
      // std::cout << "b allocation\n";
      // b_k(jj) = 1 - calculate_a(curr_norm_likelihoods,
      //     outlier_weight,
      //     curr_outlier_likelihood,
      //     k);
      
      if(fix_vec[jj] == 0){
        
        class_labels(jj) = predicted_class;
        outlier_vec(jj) = predicted_outlier;
        // std::cout << "New label\n" << class_labels(jj) << "\n";
        
      }
    }
    // std::cout << "Labels\n" << class_labels << "\n";
    // std::cout << "Generic message\n" << "\n";
    
    if (i >= burn && (i - burn) % thinning == 0) {
      // std::cout << i << "\n";
      
      record_ind = (i - burn) / thinning;
      record.col(record_ind) = class_labels;
      outlier_probs_saved.col(record_ind) = outlier_vec;
      
      // std::cout << "record accessed" << "\n";
      
      if(record_posteriors){
        for(arma::uword j = 0; j < k; j ++){
          
          // std::cout << "Recording params" << j << "\n";
          
          mu(record_ind, j) = loc_mu_variance(1).slice(j);
          variance(record_ind, j) = loc_mu_variance(0).slice(j);
        }
      }
      
    }
    // b = sum(b_k);
    b_k = arma::find(outlier_vec == 0);
    b = b_k.n_elem;
    
    outlier_weight = 1 - sample_beta(u + b, v + N - b);

  }
  
  // std::cout << "Issue is here";
  sim = similarity_mat(record);
  alloc_prob_gauss = alloc_prob_gauss / eff_count; 
  // std::cout << "Y is here";
  // return sim;
  
  if(record_posteriors){
    
    return List::create(Named("similarity") = sim,
                        Named("class_record") = record,
                        Named("mean_posterior") = mu,
                        Named("variance_posterior") = variance,
                        Named("entropy") = entropy_cw,
                        Named("allocation_mat_1") = alloc_prob_gauss,
                        Named("outliers") = outlier_probs_saved);
  }
  
  return List::create(Named("similarity") = sim,
                      Named("class_record") = record,
                      Named("entropy") = entropy_cw,
                      Named("allocation_mat_1") = alloc_prob_gauss,
                      Named("outliers") = outlier_probs_saved);
}

// === MDI Code ================================================================

// returns the count of the number of points with the same label in both contexts
arma::uword count_common_cluster(arma::uvec cluster_labels_1,
                                 arma::uvec cluster_labels_2,
                                 arma::uword n){
  arma::uword common_cluster = 0;
  for(arma::uword i = 0; i < n; i++){
    if(cluster_labels_1(i) == cluster_labels_2(i)){
      common_cluster++;
    }
  }
  return common_cluster;
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

// samples a gamma distn for the current iterations context similarity parameter
// (phi in the original 2012 paper).
double sample_phi(arma::uvec cluster_labels_1,
                  arma::uvec cluster_labels_2,
                  arma::vec cluster_weights_1,
                  arma::vec cluster_weights_2,
                  double v,
                  arma::uword n,
                  arma::uword min_num_clusters,
                  double a0,
                  double b0){
  
  // calculate the shape of the relevant gamma function
  arma::uword count_same_cluster = count_common_cluster(cluster_labels_1,
                                                        cluster_labels_2,
                                                        n);
  
  arma::vec prob_vec(count_same_cluster);
  
  // calculate the rate
  double b = observed_rate(v,
                           min_num_clusters,
                           cluster_weights_1,
                           cluster_weights_2)
    + b0;
  
  double phi = 0.0;
  
  // std::cout << "\nSampling phi\nCount same cluster: " << count_same_cluster
  //           << "\nb_n: " << b << "\nv: " << v << "\nSampling done.\n\n";
  
  // context similarity is a weighted sum of gammas
  if(count_same_cluster > 0){
    
    for(arma::uword i = 0; i < count_same_cluster; i++){
      prob_vec(i) = log_factorial(count_same_cluster)
      - log_factorial(i)
      - log_factorial(count_same_cluster - i)
      + log_factorial(i + a0 - 1)
      - (i + a0)*log(b);
    }
    
    prob_vec = exp(prob_vec - max(prob_vec));
    prob_vec = prob_vec / sum(prob_vec);
    
    // If move to general null value of 0 can use the prediction function
    double u;
    u = arma::randu<double>( );
    
    // include + 1 if labels centred on 1
    arma::uword pred_ind = 0;
    
    pred_ind = sum(u > cumsum(prob_vec)); // 1 +
    
    // std::cout << "\n\ncount same cluster:\n" <<  count_same_cluster 
    //           << "\na0:\n" << a0 << "\nbn:\n" << b << "\n";
    
    phi = arma::randg( arma::distr_param(pred_ind + a0, 1/b) );
    // phi = arma::randg( arma::distr_param(count_same_cluster + a0, 1/b) );
    
  } else {
    phi = arma::randg( arma::distr_param(count_same_cluster + a0, 1/b) );
  }
    // phi = arma::randg( arma::distr_param(count_same_cluster + a0, 1/b) );
    // phi = arma::randg( arma::distr_param(pred_ind + a0, 1/b) );
  
  return phi;
}


// Calculates the normalising constant for the posterior
double calculate_normalising_constant(arma::vec cluster_weights_1,
                                      arma::vec cluster_weights_2,
                                      double context_similarity,
                                      arma::uword num_clusters_1,
                                      arma::uword num_clusters_2){
  double Z = 0.0;
  
  for(arma::uword i = 0; i < num_clusters_1; i++){
    Z += arma::sum(cluster_weights_1(i) * cluster_weights_2);
    if(i < num_clusters_2){
      Z += cluster_weights_1(i) * cluster_weights_2(i) * context_similarity;
    }
    // * (1 + context_similarity * (i == j));
    
    // for(arma::uword j = 0; j < num_clusters_2; j++){
      // Z += (cluster_weights_1(i) * cluster_weights_2(j))
            // * (1 + context_similarity * (i == j));
    // }
  }
  return Z;
}

// Sample the cluster membership of a categorical sample for MDI
arma::vec mdi_cat_clust_prob(arma::uword row_index,
                             arma::umat data,
                             arma::field<arma::mat> class_probs,
                             arma::uword num_clusters,
                             arma::uword num_cols_cat,
                             double phi,
                             arma::vec cluster_weights_categorical,
                             arma::uvec clust_labels,
                             arma::uvec clust_labels_comp){
  
  // cluster_labels_comparison is the labels of the data in the other context
  
  arma::vec prob_vec = arma::zeros<arma::vec>(num_clusters);
  arma::urowvec point = data.row(row_index);
  
  // pretty much this is the product of probabilities possibly up-weighted by
  // being in the same cluster in a different context and weighted by the cluster
  // weight in the current context
  double curr_weight = 0.0;
  
  // Upweight for similarity of contexts
  double similarity_upweight = 0.0;
  
  arma::uword common_cluster = 0;
  
  // std::cout << "Reached first for loop in cat function \n";
  
  for(arma::uword i = 0; i < num_clusters; i++){
    
    // calculate the log-weights for the context specific cluster and the across
    // context similarity
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
  
  // to handle overflowing
  prob_vec = exp(prob_vec - max(prob_vec));
  
  // normalise
  prob_vec = prob_vec / sum(prob_vec);
  
  return prob_vec;
}


// sample calculate probabilties for cluster allocation in the gaussian data in 
// MDI (hence the presence of the context similarity parameter)
arma::vec mdi_gauss_clust_probs(arma::uword row_index,
                                arma::mat data,
                                arma::uword k,
                                arma::cube mu,
                                arma::cube variance,
                                double context_similarity,
                                arma::vec cluster_weights,
                                arma::uvec cluster_labels,
                                arma::uvec cluster_labels_comp){
                                // bool outlier = false,
                                // arma::vec global_mean = arma::zeros<arma::vec>(1),
                                // arma::mat global_variance = arma::zeros<arma::mat>(1, 1),
                                // double t_df = 4.0){
  
  double curr_weight;
  // double exponent;
  double log_likelihood;
  
  // Upweight for similarity of contexts
  double similarity_upweight = 0.0;
  
  // double log_det;
  arma::vec prob_vec(k);
  
  arma::uvec count_probs;
  
  arma::uword d = data.n_cols;
  
  arma::uword common_cluster = 0;
  
  arma::vec point = arma::trans(data.row(row_index));
  
  common_cluster = 1 * (cluster_labels(row_index) == cluster_labels_comp(row_index));
  
  for(arma::uword i = 0; i < k ; i++){
    curr_weight = log(cluster_weights(i));
    
    // If in the cluster that the point is in in the comparison context, upweight
    if((i + 1) == cluster_labels(row_index)){
      similarity_upweight = log(1 + context_similarity * common_cluster);
    }
    
    // if(outlier && i == k){
    // 
    //   log_likelihood = t_likelihood(point, global_mean, global_variance, d, t_df);
    //                                
    // }
    // else {
      
    log_likelihood = normal_likelihood(point, mu.slice(i), variance.slice(i), d);

    // }
    
    prob_vec(i) = curr_weight + log_likelihood + similarity_upweight;
    
    similarity_upweight = 0.0;

  }
  
  // prob_vec = exp(prob_vec - max(prob_vec));
  // prob_vec = prob_vec / sum(prob_vec);
  
  return prob_vec;
}


// In a vector changes all values of ``label_1'' to ``label_2'' and vice versa
arma::uvec swap_labels(arma::uvec cluster_labels, 
                       arma::uword label_1, 
                       arma::uword label_2){
  
  arma::uvec label_1_ind = find(cluster_labels == label_1);
  arma::uvec label_2_ind = find(cluster_labels == label_2);
  
  cluster_labels.elem(label_1_ind).fill(label_2);
  cluster_labels.elem(label_2_ind).fill(label_1);
  return cluster_labels;
}

// Swap cluster weights
arma::vec swap_cluster_weights(arma::vec cluster_weights,
                               arma::uword label_1, 
                               arma::uword label_2){
  
  arma::vec decoy_weights = cluster_weights;
  
  cluster_weights(label_1) = decoy_weights(label_2);
  cluster_weights(label_2) = decoy_weights(label_1);
  
  return cluster_weights;
}

double comp(arma::uword n,
            double context_similarity,
            arma::uvec cluster_labels_1,
            arma::uvec cluster_labels_2){
  
  double score = 0.0;
  
  for(arma::uword i = 0; i < n; i++){
    score += log(1 + context_similarity * (cluster_labels_1(i) == cluster_labels_2(i)));
  }
  
  return score;
}

// The loglikelihood of the current clustering of MDI
// Eqn 4 of the supplement to the original MDI paper
double MDI_log_likelihood(double v,
                          arma::uword n,
                          double Z,
                          double phi,
                          arma::vec cluster_weights_1,
                          arma::vec cluster_weights_2,
                          arma::uvec cluster_labels_1,
                          arma::uvec cluster_labels_2){
  double log_likelihood = 0.0;
  
  // The coefficient outside the product
  log_likelihood += (n - 1) * log(v) - (v*Z) -log_factorial(n - 1);
  for(arma::uword i = 0; i < n; i++){
    // The observation specific update
    log_likelihood += log(cluster_weights_1(cluster_labels_1(i) - 1))
      + log(cluster_weights_2(cluster_labels_2(i) - 1))
      + log(1 + phi * (cluster_labels_1(i) == cluster_labels_2(i)));
  }
  
  return log_likelihood;
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

arma::vec cluster_label_update(arma::uvec cluster_labels_1,
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
  
  // double new_phi = 0.0;
  
  arma::uvec new_labels(n);
  arma::vec new_weights(num_clusters_2);
  
  double log_accept = 0.0;
  double accept = 0.0;
  
  double old_score = 0.0;
  double new_score = 0.0;
  
  // arma::uword curr_count = 0;
  // arma::uword new_count = 0;
  
  old_score = comp(n, phi, cluster_labels_1, cluster_labels_2);
  
  // arma::uvec occupied_clusters;
  // occupied_clusters = arma::unique(cluster_labels_1);
  
  // std::cout << occupied_clusters << "\n";
  
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
    
    // std::cout << "\nTruth statement: \n"
    //           << (arma::any(occupied_clusters) == i || arma::any(occupied_clusters) == new_pos) << "\n";
    
    // if(arma::any(occupied_clusters) == i || arma::any(occupied_clusters) == new_pos){
    // if(i < min_num_clusters || new_pos < min_num_clusters){
  
    new_labels = swap_labels(cluster_labels_2, i + 1, new_pos + 1);
    
    new_weights = swap_cluster_weights(cluster_weights_2, i, new_pos);
    
    new_score = comp(n, phi, cluster_labels_1, new_labels);
    
    log_accept = new_score - old_score;
    
    accept = 1.0;
    
    if(log_accept < 0){
      // std::cout << "\nOld score: " << old_score << "\nNew score: " << new_score << "\n";
      
      accept = exp(log_accept);
    }
    
    // curr_count = count_common_cluster(cluster_labels_1,
    //                                   cluster_labels_2,
    //                                   n);
    // 
    // new_count = count_common_cluster(cluster_labels_1,
    //                                  new_labels,
    //                                  n);

    // std::cout << "Current position: " << i << "\nNew position: " << new_pos
    //   << "\nCurrent phi: " << phi << "\nNew phi: " << new_phi
    //   << "\nCurrent common cluster count: " << curr_count
    //   << "\nNew common cluster count: " << new_count
    //   << "\n\nCluster 1 weights:\n" << cluster_weights_1
    //   << "\n\nCurrent weights:\n" << cluster_weights_2 << "\n\nNew weights:\n"
    //   << new_weights << "\n\n";
    
    if(arma::randu<double>( ) < accept){
      // std::cout << "\nOld score: " << old_score << "\nNew score: " << new_score << "\n";
      // std::cout << "Labels\n";
      cluster_labels_2 = new_labels;
      
      // std::cout << "Weights\n";
      cluster_weights_2 = new_weights;
      
      old_score = comp(n, phi, cluster_labels_1, cluster_labels_2);
  
    }
    // }
  }
  
  arma::vec output = arma::zeros<arma::vec>(n + num_clusters_2 + 1);
  output.subvec(0, n - 1) = arma::conv_to<arma::vec>::from(cluster_labels_2);
  
  output.subvec(n, n + num_clusters_2 - 1) = cluster_weights_2;
  output(n + num_clusters_2) = phi;
  return output;
}

// Declare object to save class_probabilities in for categorical clustering
arma::field<arma::field<arma::mat>> class_probs_obj(arma::uword n_clust,
                                                    arma::uword n_col,
                                                    arma::uword n_iter,
                                                    arma::uvec cat_count){
  arma::field<arma::field<arma::mat>> class_probabilities_saved(n_clust);
  arma::field<arma::rowvec> comp_class_probs(n_col);
  
  for(arma::uword i = 0; i < n_clust; i++){
    arma::field<arma::mat> loc_level(n_col);
    
    for(arma::uword j = 0; j < n_col; j++){
      arma::mat iter_level(n_iter, cat_count(j));
      iter_level.zeros();
      loc_level(j) = iter_level;
    }
    class_probabilities_saved(i) = loc_level;
  }
  return class_probabilities_saved;
}
// 
// // MDI clustering for a gaussian and cateogrical dataset
// // [[Rcpp::export]]
// Rcpp::List mdi_gauss_cat(arma::mat gaussian_data,
//                          arma::umat categorical_data,
//                          arma::vec mu_0,
//                          double lambda_0,
//                          arma::mat scale_0,
//                          int df_0,
//                          double a0,
//                          double b0,
//                          arma::vec cluster_weight_priors_gaussian,
//                          arma::vec cluster_weight_priors_categorical,
//                          arma::field<arma::vec> phi_prior,
//                          arma::uvec cluster_labels_gaussian,
//                          arma::uvec cluster_labels_categorical,
//                          arma::uword num_clusters_gaussian,
//                          arma::uword num_clusters_categorical,
//                          arma::uvec fix_vec_1,
//                          arma::uvec fix_vec_2,
//                          arma::uword num_iter,
//                          arma::uword burn,
//                          arma::uword thinning,
//                          bool outlier = false,
//                          double t_df = 4.0,
//                          bool record_posteriors = false,
//                          bool normalise = false,
//                          double u_1 = 2,
//                          double v_1 = 10
// ){
//   
//   // Declare the sample size and dimensionality of the continuous and 
//   // categorical data
//   
//   // std::cout << "In function \n";
//   arma::uword n = gaussian_data.n_rows;
//   arma::uword num_cols_cont = gaussian_data.n_cols;
//   
//   // arma::uword n_cat = categorical_data.n_rows;
//   arma::uword num_cols_cat = categorical_data.n_cols;
//   
// 
//   // Frequently will compare clusters across contexts so this is a useful bound
//   // to iterate to
//   arma::uword min_num_clusters = std::min(num_clusters_gaussian,
//                                           num_clusters_categorical);
// 
//   //  Normalise the continuous data
//   if(normalise){
//     gaussian_data = arma::normalise(gaussian_data);
//   }
// 
//   // Declare global variance and mean - used in outlier t-distribution 
//   // (if included)
//   arma::mat global_variance(num_cols_cont, num_cols_cont);
//   global_variance = 0.5 * arma::cov(gaussian_data); // Olly's rec
//   
//   arma::vec global_mean(num_cols_cont);
//   global_mean = arma::trans(arma::mean(gaussian_data, 0));
//   
//   
//   double v = 0.0; // strategic latent variable
//   
//   // Cluster weights for each dataset
//   arma::vec cluster_weights_gaussian(num_clusters_gaussian);
//   arma::vec cluster_weights_categorical(num_clusters_categorical);
//   
//   // These are the local cubes of posterior mean and variance overwritten each
//   // iteration
//   arma::field<arma::cube> loc_mu_variance(2);
//   
//     // std::cout << "Declared to mean/variance thing \n";
//   
//   // Declare the field for the phi variable for the categorical data
//   arma::uvec cat_count(num_cols_cat);
//   cat_count = cat_counter(categorical_data);
//   arma::field<arma::mat> class_probabilities(num_cols_cat);
//   
//   class_probabilities = declare_class_probs_field(cat_count,
//                                                   num_cols_cat,
//                                                   num_clusters_categorical);
//   
//   // Initialise the context similarity parameter based on priors
//   double context_similarity = arma::randg(arma::distr_param(a0, 1/b0) );
// 
//   // Declare the normalising constant
//   double Z = 1.0;
//  
//   arma::vec curr_gaussian_prob_vec(num_clusters_gaussian);
//   arma::vec curr_categorical_prob_vec(num_clusters_categorical);
// 
//   // Various objects to record values for posterior distributions and clustering
//   // the record for similarity in each clustering
//   arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
//   
//   // These are the lists recording the posterior mean and 
//   // variance for each class for each recorded iteration
//   ListMatrix variance(eff_count, num_clusters_gaussian);
//   ListMatrix mu(eff_count, num_clusters_gaussian);
//   
//   // Objects required to save the categorical variable for each component
//   arma::field<arma::field<arma::mat>> class_probabilities_saved(num_clusters_categorical);
//   // arma::field<arma::rowvec> comp_class_probs(num_cols_cat);
//   
//   class_probabilities_saved = class_probs_obj(num_clusters_categorical,
//                                               num_cols_cat,
//                                               eff_count,
//                                               cat_count);
//   
//   // Records of allocation
//   arma::umat gaussian_record(n, eff_count);
//   gaussian_record.zeros();
//   
//   arma::umat categorical_record(n, eff_count);
//   categorical_record.zeros();
//   
//   arma::vec context_similarity_record(eff_count);
//   
//   // Filenames for saving posteriors
//   // std::string gauss_lab_folder = "gaussian_allocation";
//   // std::string cat_lab_folder = "categorical_allocation/";
//   // std::string out_lab_folder = "outlier_allocation/";
//   
//   std::string gauss_lab_file = "gaussian_allocation/gaussian_labels_iter_";
//   std::string cat_lab_file = "categorical_allocation/categorical_labels_iter_";
//   std::string out_lab_file = "outlier_allocation/outlier_labels_iter_";
//   std::string i_str;
// 
//   std::string gauss_lab_file_loc;
//   std::string cat_lab_file_loc;
//   std::string out_lab_file_loc;
//   
//   arma::uword record_iter = 0;
//   
//   // To hold output of label flipping function
//   arma::vec labels_weights_phi(n + num_clusters_categorical + 1);
//   
//   arma::vec entropy_cw(num_iter);
// 
//   // Rate priors for weight smapling from Gamma distn
//   arma::vec rate_0_gauss(num_clusters_gaussian);
//   arma::vec rate_0_cat(num_clusters_categorical);
//   
//   // Placeholder prior
//   rate_0_gauss.fill(1);
//   rate_0_cat.fill(1);
//   
//   // Not sure this is sensible
//   cluster_weights_gaussian = cluster_weight_priors_gaussian;
//   cluster_weights_categorical = cluster_weight_priors_categorical;
//   
//   // OUTLIER COMPONENT
//   // Variables to handle outlier component from tagm
//   // Vector of 0 and 1 for points assigned to outlier group or not
//   arma::uvec outlier_vec(n);
//   
//   
//   arma::uvec b_k(n); // the vector of b_k's the sum of which gives b
//   double b = 0.0; // a will be N - b, no need to declare
//   
//   // for recording the probabilities for each class
//   arma::cube gaussian_class_probs(eff_count, num_clusters_gaussian, n);
//   arma::cube cat_class_probs(eff_count, num_clusters_categorical, n); 
//   
//   // Begin with all non-fixed points (i.e. unlabelled) to outlier component
//   outlier_vec.fill(0);
//   if(outlier && arma::any(fix_vec_1)){
//     outlier_vec = 1 - fix_vec_1;
//   }
//   
//   // the current iterations probabilities (overwritten - saved to the above cube
//   // after burn in and as defined by thinning)
//   arma::vec curr_class_probs(num_clusters_gaussian);
//   arma::vec curr_norm_likelihoods(num_clusters_gaussian);
//   
//   // Class labels of points not currently assigned as outliers
//   arma::uvec relevant_labels(n);
//   
//   // the current iterations probabilities of being an outlier (non-outlier prob
//   // is 1 - curr_outlier_prob)
//   arma::vec curr_outlier_prob(2);
//   double curr_outlier_likelihood = 0.0;
//   arma::uword predicted_outlier = 0;
//   double outlier_weight = 1 - sample_beta(u_1, v_1);
// 
//   // the predicted class assuming the point is not an outlier
//   arma::uword predicted_class = 0;
//   
//   // This is where we save the outlier labels
//   arma::umat outlier_probs_saved(n, eff_count);
//   
//   // Calculate the current normalising constant (consider being more clever 
//   // about this) 
//   // Z = calculate_normalising_constant(cluster_weights_gaussian,
//   //                                    cluster_weights_categorical,
//   //                                    context_similarity,
//   //                                    num_clusters_gaussian,
//   //                                    num_clusters_categorical);
//   
//   double predicted_norm_likelihood = 0.0;
//   
//   for(arma::uword i = 0; i < num_iter + 1; i++){
//     
//     // Consider only the labels of points not considered to be outliers
//     relevant_labels = cluster_labels_gaussian % (1 - outlier_vec);
// 
//     // sample the strategic latent variable, v
//     v = arma::randg( arma::distr_param(n, 1/Z) );
//     
//     // Entropy for graphing convergence
//     entropy_cw(i) = entropy(cluster_weights_gaussian);
//     
//     // Sample the posterior mean and variance for the gaussian data
//     loc_mu_variance = mean_variance_sampling(gaussian_data,
//                                              relevant_labels,
//                                              num_clusters_gaussian,
//                                              df_0,
//                                              num_cols_cont,
//                                              scale_0,
//                                              lambda_0,
//                                              mu_0);
// 
//     // For the categorical data, sample the probabilities for each class
//     class_probabilities = sample_class_probabilities(categorical_data,
//                                                      class_probabilities,
//                                                      phi_prior,
//                                                      cluster_labels_categorical,
//                                                      cat_count,
//                                                      num_clusters_categorical,
//                                                      num_cols_cat);
//   
//     // sample cluster weights for the two datasets]
//     // Gaussian weights
//     cluster_weights_gaussian = mdi_cluster_weights(cluster_weight_priors_gaussian,
//                                                    rate_0_gauss,
//                                                    v,
//                                                    num_clusters_gaussian,
//                                                    num_clusters_categorical,
//                                                    cluster_weights_categorical,
//                                                    cluster_labels_gaussian,
//                                                    // relevant_labels,
//                                                    cluster_labels_categorical,
//                                                    context_similarity);
// 
//     // Categorical weights
//     cluster_weights_categorical = mdi_cluster_weights(cluster_weight_priors_categorical,
//                                                       rate_0_cat,
//                                                       v,
//                                                       num_clusters_categorical,
//                                                       num_clusters_gaussian,
//                                                       cluster_weights_gaussian,
//                                                       cluster_labels_categorical,
//                                                       cluster_labels_gaussian,
//                                                       context_similarity);
//     
//     // Calculate the current normalising constant (consider being more clever
//     // about this)
//     Z = calculate_normalising_constant(cluster_weights_gaussian,
//                                        cluster_weights_categorical,
//                                        context_similarity,
//                                        num_clusters_gaussian,
//                                        num_clusters_categorical);
// 
//     // sample the context similarity parameter (as only two contexts this is a
//     // single double - number not a drink)
//     context_similarity = sample_phi(cluster_labels_gaussian,
//                                     cluster_labels_categorical,
//                                     cluster_weights_gaussian,
//                                     cluster_weights_categorical,
//                                     v,
//                                     n,
//                                     min_num_clusters,
//                                     a0,
//                                     b0);
// 
//     // sample class for each observation
//     for(arma::uword j = 0; j < n; j++){
//       
//       // for each point create the vector of probabilities associated with 
//       // assignment to each cluster
//       curr_norm_likelihoods = mdi_gauss_clust_probs(j,
//                                                     gaussian_data,
//                                                     num_clusters_gaussian,
//                                                     loc_mu_variance(1),
//                                                     loc_mu_variance(0),
//                                                     context_similarity,
//                                                     cluster_weights_gaussian,
//                                                     relevant_labels,
//                                                     cluster_labels_categorical);
//       
//       curr_gaussian_prob_vec = over_flow_handling(curr_norm_likelihoods);
//       predicted_class = cluster_predictor(curr_gaussian_prob_vec);
//       
//       // The various probabilities to determine if the observation is considered 
//       // an outlier or not
//       if(outlier){
//         
//         curr_outlier_likelihood = sample_outlier(arma::trans(gaussian_data.row(j)),
//                                                  gaussian_data,
//                                                  outlier_weight,
//                                                  global_mean,
//                                                  global_variance,
//                                                  t_df,
//                                                  u_1 + b,
//                                                  v_1 + n - b);
//         
//         predicted_norm_likelihood = normal_likelihood(arma::trans(gaussian_data.row(j)),
//                                                       loc_mu_variance(1).slice(predicted_class - 1),
//                                                       loc_mu_variance(0).slice(predicted_class - 1),
//                                                       num_cols_cont);
//         
//         predicted_norm_likelihood += log(1 - outlier_weight);
//         curr_outlier_prob(0) = predicted_norm_likelihood;
//         
//         // Overflow handling
//         curr_outlier_prob = exp(curr_outlier_prob - max(curr_outlier_prob));
//         
//         // Normalise the vector
//         curr_outlier_prob = curr_outlier_prob / sum(curr_outlier_prob);
//         
//         predicted_outlier = cluster_predictor(curr_outlier_prob) - 1; // as +1 to handle R
//       }
//       
//       curr_categorical_prob_vec = mdi_cat_clust_prob(j,
//                                                      categorical_data,
//                                                      class_probabilities,
//                                                      num_clusters_categorical,
//                                                      num_cols_cat,
//                                                      context_similarity,
//                                                      cluster_weights_categorical,
//                                                      cluster_labels_gaussian,
//                                                      cluster_labels_categorical);
//       
//       
//       
//       if (i >= burn && (i - burn) % thinning == 0 && record_posteriors) {
//         record_iter = (i - burn) / thinning;
//         gaussian_class_probs.slice(j).row(record_iter) = arma::trans(curr_gaussian_prob_vec);
//         cat_class_probs.slice(j).row(record_iter) = arma::trans(curr_categorical_prob_vec);
//       }
// 
//       // update labels - in gaussian data this is only if the current point is 
//       // not fixed
//       if(fix_vec_1[j] == 0){
//         cluster_labels_gaussian(j) = predicted_class; // cluster_predictor(curr_gaussian_prob_vec);
//         if(outlier){
//           outlier_vec(j) = predicted_outlier;
//         }
//       }
//       
//       if (fix_vec_2[j] == 0){
//         cluster_labels_categorical(j) = cluster_predictor(curr_categorical_prob_vec);
//       }
//       
//     }
// 
//     if(outlier){
//       
//       // Components of outlier weight
//       b_k = arma::find(outlier_vec == 0);
//       b = b_k.n_elem;
//       
//       // Sample outlier weight
//       outlier_weight = 1 - sample_beta(u_1 + b, v_1 + n - b);
// 
//     }
// 
//     // Update cluster labels in second dataset
//     // Return the new labels, weights and similarity in a single vector
//     labels_weights_phi = cluster_label_update(cluster_labels_gaussian,
//                                               cluster_labels_categorical,
//                                               cluster_weights_gaussian,
//                                               cluster_weights_categorical,
//                                               num_clusters_gaussian,
//                                               num_clusters_categorical,
//                                               context_similarity,
//                                               min_num_clusters,
//                                               v,
//                                               n,
//                                               a0,
//                                               b0,
//                                               Z);
// 
//     // Separate the output into the relevant components
//     cluster_labels_categorical = arma::conv_to<arma::uvec>::from(labels_weights_phi.subvec(0, n - 1));
//     cluster_weights_categorical = labels_weights_phi.subvec(n, n + num_clusters_categorical - 1);
//     context_similarity = arma::as_scalar(labels_weights_phi(n + num_clusters_categorical));
// 
//     // if current iteration is a recorded iteration, save the labels
//     if (i >= burn && (i - burn) % thinning == 0) {
//       
//       record_iter = (i - burn) / thinning;
//       
//       // gaussian_record.col(record_iter) = cluster_labels_gaussian;
//       // categorical_record.col(record_iter) = cluster_labels_categorical;
//       context_similarity_record(record_iter) = context_similarity;
//       // outlier_probs_saved.col(record_iter) = outlier_vec;
//       
//       
//       // Save to file
//       i_str = std::to_string(record_iter);
// 
//       gauss_lab_file_loc = gauss_lab_file + i_str;
//       cat_lab_file_loc = cat_lab_file + i_str;
//       out_lab_file_loc = out_lab_file + i_str;
// 
//       cluster_labels_gaussian.save(gauss_lab_file_loc);
//       cluster_labels_categorical.save(cat_lab_file_loc);
//       outlier_vec.save(out_lab_file_loc);
//       
//       // Record posteriors of parameters for Gaussian and Categorical
//       // distributions
//       if(record_posteriors){
//         for(arma::uword j = 0; j < num_clusters_gaussian; j++){
//           mu(record_iter, j) = loc_mu_variance(1).slice(j);
//           variance(record_iter, j) = loc_mu_variance(0).slice(j);
//         }
//         
//         for(arma::uword j = 0; j < num_clusters_categorical; j++){
//           for(arma::uword k = 0; k < num_cols_cat; k++){
//             class_probabilities_saved(j)(k).row(record_iter) = class_probabilities(k).row(j);
//           }
//         }
//       }
//     }
//   
//   }
//   
//   // Loading posterior objects
//   for(arma::uword i = 0; i < eff_count; i++){
//     i_str = std::to_string(record_iter);
//     
//     gauss_lab_file_loc = gauss_lab_file + i_str;
//     cat_lab_file_loc = cat_lab_file + i_str;
//     out_lab_file_loc = out_lab_file + i_str;
//     
//     cluster_labels_gaussian.load(gauss_lab_file_loc);
//     cluster_labels_categorical.load(cat_lab_file_loc);
//     outlier_vec.load(out_lab_file_loc);
//     
//     gaussian_record.col(i) = cluster_labels_gaussian;
//     categorical_record.col(i) = cluster_labels_categorical;
//     outlier_probs_saved.col(i) = outlier_vec;
// 
//     // if(record_posteriors){
//     //   for(arma::uword j = 0; j < num_clusters_gaussian; j++){
//     //     mu(record_iter, j) = loc_mu_variance(1).slice(j);
//     //     variance(record_iter, j) = loc_mu_variance(0).slice(j);
//     //   }
//     //   
//     //   for(arma::uword j = 0; j < num_clusters_categorical; j++){
//     //     for(arma::uword k = 0; k < num_cols_cat; k++){
//     //       class_probabilities_saved(j)(k).row(record_iter) = class_probabilities(k).row(j);
//     //     }
//     //   }
//     // }
//     
//   }
//   
//   // construct similarity matrix
//   arma::mat sim(n, n); 
//   arma::mat cat_sim(n, n);
//   sim = similarity_mat(gaussian_record);
//   cat_sim = similarity_mat(categorical_record);
//   
//   return List::create(Named("similarity_1") = sim,
//                       Named("similarity_2") = cat_sim,
//                       Named("class_record_1") = gaussian_record,
//                       Named("class_record_2") = categorical_record,
//                       Named("mean_posterior") = mu,
//                       Named("variance_posterior") = variance,
//                       Named("class_prob_posterior") = class_probabilities_saved,
//                       Named("context_similarity") = context_similarity_record,
//                       Named("entropy") = entropy_cw,
//                       Named("outlier") = outlier_probs_saved);
// }
// 
// 
// ///////////////////////////////////////////////////////////////////////////////
// 
// // MDI for different types
// 
// ///////////////////////////////////////////////////////////////////////////////
// 
// // MDI clustering for two gaussian datasets
// // [[Rcpp::export]]
// Rcpp::List mdi_gauss_gauss(arma::mat data_1,
//                            arma::mat data_2,
//                            arma::vec mu_0_1,
//                            double lambda_0_1,
//                            arma::mat scale_0_1,
//                            int df_0_1,
//                            arma::vec mu_0_2,
//                            double lambda_0_2,
//                            arma::mat scale_0_2,
//                            int df_0_2,
//                            arma::vec clust_weight_priors_1,
//                            arma::vec clust_weight_priors_2,
//                            arma::uvec clust_labels_1,
//                            arma::uvec clust_labels_2,
//                            arma::uword n_clust_1,
//                            arma::uword n_clust_2,
//                            arma::uvec fix_vec_1,
//                            arma::uvec fix_vec_2,
//                            double a0,
//                            double b0,
//                            arma::uword num_iter,
//                            arma::uword burn,
//                            arma::uword thinning,
//                            bool outlier_1 = false,
//                            double t_df_1 = 4.0,
//                            bool outlier_2 = false,
//                            double t_df_2 = 4.0,
//                            bool record_posteriors = false,
//                            bool normalise_1 = false,
//                            bool normalise_2 = false,
//                            double u_1 = 2,
//                            double v_1 = 10,
//                            double u_2 = 2,
//                            double v_2 = 10
// ){
//   
//   // Declare the sample size and dimensionality of the continuous and 
//   // categorical data
//   
//   // std::cout << "In function \n";
//   
//   arma::uword n = data_1.n_rows;
//   arma::uword n_cols_1 = data_1.n_cols;
//   
//   // arma::uword n_cat = categorical_data.n_rows;
//   arma::uword n_cols_2 = data_2.n_cols;
//   
//   // Frequently will compare clusters across contexts so this is a useful bound
//   // to iterate to
//   arma::uword min_n_clust = std::min(n_clust_1, n_clust_2);
//   
//   // std::cout << "Nomalising data\n";
//   
//   //  Normalise the continuous data
//   if(normalise_1){
//     data_1 = arma::normalise(data_1);
//   }
//   
//   if(normalise_2){
//     data_2 = arma::normalise(data_2);
//   }
//   
//   // Declare global variance and mean - used in outlier t-distribution 
//   // (if included)
//   arma::mat global_variance_1(n_cols_1, n_cols_1);
//   global_variance_1 = 0.5 * arma::cov(data_1); // Olly's rec
//   
//   arma::vec global_mean_1(n_cols_1);
//   global_mean_1 = arma::trans(arma::mean(data_1, 0));
//   
//   arma::mat global_variance_2(n_cols_2, n_cols_2);
//   global_variance_2 = 0.5 * arma::cov(data_2); // Olly's rec
//   
//   arma::vec global_mean_2(n_cols_2);
//   global_mean_2 = arma::trans(arma::mean(data_2, 0));
//   
//   
//   double v = 0.0; // strategic latent variable
//   
//   // Cluster weights for each dataset
//   arma::vec clust_weights_1(n_clust_1);
//   arma::vec clust_weights_2(n_clust_2);
//   
//   // These are the local cubes of posterior mean and variance overwritten each
//   // iteration
//   arma::field<arma::cube> loc_mu_variance_1(2);
//   arma::field<arma::cube> loc_mu_variance_2(2);
//   
//   // std::cout << "Generic message. Parameter fields declared.\n";
//   
//   // Context similarity - smaple prior
//   double phi = arma::randg(arma::distr_param(a0, 1/b0) );
//   
//   // Declare normalising constant
//   double Z = 1.0;
//   
//   // Used in each iteration
//   arma::vec curr_prob_vec_1(n_clust_1);
//   arma::vec curr_prob_vec_2(n_clust_2);
//   
//   // std::cout << "Declared prob vectors \n";
//   
//   // the record for similarity in each clustering
//   arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
//   
//   arma::umat record_1(n, eff_count);
//   record_1.zeros();
//   
//   arma::umat record_2(n, eff_count);
//   record_2.zeros();
//   
//   // The vector to record the context similarity parameter
//   arma::vec phi_record(eff_count);
//   
//   arma::vec labels_weights_phi(n + n_clust_2 + 1);
//   
//   arma::vec entropy_cw(num_iter);
//   
//   // std::cout << "Entropy vec declared.\n";
//   
//   arma::vec rate_0_1(n_clust_1);
//   arma::vec rate_0_2(n_clust_2);
//   
//   // Placeholder prior
//   rate_0_1.fill(1);
//   rate_0_2.fill(1);
//   
//   // std::cout << "Rates rated.\n";
//   
//   // Not sure this is sensible - simply to have non-zero numbers here
//   clust_weights_1 = clust_weight_priors_1;
//   clust_weights_2 = clust_weight_priors_2;
//   
//   // OUTLIER COMPONENT
//   // Variables to handle outlier component from tagm
//   // Vector of 0 and 1 for points assigned to outlier group or not
//   arma::uvec outlier_vec_1(n);
//   arma::uvec outlier_vec_2(n);
//   
//   outlier_vec_1.fill(0);
//   outlier_vec_2.fill(0);
//   
//   // Begin with all non-fixed points (i.e. unlabelled) to outlier component
//   if(outlier_1 && arma::any(fix_vec_1)){
//     outlier_vec_1 = 1 - fix_vec_1;
//   }
//   
//   if(outlier_2 && arma::any(fix_vec_2)){
//     outlier_vec_2 = 1 - fix_vec_2;
//   }
//   
//   // Declare object for counting the number of items in the outlier component
//   arma::uvec b_k_1(n);
//   double b_1 = 0.0; 
//   
//   arma::uvec b_k_2(n);
//   double b_2 = 0.0;
//   
//   // for recording the probabilities for each class
//   arma::cube class_probs_1(eff_count, n_clust_1, n);
//   arma::cube class_probs_2(eff_count, n_clust_2, n); 
//   
//   // the current iterations probabilities (overwritten - saved to the above cube
//   // after burn in and as defined by thinning)
//   arma::vec curr_class_probs_1(n_clust_1);
//   arma::vec norm_likelihoods_1(n_clust_1);
//   
//   arma::vec curr_class_probs_2(n_clust_2);
//   arma::vec norm_likelihoods_2(n_clust_2);
//   
//   // Class labels of points not currently assigned as outliers
//   arma::uvec relevant_labels_1(n);
//   arma::uvec relevant_labels_2(n);
//   
//   // the current iterations probabilities of being an outlier (non-outlier prob
//   // is 1 - curr_outlier_prob)
//   arma::vec curr_outlier_prob_1(2);
//   double outlier_likelihood_1 = 0.0;
//   arma::uword predicted_outlier_1 = 0;
//   double outlier_weight_1 = 1 - sample_beta(u_1, v_1);
//   
//   arma::vec curr_outlier_prob_2(2);
//   double outlier_likelihood_2 = 0.0;
//   arma::uword predicted_outlier_2 = 0;
//   double outlier_weight_2 = 1 - sample_beta(u_2, v_2);
//   
//   // the predicted class assuming the point is not an outlier for the two contexts
//   arma::uword predicted_class_1 = 0;
//   arma::uword predicted_class_2 = 0;
//   
//   // This is where we save the outlier labels
//   arma::umat outlier_probs_saved_1(n, eff_count);
//   arma::umat outlier_probs_saved_2(n, eff_count);
//   
//   
//   // std::cout << "All declared \n";
//   
//   // Z = calculate_normalising_constant(clust_weights_1,
//   //                                    clust_weights_2,
//   //                                    phi,
//   //                                    n_clust_1,
//   //                                    n_clust_2);
//   
//   for(arma::uword i = 0; i < num_iter; i++){
//     
//     // sample cluster weights for the two datasets
//     // clust_weights_1 = gamma_posterior(clust_weight_priors_1,
//     //                                   clust_labels_1,
//     //                                   n_clust_1);
//     
//     
//     // Consider only the labels of points not considered to be outliers
//     relevant_labels_1 = clust_labels_1 % (1 - outlier_vec_1);
//     relevant_labels_2 = clust_labels_2 % (1 - outlier_vec_2);
//     
//     // std::cout << "\nStrategic latent variable sampling\n";
//     
//     // sample the strategic latent variable, v
//     v = arma::randg( arma::distr_param(n, 1/Z) );
//     
//     arma::uvec wert = find(outlier_vec_1);
//     arma::uvec wart = find(outlier_vec_2);
//     
//     // std::cout << "\nStrategic latent variable sampled: " << v << "\n";
//     // std::cout << "Context similarity parameter: " << phi << "\n";
//     // std::cout << "Normalising constant: " << Z << "\n";
//     // std::cout << "\nNum Relevant labels (context 1): " << n - wert.n_elem << "\n";
//     // std::cout << "b_1: " << b_1 << "\n";
//     // 
//     // std::cout << "\nNum Relevant labels (context 2): " << n - wart.n_elem << "\n";
//     // std::cout << "b_2: " << b_2 << "\n";
//   
//     
//     // Entropy for graphing convergence
//     entropy_cw(i) = entropy(clust_weights_1);
//     
//     // Sample the posterior mean and variance for the first dataset
//     loc_mu_variance_1 = mean_variance_sampling(data_1,
//                                                relevant_labels_1,
//                                                n_clust_1,
//                                                df_0_1,
//                                                n_cols_1,
//                                                scale_0_1,
//                                                lambda_0_1,
//                                                mu_0_1);
//     
//     // std::cout << "Variance sampled\n";
//     
//     // Sample the posterior mean and variance for the second dataset
//     loc_mu_variance_2 = mean_variance_sampling(data_2,
//                                                relevant_labels_2,
//                                                n_clust_2,
//                                                df_0_2,
//                                                n_cols_2,
//                                                scale_0_2,
//                                                lambda_0_2,
//                                                mu_0_2);
//     
//     // std::cout << "Sampled parameters for both datasets\n";
//     
//     // Sample cluster weights within each context
//     clust_weights_1 = mdi_cluster_weights(clust_weight_priors_1,
//                                           rate_0_1,
//                                           v,
//                                           n_clust_1,
//                                           n_clust_2,
//                                           clust_weights_2,
//                                           relevant_labels_1,
//                                           relevant_labels_2,
//                                           phi);
//     
//     // std::cout << "\nCluster weights 1:\n" << clust_weights_1 << "\n";
// 
//     
//     clust_weights_2 = mdi_cluster_weights(clust_weight_priors_2,
//                                           rate_0_2,
//                                           v,
//                                           n_clust_2,
//                                           n_clust_1,
//                                           clust_weights_1,
//                                           relevant_labels_2,
//                                           relevant_labels_1,
//                                           phi);
//     
//     // std::cout << "\nCluster weights 2:\n" << clust_weights_2 << "\n";
//     
//     // std::cout << "Sampled weights for both datasets\n";
//     
//     
//     
//     // sample the context similarity parameter (as only two contexts this is a
//     // single double - number not a drink)
//     phi = sample_phi(clust_labels_1,
//                      clust_labels_2,
//                      clust_weights_1,
//                      clust_weights_2,
//                      v,
//                      n,
//                      min_n_clust,
//                      a0,
//                      b0);
//     
//     // std::cout << "Sampled phi\n";
//     
//     // Calculate the current normalising constant (consider being more clever 
//     // about this) 
//     Z = calculate_normalising_constant(clust_weights_1,
//                                        clust_weights_2,
//                                        phi,
//                                        n_clust_1,
//                                        n_clust_2);
//     
//     // std::cout << "Z calculated \n";
//     
//     // sample the strategic latent variable, v
//     // v = arma::randg( arma::distr_param(n, 1/Z) );
//     
//     // sample 
//     for(arma::uword j = 0; j < n; j++){
//       
//       // Calculate the log-likelihoods for each cluster in each context
//       norm_likelihoods_1 = mdi_gauss_clust_probs(j,
//                                                  data_1,
//                                                  n_clust_1,
//                                                  loc_mu_variance_1(1),
//                                                  loc_mu_variance_1(0),
//                                                  phi,
//                                                  clust_weights_1,
//                                                  relevant_labels_1,
//                                                  relevant_labels_2);
//       
//       norm_likelihoods_2 = mdi_gauss_clust_probs(j,
//                                                  data_2,
//                                                  n_clust_2,
//                                                  loc_mu_variance_2(1),
//                                                  loc_mu_variance_2(0),
//                                                  phi,
//                                                  clust_weights_2,
//                                                  relevant_labels_2,
//                                                  relevant_labels_1);
//       
//       // std::cout << "Normal likelihoods calculated.\n";
//       
//       // Convert to likelihoods, handle overflow and normalise
//       curr_prob_vec_1 = over_flow_handling(norm_likelihoods_1);
//       curr_prob_vec_2 = over_flow_handling(norm_likelihoods_2);
//       
//       // Predict the component to which the obersation belongs
//       predicted_class_1 = cluster_predictor(curr_prob_vec_1);
//       predicted_class_2 = cluster_predictor(curr_prob_vec_2);
//       
//       // std::cout << "Classes predicted.\n";
//       
//       if(outlier_1){
//         
//         // std::cout << "Outlier 1.\n";
//         
//         // Calculate the likelihood associated with the outlier component
//         outlier_likelihood_1 = sample_outlier(arma::trans(data_1.row(j)),
//                                               data_1,
//                                               outlier_weight_1,
//                                               global_mean_1,
//                                               global_variance_1,
//                                               t_df_1,
//                                               u_1 + b_1,
//                                               v_1 + n - b_1);
//         
//         // std::cout << "Outlier likelihood calculated.\n";
//         
//         // Put it into a vector (for normalising and whatnot)
//         curr_outlier_prob_1(1) = outlier_likelihood_1;
//         
//         // std::cout << "t likelihood OK!\n";
//         
//         // Put in the likelihood of the predicted class, with the class weight
//         // dropped and the non-outlier weight added
//         curr_outlier_prob_1(0)= norm_likelihoods_1(predicted_class_1 - 1) 
//           + log(1 - outlier_weight_1)
//           - log(clust_weights_1(predicted_class_1 - 1));
//           
//         // Convert to likelihoods, handle overflow and normalise
//         curr_outlier_prob_1 = over_flow_handling(curr_outlier_prob_1);
//         
//         
//         // std::cout << "\nOutlier prob (context 1):\n" << curr_outlier_prob_1 << "\n";
//         
//         // Predict if the current observation is an outlier or not  
//         predicted_outlier_1 = cluster_predictor(curr_outlier_prob_1) - 1; // as +1 to handle R
//         
//         // std::cout << "Outlier prediction done.\n";
//       }
//       
//       // Similarly for context 2
//       if(outlier_2){
//         
//         // std::cout << "In outlier 2 caculations.\n";
//         
//         // Calculate the likelihood associated with the outlier component
//         outlier_likelihood_2 = sample_outlier(arma::trans(data_2.row(j)),
//                                               data_2,
//                                               outlier_weight_2,
//                                               global_mean_2,
//                                               global_variance_2,
//                                               t_df_2,
//                                               u_2 + b_2,
//                                               v_2 + n - b_2);
//         
//         // std::cout << "Outlier likelihood calculated.\n";
//         
//         // Put it into a vector (for normalising and whatnot)
//         curr_outlier_prob_2(1) = outlier_likelihood_2;
//         
//         // std::cout << "t likelihood OK!\n";
//         
//         // Put in the likelihood of the predicted class, with the class weight
//         // dropped and the non-outlier weight added
//         curr_outlier_prob_2(0)= norm_likelihoods_2(predicted_class_2 - 1) 
//           + log(1 - outlier_weight_2)
//           - log(clust_weights_2(predicted_class_2 - 1));
//           
//         // Convert to likelihoods, handle overflow and normalise
//         curr_outlier_prob_2 = over_flow_handling(curr_outlier_prob_2);
//         
//         // std::cout << "Outlier prob (context 2):\n" << curr_outlier_prob_2 << "\n";
//         
//         // Predict if the current observation is an outlier or not  
//         predicted_outlier_2 = cluster_predictor(curr_outlier_prob_2) - 1; // as +1 to handle R
//         
//         // std::cout << "Outlier prediction done.\n";
//       }
// 
//       // std::cout << "\nContext 1 prediction: " << predicted_class_1
//       //           << "\nContext 2 prediction: " << predicted_class_2 << "\n";
//       
//       // update labels if current point's label is not fixed
//       if(fix_vec_1[j] == 0){
//         // std::cout << "Gaussian 1\n";
//         clust_labels_1(j) = predicted_class_1;
//         if(outlier_1){
//           // std::cout << "Gaussian 1 outlier\n";
//           outlier_vec_1(j) = predicted_outlier_1;
//         }
//       }
//       
//       if(fix_vec_2[j] == 0){
//         // std::cout << "Gaussian 2\n";
//         clust_labels_2(j) = predicted_class_2;
//         if(outlier_2){
//           // std::cout << "Gaussian 2 outlier\n";
//           outlier_vec_2(j) = predicted_outlier_2;
//         }
//       }
//       
//     }
//     
//     // std::cout << "All the context update stuff\n";
//     
//     // std::cout << cluster_labels_categorical.n_elem << "\n"; 
//     
//     // Update cluster labels in second dataset
//     // Return the new labels, weights and similarity in a single vector
//     
//     // Do not allow label flipping if any of context 2 have fixed labels
//     
//     // if(data_2_unsupervised){
//     labels_weights_phi = cluster_label_update(clust_labels_1,
//                                               clust_labels_2,
//                                               clust_weights_1,
//                                               clust_weights_2,
//                                               n_clust_1,
//                                               n_clust_2,
//                                               phi,
//                                               min_n_clust,
//                                               v,
//                                               n,
//                                               a0,
//                                               b0,
//                                               Z);
// 
//     // Separate the output into the relevant components
//     
//     // std::cout << "Values calculated now sharing out\n";
//     clust_labels_2 = arma::conv_to<arma::uvec>::from(labels_weights_phi.subvec(0, n - 1));
//     clust_weights_2 = labels_weights_phi.subvec(n, n + n_clust_2 - 1);
//     phi = arma::as_scalar(labels_weights_phi(n + n_clust_2));
//     
//     // std::cout <<"phi updated \n";
//     
//     // std::cout <<"\nContext similarity after label swapping:\n" << context_similarity << "\n\n";
//     
//     // }
//     
//     // if current iteration is a recorded iteration, save the labels
//     if (i >= burn && (i - burn) % thinning == 0) {
//       record_1.col((i - burn) / thinning) = clust_labels_1;
//       record_2.col((i - burn) / thinning) = clust_labels_2;
//       phi_record((i - burn) / thinning) = phi;
//       outlier_probs_saved_1.col((i - burn) / thinning) = outlier_vec_1;
//       outlier_probs_saved_2.col((i - burn) / thinning) = outlier_vec_2;
//     }
//     
//     if(outlier_1){
//       b_k_1 = arma::find(outlier_vec_1 == 0);
//       b_1 = b_k_1.n_elem;
//       
//       // std::cout << "Outlier weight:\n";
//       outlier_weight_1 = 1 - sample_beta(u_1 + b_1, v_1 + n - b_1);
//       // std::cout << "Outlier weight success!\n";
//     }
//     
//     if(outlier_2){
//       b_k_2 = arma::find(outlier_vec_2 == 0);
//       b_2 = b_k_2.n_elem;
//       
//       // std::cout << "Outlier weight:\n";
//       outlier_weight_2 = 1 - sample_beta(u_2 + b_2, v_2 + n - b_2);
//       // std::cout << "Outlier weight success!\n";
//     }
//   }
//   
//   // construct similarity matrix
//   arma::mat sim_1(n, n); 
//   arma::mat sim_2(n, n);
//   sim_1 = similarity_mat(record_1);
//   sim_2 = similarity_mat(record_2);
//   
//   // std::cout << "Context similarity: " << context_similarity << "\n";
//   
//   return List::create(Named("similarity_1") = sim_1,
//                       Named("similarity_2") = sim_2,
//                       Named("class_record_1") = record_1,
//                       Named("class_record_2") = record_2,
//                       Named("context_similarity") = phi_record,
//                       Named("entropy") = entropy_cw,
//                       Named("outliers_1") = outlier_probs_saved_1,
//                       Named("outliers_2") = outlier_probs_saved_2);
//   
// }

// MDI clustering for two categorical datasets
// [[Rcpp::export]]
Rcpp::List mdi_cat_cat(arma::umat data_1,
                       arma::umat data_2,
                       arma::field<arma::vec> class_dist_prior_1,
                       arma::field<arma::vec> class_dist_prior_2,
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
                       arma::uword thinning
){
  // As we tend to iterate to less than num_iter, add 1 to itself
  // Done here as many objects are sized based on this.
  num_iter++;
  
  // Declare the sample size and dimensionality of the continuous and 
  // categorical data
  
  // std::cout << "In function \n";
  
  arma::uword n = data_1.n_rows;
  arma::uword n_cols_1 = data_1.n_cols;
  
  // arma::uword n_cat = categorical_data.n_rows;
  arma::uword n_cols_2 = data_2.n_cols;
  
  // Frequently will compare clusters across contexts so this is a useful bound
  // to iterate to
  arma::uword min_n_clust = std::min(n_clust_1, n_clust_2);
  
  double v = 0.0; // strategic latent variable
  
  // Cluster weights for each dataset
  arma::vec clust_weights_1(n_clust_1);
  clust_weights_1.zeros();
  arma::vec clust_weights_2(n_clust_2);
  clust_weights_2.zeros();
  
  // Declare the field for the phi variable for the categorical data
  arma::uvec cat_count_1(n_cols_1);
  cat_count_1 = cat_counter(data_1);
  arma::field<arma::mat> class_prob_1(n_cols_1);
  
  class_prob_1 = declare_class_probs_field(cat_count_1,
                                           n_cols_1,
                                           n_clust_1);
  
  
  arma::uvec cat_count_2(n_cols_2);
  cat_count_2 = cat_counter(data_2);
  arma::field<arma::mat> class_prob_2(n_cols_2);
  
  class_prob_2 = declare_class_probs_field(cat_count_2,
                                           n_cols_2,
                                           n_clust_2);
  
  // Context similarity - smaple prior
  double phi = arma::randg(arma::distr_param(a0, 1/b0) );
  
  // Declare normalising constant
  double Z = 1.0;
  
  // Used in each iteration
  arma::vec curr_prob_vec_1(n_clust_1);
  arma::vec curr_prob_vec_2(n_clust_2);
  
  // std::cout << "Declared prob vectors \n";
  
  // the record for similarity in each clustering
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  
  arma::umat record_1(n, eff_count);
  record_1.zeros();
  
  arma::umat record_2(n, eff_count);
  record_2.zeros();
  
  arma::vec phi_record(eff_count);
  
  arma::vec labels_weights_phi(n + n_clust_2 + 1);
  
  arma::vec entropy_cw(num_iter);
  
  // std::cout << "All declared \n";
  
  
  arma::vec rate_0_1(n_clust_1);
  arma::vec rate_0_2(n_clust_2);
  
  // Placeholder prior
  rate_0_1.fill(1);
  rate_0_2.fill(1);
  
  // Not sure this is sensible - simply to have non-zero numbers here
  // clust_weights_1 = clust_weight_priors_1;
  // clust_weights_2 = clust_weight_priors_2;
  
  // Initialise v based on the prior
  arma::uword v_a_0 = 1;
  arma::uword v_b_0 = 1;
  
  v = arma::randg( arma::distr_param(v_a_0, 1.0/(v_b_0) ) );
  
  
  for(arma::uword i = 0; i < num_iter; i++){
    
    // sample cluster weights for the two datasets
    // clust_weights_1 = gamma_posterior(clust_weight_priors_1,
    //                                   clust_labels_1,
    //                                   n_clust_1);
    // 
    // clust_weights_2 = gamma_posterior(clust_weight_priors_2,
    //                                   clust_labels_2,
    //                                   n_clust_2);
    
    
    
    // std::cout << "Latent variable calculated.\n";
    
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
    
    // std::cout << "Cluster weights calculated.\n";
    
    // Entropy for graphing convergence
    entropy_cw(i) = entropy(clust_weights_1);
    
    // For the categorical data, sample the probabilities for each class
    class_prob_1 = sample_class_probabilities(data_1,
                                              class_prob_1,
                                              class_dist_prior_1,
                                              clust_labels_1,
                                              cat_count_1,
                                              n_clust_1,
                                              n_cols_1);
    
    class_prob_2 = sample_class_probabilities(data_2,
                                              class_prob_2,
                                              class_dist_prior_2,
                                              clust_labels_2,
                                              cat_count_2,
                                              n_clust_2,
                                              n_cols_2);
    
    // std::cout << "Context parameters calculated.\n";
    
    // Calculate the current normalising constant (consider being more clever 
    // about this) 
    Z = calculate_normalising_constant(clust_weights_1,
                                       clust_weights_2,
                                       phi,
                                       n_clust_1,
                                       n_clust_2);
    
    // sample the strategic latent variable, v
    v = arma::randg( arma::distr_param(v_a_0 + n, 1.0/(v_b_0 + Z) ) );
    
    // std::cout << "Z calculated \n";
    
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
    
    // std::cout << "Sampled phi\n";
    
    
    // sample the strategic latent variable, v
    // v = arma::randg( arma::distr_param(n, 1/Z) );
    
    // sample 
    for(arma::uword j = 0; j < n; j++){
      
      // for each point create the vector of probabilities associated with 
      // assignment to each cluster
      // std::cout << "Gaussian cluser prob vec next\n";
      
      
      curr_prob_vec_1 = mdi_cat_clust_prob(j, 
                                           data_1,
                                           class_prob_1,
                                           n_clust_1,
                                           n_cols_1,
                                           phi,
                                           clust_weights_1,
                                           clust_labels_1,
                                           clust_labels_2);
      
      curr_prob_vec_2 = mdi_cat_clust_prob(j, 
                                           data_2,
                                           class_prob_2,
                                           n_clust_2,
                                           n_cols_2,
                                           phi,
                                           clust_weights_2,
                                           clust_labels_2,
                                           clust_labels_1);
      
      // std::cout << "Predict label per point\n";
      
      // update labels - in gaussian data this is only if the current point is 
      // not fixed
      if(fix_vec_1[j] == 0){
        // std::cout << "Gaussian\n";
        clust_labels_1(j) = cluster_predictor(curr_prob_vec_1);
      }
      
      if(fix_vec_2[j] == 0){
        // std::cout << "Categorical\n";
        clust_labels_2(j) = cluster_predictor(curr_prob_vec_2);
      }
    }
    
    // std::cout << "All the context update stuff\n";
    
    // std::cout << cluster_labels_categorical.n_elem << "\n"; 
    
    // Update cluster labels in second dataset
    // Return the new labels, weights and similarity in a single vector
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
    
    // labels_weights_phi = cluster_label_update(cluster_labels_gaussian,
    //                                           cluster_labels_categorical,
    //                                           cluster_weights_gaussian,
    //                                           cluster_weights_categorical,
    //                                           num_clusters_categorical,
    //                                           context_similarity,
    //                                           min_num_clusters,
    //                                           v,
    //                                           n,
    //                                           a0,
    //                                           b0);
    
    // Separate the output into the relevant components
    
    // std::cout << "Values calculated now sharing out\n";
    clust_labels_2 = arma::conv_to<arma::uvec>::from(labels_weights_phi.subvec(0, n - 1));
    
    // std::cout << "cluster labels updated \n";
    
    // std::cout <<"\nCluster weights before:\n" << cluster_weights_categorical << "\n";
    
    clust_weights_2 = labels_weights_phi.subvec(n, n + n_clust_2 - 1);
    
    // std::cout <<"\nCluster weights after:\n" << cluster_weights_categorical << "\n\n";
    
    // std::cout <<"cluster weights updated \n";
    // std::cout <<"\nContext similarity before checking label swapping:\n" << context_similarity << "\n";
    
    phi = arma::as_scalar(labels_weights_phi(n + n_clust_2));
    // std::cout <<"phi updated \n";
    
    // std::cout <<"\nContext similarity after label swapping:\n" << context_similarity << "\n\n";
    // }
      
    // if current iteration is a recorded iteration, save the labels
    if (i >= burn && (i - burn) % thinning == 0) {
      record_1.col((i - burn) / thinning) = clust_labels_1;
      record_2.col((i - burn) / thinning) = clust_labels_2;
      phi_record((i - burn) / thinning) = phi;
    }
  }
  
  // construct similarity matrix
  arma::mat sim_1(n, n); 
  arma::mat sim_2(n, n);
  sim_1 = similarity_mat(record_1);
  sim_2 = similarity_mat(record_2);
  
  // std::cout << "Context similarity: " << context_similarity << "\n";
  
  return List::create(Named("similarity_1") = sim_1,
                      Named("similarity_2") = sim_2,
                      Named("class_record_1") = record_1,
                      Named("class_record_2") = record_2,
                      Named("context_similarity") = phi_record,
                      Named("entropy") = entropy_cw);
  
}


// // Checks if file_name exsits
// bool check_file_exists(char file_name[]){
//   std::ifstream file;
//   
//   file.open(file_name, std::ofstream::out);
//   return file.good();
// }
// 
// // Returns a bool based on the user's instruction to delete file_name
// bool delete_file_user(char file_name[]){
//   bool delete_choice = false;
//   std::string user_choice;
//   
//   if(check_file_exists(file_name)){
//     std::cout << "File already exists, delete? (y/n)\n"; 
//     std::cin >> user_choice;
//     while(user_choice != "y" |  user_choice != "n"){
//       std::cout << "File already exists, delete? Please select one of 'y' or 'n'"
//                 << "(y/n)\n"; 
//       std::cin >> user_choice;
//     }
//     if(user_choice == "y"){
//       delete_choice = true;
//     } else {
//       delete_choice = false;
//     }
//   }
//   return delete_choice;
// }
// 
// // Deletes input file if possible
// int delete_file(char file_name[])
// {
//   
//   /*	Deletes the file if exists */
//   if (remove(file_name) != 0)
//     perror("File deletion failed");
//   else
//     std::cout << "File deleted successfully";
//   
//   return 0;
// }
// 
// // Function to write to a file - gives the option of deleting pre-existing file
// // [[Rcpp::export]]
// int append_file (char file_name[]) {
//   
//   bool delete_choice = false;
//   std::ofstream ofs;
//   
//   delete_choice = delete_file_user(file_name);
//   
//   if(delete_choice) delete_file(file_name);
//   
//   
//   ofs.open (file_name, std::ofstream::out | std::ofstream::app);
//   
//   
//   ofs << " more lorem ipsum";
//   
//   ofs.close();
//   
//   return 0;
// }


// MDI clustering for a gaussian and cateogrical dataset
// [[Rcpp::export]]
Rcpp::List mdi_gauss_cat(arma::mat gaussian_data,
                         arma::umat categorical_data,
                         arma::vec mu_0,
                         double lambda_0,
                         arma::mat scale_0,
                         int df_0,
                         double a0,
                         double b0,
                         arma::vec cluster_weight_priors_gaussian,
                         arma::vec cluster_weight_priors_categorical,
                         arma::field<arma::vec> phi_prior,
                         arma::uvec cluster_labels_gaussian,
                         arma::uvec cluster_labels_categorical,
                         arma::uword num_clusters_gaussian,
                         arma::uword num_clusters_categorical,
                         arma::uvec fix_vec_1,
                         arma::uvec fix_vec_2,
                         arma::uword num_iter,
                         arma::uword burn,
                         arma::uword thinning,
                         bool outlier = false,
                         double t_df = 4.0,
                         bool record_posteriors = false,
                         bool normalise = false,
                         double u_1 = 2,
                         double v_1 = 10,
                         arma::uword rate_gauss_0 = 1,
                         arma::uword rate_cat_0 = 1,
                         bool save_results = false
){
  
  // As we tend to iterate to less than num_iter, add 1 to itself
  // Done here as many objects are sized based on this.
  num_iter++;
  
  // Declare the sample size and dimensionality of the continuous and 
  // categorical data
  
  // std::cout << "In function \n";
  arma::uword n = gaussian_data.n_rows;
  arma::uword num_cols_cont = gaussian_data.n_cols;
  
  // arma::uword n_cat = categorical_data.n_rows;
  arma::uword num_cols_cat = categorical_data.n_cols;
  
  
  // Frequently will compare clusters across contexts so this is a useful bound
  // to iterate to
  arma::uword min_num_clusters = std::min(num_clusters_gaussian,
                                          num_clusters_categorical);
  
  //  Normalise the continuous data
  if(normalise){
    gaussian_data = arma::normalise(gaussian_data);
  }
  
  // Declare global variance and mean - used in outlier t-distribution 
  // (if included)
  arma::mat global_variance(num_cols_cont, num_cols_cont);
  global_variance = 0.5 * arma::cov(gaussian_data); // Olly's rec
  
  arma::vec global_mean(num_cols_cont);
  global_mean = arma::trans(arma::mean(gaussian_data, 0));
  
  
  double v = 1.0; // strategic latent variable
  
  // Cluster weights for each dataset
  arma::vec cluster_weights_gaussian(num_clusters_gaussian);
  arma::vec cluster_weights_categorical(num_clusters_categorical);
  
  cluster_weights_gaussian.zeros();
  cluster_weights_categorical.zeros();
  
  // These are the local cubes of posterior mean and variance overwritten each
  // iteration
  arma::field<arma::cube> loc_mu_variance(2);
  
  // std::cout << "Declared to mean/variance thing \n";
  
  // Declare the field for the phi variable for the categorical data
  arma::uvec cat_count(num_cols_cat);
  cat_count = cat_counter(categorical_data);
  arma::field<arma::mat> class_probabilities(num_cols_cat);
  
  class_probabilities = declare_class_probs_field(cat_count,
                                                  num_cols_cat,
                                                  num_clusters_categorical);
  
  // Initialise the context similarity parameter based on priors
  double context_similarity = arma::randg(arma::distr_param(a0, 1/b0) );
  
  // Declare the normalising constant
  double Z = 1.0;
  
  arma::vec curr_gaussian_prob_vec(num_clusters_gaussian);
  arma::vec curr_categorical_prob_vec(num_clusters_categorical);
  
  // Various objects to record values for posterior distributions and clustering
  // the record for similarity in each clustering
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  
  // These are the lists recording the posterior mean and 
  // variance for each class for each recorded iteration
  ListMatrix variance(eff_count, num_clusters_gaussian);
  ListMatrix mu(eff_count, num_clusters_gaussian);
  
  // Objects required to save the categorical variable for each component
  arma::field<arma::field<arma::mat>> class_probabilities_saved(num_clusters_categorical);
  // arma::field<arma::rowvec> comp_class_probs(num_cols_cat);
  
  class_probabilities_saved = class_probs_obj(num_clusters_categorical,
                                              num_cols_cat,
                                              eff_count,
                                              cat_count);
  
  // Records of allocation
  arma::umat gaussian_record(n, eff_count);
  gaussian_record.zeros();
  
  arma::umat categorical_record(n, eff_count);
  categorical_record.zeros();
  
  arma::vec context_similarity_record(eff_count);
  
  // ## Filenames for saving posteriors ##
  std::string gauss_lab_file = "gaussian_allocation/gaussian_labels_iter_";
  std::string cat_lab_file = "categorical_allocation/categorical_labels_iter_";
  std::string out_lab_file = "outlier_allocation/outlier_labels_iter_";
  
  std::string mean_file = "mean/mean_";
  std::string var_file = "variance/var_";
  
  std::string i_str;
  std::string comp_str;
  
  std::string gauss_lab_file_loc;
  std::string cat_lab_file_loc;
  std::string out_lab_file_loc;
  
  std::string mu_file_loc;
  std::string var_file_loc;

  // A positive integer to hold the current count of recorded iterations
  arma::uword record_iter = 0;
  
  // To hold output of label flipping function
  arma::vec labels_weights_phi(n + num_clusters_categorical + 1);
  
  // A vector to record the entropy at each iteration
  arma::vec entropy_cw(num_iter);
  
  // Rate priors for weight smapling from Gamma distn
  arma::vec rate_0_gauss(num_clusters_gaussian);
  arma::vec rate_0_cat(num_clusters_categorical);
  
  // Placeholder prior
  rate_0_gauss.fill(rate_gauss_0);
  rate_0_cat.fill(rate_cat_0);

  // OUTLIER COMPONENT
  // Variables to handle outlier component from tagm
  // Vector of 0 and 1 for points assigned to outlier group or not
  arma::uvec outlier_vec(n);
  
  
  arma::uvec b_k(n); // the vector of b_k's the sum of which gives b
  double b = 0.0; // a will be N - b, no need to declare
  
  // for recording the probabilities for each class
  arma::cube gaussian_class_probs(eff_count, num_clusters_gaussian, n);
  arma::cube cat_class_probs(eff_count, num_clusters_categorical, n); 
  
  // Begin with all non-fixed points (i.e. unlabelled) to outlier component
  outlier_vec.zeros();
  if(outlier && arma::any(fix_vec_1)){
    outlier_vec = 1 - fix_vec_1;
  }
  
  // the current iterations probabilities (overwritten - saved to the above cube
  // after burn in and as defined by thinning)
  arma::vec curr_class_probs(num_clusters_gaussian);
  arma::vec curr_norm_likelihoods(num_clusters_gaussian);
  
  // Class labels of points not currently assigned as outliers
  arma::uvec relevant_labels(n);
  
  // the current iterations probabilities of being an outlier (non-outlier prob
  // is 1 - curr_outlier_prob)
  arma::vec curr_outlier_prob(2);
  double curr_outlier_likelihood = 0.0;
  arma::uword predicted_outlier = 0;
  double outlier_weight = 1 - sample_beta(u_1, v_1);
  
  // the predicted class assuming the point is not an outlier
  arma::uword predicted_class = 0;
  
  // This is where we save the outlier labels
  arma::umat outlier_probs_saved(n, eff_count);
  
  // Declare the variable to hold the likelihood of bein non-outlier
  double predicted_norm_likelihood = 0.0;
  
  // For holding the log-likelihood of each recorded iteration
  arma::vec mdi_recorded_likelihood(eff_count);
  mdi_recorded_likelihood.zeros();
  
  // The total likelihood of the gaussian model in the current iteration
  double gaussian_score = 0.0;
  
  // The total likelihood of the categorical model in the current iteration
  double cat_score = 0.0;
  
  double mdi_likelihood = 0.0;
  double log_gamma_n = log_factorial(n - 1);

  // The matrices to hold the allocation probabilities for each class
  arma::mat alloc_prob_gauss(n, num_clusters_gaussian);
  arma::mat alloc_prob_cat(n, num_clusters_categorical);
  
  // Set all entries to zero as we will be using += rather than = (i.e will 
  // never over-write the initial value held, only add to it)
  alloc_prob_gauss.zeros();
  alloc_prob_cat.zeros();
  
  for(arma::uword i = 0; i < num_iter; i++){
    
    // Consider only the labels of points not considered to be outliers
    relevant_labels = cluster_labels_gaussian % (1 - outlier_vec);
    
    // Entropy for graphing convergence
    entropy_cw(i) = entropy(cluster_weights_gaussian);
    
    // ## Sample the parameters for the two datasets ##
    
    // Sample the posterior mean and variance for the gaussian data
    loc_mu_variance = mean_variance_sampling(gaussian_data,
                                             relevant_labels,
                                             num_clusters_gaussian,
                                             df_0,
                                             num_cols_cont,
                                             scale_0,
                                             lambda_0,
                                             mu_0);
    
    // For the categorical data, sample the probabilities for each class
    class_probabilities = sample_class_probabilities(categorical_data,
                                                     class_probabilities,
                                                     phi_prior,
                                                     cluster_labels_categorical,
                                                     cat_count,
                                                     num_clusters_categorical,
                                                     num_cols_cat);

    // ## Sample cluster weights for the two datasets ##
    
    // Gaussian weights
    cluster_weights_gaussian = mdi_cluster_weights(cluster_weight_priors_gaussian,
                                                   rate_0_gauss,
                                                   v,
                                                   num_clusters_gaussian,
                                                   num_clusters_categorical,
                                                   cluster_weights_categorical,
                                                   cluster_labels_gaussian,
                                                   // relevant_labels,
                                                   cluster_labels_categorical,
                                                   context_similarity);
    
    // Categorical weights
    cluster_weights_categorical = mdi_cluster_weights(cluster_weight_priors_categorical,
                                                      rate_0_cat,
                                                      v,
                                                      num_clusters_categorical,
                                                      num_clusters_gaussian,
                                                      cluster_weights_gaussian,
                                                      cluster_labels_categorical,
                                                      cluster_labels_gaussian,
                                                      context_similarity);
    
    // Calculate the current normalising constant (consider being more clever
    // about this)
    Z = calculate_normalising_constant(cluster_weights_gaussian,
                                       cluster_weights_categorical,
                                       context_similarity,
                                       num_clusters_gaussian,
                                       num_clusters_categorical);
    
    // sample the strategic latent variable, v
    v = arma::randg( arma::distr_param(n, 1.0/Z ) );
    
    // sample the context similarity parameter (as only two contexts this is a
    // single double - number not a drink)
    context_similarity = sample_phi(cluster_labels_gaussian,
                                    cluster_labels_categorical,
                                    cluster_weights_gaussian,
                                    cluster_weights_categorical,
                                    v,
                                    n,
                                    min_num_clusters,
                                    a0,
                                    b0);
    
    // sample class for each observation
    for(arma::uword j = 0; j < n; j++){
      
      // for each point create the vector of probabilities associated with 
      // assignment to each cluster within the gaussian dataset
      curr_norm_likelihoods = mdi_gauss_clust_probs(j,
                                                    gaussian_data,
                                                    num_clusters_gaussian,
                                                    loc_mu_variance(1),
                                                    loc_mu_variance(0),
                                                    context_similarity,
                                                    cluster_weights_gaussian,
                                                    relevant_labels,
                                                    cluster_labels_categorical);
      
      // Normalise this vector
      curr_gaussian_prob_vec = over_flow_handling(curr_norm_likelihoods);
      
      // Using the rejection method sample a single class from this vector
      predicted_class = cluster_predictor(curr_gaussian_prob_vec);
      
      // The various probabilities to determine if the observation is considered 
      // an outlier or not (if using TAGM rather than traditional mixtures)
      if(outlier){
        
        // The likelihood associated with the outlier t-distribution
        curr_outlier_likelihood = sample_outlier(arma::trans(gaussian_data.row(j)),
                                                 gaussian_data,
                                                 outlier_weight,
                                                 global_mean,
                                                 global_variance,
                                                 t_df,
                                                 u_1 + b,
                                                 v_1 + n - b);
          
        // The normal likelihood for the current class allocation
        predicted_norm_likelihood = normal_likelihood(arma::trans(gaussian_data.row(j)),
                                                      loc_mu_variance(1).slice(predicted_class - 1),
                                                      loc_mu_variance(0).slice(predicted_class - 1),
                                                      num_cols_cont);
        
        predicted_norm_likelihood += log(1 - outlier_weight);
        curr_outlier_prob(0) = predicted_norm_likelihood;
        
        // Overflow handling
        curr_outlier_prob = exp(curr_outlier_prob - max(curr_outlier_prob));
        
        // Normalise the vector
        curr_outlier_prob = curr_outlier_prob / sum(curr_outlier_prob);
        
        // Check if we predict the current individual to be assigned as an 
        // outlier using the rejection method
        predicted_outlier = cluster_predictor(curr_outlier_prob) - 1; // as +1 to handle R
        
        gaussian_score += predicted_norm_likelihood * predicted_outlier 
          + curr_outlier_likelihood * (1 - predicted_outlier);
        // gaussian_score = gaussian_score / sum(cluster_weights_gaussian)
      }
      
      curr_categorical_prob_vec = mdi_cat_clust_prob(j,
                                                     categorical_data,
                                                     class_probabilities,
                                                     num_clusters_categorical,
                                                     num_cols_cat,
                                                     context_similarity,
                                                     cluster_weights_categorical,
                                                     cluster_labels_gaussian,
                                                     cluster_labels_categorical);
      
      
      
      if (i >= burn && (i - burn) % thinning == 0) { // include  && record_posteriors in if?
        record_iter = (i - burn) / thinning;
        // gaussian_class_probs.slice(j).row(record_iter) = arma::trans(curr_gaussian_prob_vec);
        // cat_class_probs.slice(j).row(record_iter) = arma::trans(curr_categorical_prob_vec);
        
        alloc_prob_gauss.row(j) += arma::trans(curr_gaussian_prob_vec);
        alloc_prob_cat.row(j) += arma::trans(curr_categorical_prob_vec);
      }
      
      // update labels - in gaussian data this is only if the current point is 
      // not fixed
      if(fix_vec_1[j] == 0){
        cluster_labels_gaussian(j) = predicted_class; // cluster_predictor(curr_gaussian_prob_vec);
        if(outlier){
          outlier_vec(j) = predicted_outlier;
        }
      }
      
      if (fix_vec_2[j] == 0){
        cluster_labels_categorical(j) = cluster_predictor(curr_categorical_prob_vec);
      }
      
      // mdi_likelihood += log(cluster_weights_gaussian(cluster_labels_gaussian(i) - 1))
      //   + log(cluster_weights_categorical(cluster_labels_categorical(i) - 1))
      //   + log(1 + context_similarity * (cluster_labels_gaussian(i) == cluster_labels_categorical(i)));
      
    }
    
    if(outlier){
      
      // Components of outlier weight
      b_k = arma::find(outlier_vec == 0);
      b = b_k.n_elem;
      
      // Sample outlier weight
      outlier_weight = 1 - sample_beta(u_1 + b, v_1 + n - b);
      
    }
    
    // Update cluster labels in second dataset
    // Return the new labels, weights and similarity in a single vector
    labels_weights_phi = cluster_label_update(cluster_labels_gaussian,
                                              cluster_labels_categorical,
                                              cluster_weights_gaussian,
                                              cluster_weights_categorical,
                                              num_clusters_gaussian,
                                              num_clusters_categorical,
                                              context_similarity,
                                              min_num_clusters,
                                              v,
                                              n,
                                              a0,
                                              b0,
                                              Z);
    
    // Separate the output into the relevant components
    cluster_labels_categorical = arma::conv_to<arma::uvec>::from(labels_weights_phi.subvec(0, n - 1));
    cluster_weights_categorical = labels_weights_phi.subvec(n, n + num_clusters_categorical - 1);
    context_similarity = arma::as_scalar(labels_weights_phi(n + num_clusters_categorical));
    
    
    
    // if current iteration is a recorded iteration, save the labels
    if (i >= burn && (i - burn) % thinning == 0) {
      
      // The log likelihood for the current iteration
      mdi_likelihood = MDI_log_likelihood(v,
                                          n,
                                          Z,
                                          context_similarity,
                                          cluster_weights_gaussian,
                                          cluster_weights_categorical,
                                          cluster_labels_gaussian,
                                          cluster_labels_categorical);
      
      record_iter = (i - burn) / thinning;
      
      // gaussian_record.col(record_iter) = cluster_labels_gaussian;
      // categorical_record.col(record_iter) = cluster_labels_categorical;
      context_similarity_record(record_iter) = context_similarity;
      // outlier_probs_saved.col(record_iter) = outlier_vec;
      
      // Record the current iteration's score (note this does not account for
      // label flipping)
      mdi_recorded_likelihood(record_iter) = mdi_likelihood;
      
      
      if(save_results){
        // Save to file
        i_str = std::to_string(record_iter);
        
        gauss_lab_file_loc = gauss_lab_file + i_str;
        cat_lab_file_loc = cat_lab_file + i_str;
        out_lab_file_loc = out_lab_file + i_str;
        
        cluster_labels_gaussian.save(gauss_lab_file_loc);
        cluster_labels_categorical.save(cat_lab_file_loc);
        outlier_vec.save(out_lab_file_loc);
      } else {
        gaussian_record.col(record_iter) = cluster_labels_gaussian;
        categorical_record.col(record_iter) = cluster_labels_categorical;
        outlier_probs_saved.col(record_iter) = outlier_vec;
      }
      // Record posteriors of parameters for Gaussian and Categorical
      // distributions
      if(record_posteriors){
        for(arma::uword j = 0; j < num_clusters_gaussian; j++){
          if(save_results){
            // Save to file
            
            // The current component
            comp_str = std::to_string(j);
            mu_file_loc = mean_file + comp_str + "/iter_" + i_str;
            var_file_loc = var_file + comp_str + "/iter_" + i_str;
            
            loc_mu_variance(1).slice(j).save(mu_file_loc);
            loc_mu_variance(0).slice(j).save(var_file_loc);
          } else {
            mu(record_iter, j) = loc_mu_variance(1).slice(j);
            variance(record_iter, j) = loc_mu_variance(0).slice(j);
          }
        }
        
        for(arma::uword j = 0; j < num_clusters_categorical; j++){
          for(arma::uword k = 0; k < num_cols_cat; k++){
            // if(save_results){
            //   
            //   comp_str = std::to_string(j);
            //   clust_str = std::to_string(k);
            //   
            //   class_probs_param_file = phi_file + comp_str + "m";
            //   
            //   
            // } else{
              class_probabilities_saved(j)(k).row(record_iter) = class_probabilities(k).row(j);
            // }
          }
        }
      }
    }
    
  }
  
  // Loading posterior objects
  for(arma::uword i = 0; i < eff_count; i++){
    
    if(save_results){
      i_str = std::to_string(record_iter);
      
      gauss_lab_file_loc = gauss_lab_file + i_str;
      cat_lab_file_loc = cat_lab_file + i_str;
      out_lab_file_loc = out_lab_file + i_str;
      
      cluster_labels_gaussian.load(gauss_lab_file_loc);
      cluster_labels_categorical.load(cat_lab_file_loc);
      outlier_vec.load(out_lab_file_loc);
      
      gaussian_record.col(i) = cluster_labels_gaussian;
      categorical_record.col(i) = cluster_labels_categorical;
      outlier_probs_saved.col(i) = outlier_vec;
    }
    
    // if(record_posteriors){
    //   for(arma::uword j = 0; j < num_clusters_gaussian; j++){
    //     mu(record_iter, j) = loc_mu_variance(1).slice(j);
    //     variance(record_iter, j) = loc_mu_variance(0).slice(j);
    //   }
    //   
    //   for(arma::uword j = 0; j < num_clusters_categorical; j++){
    //     for(arma::uword k = 0; k < num_cols_cat; k++){
    //       class_probabilities_saved(j)(k).row(record_iter) = class_probabilities(k).row(j);
    //     }
    //   }
    // }
    
  }
  
  alloc_prob_gauss = alloc_prob_gauss / eff_count;
  alloc_prob_cat = alloc_prob_cat / eff_count;
  
  // construct similarity matrices
  arma::mat sim(n, n); 
  arma::mat cat_sim(n, n);
  sim = similarity_mat(gaussian_record);
  cat_sim = similarity_mat(categorical_record);
  
  return List::create(Named("similarity_1") = sim,
                      Named("similarity_2") = cat_sim,
                      Named("allocation_mat_1") = alloc_prob_gauss,
                      Named("allocation_mat_2") = alloc_prob_cat,
                      Named("class_record_1") = gaussian_record,
                      Named("class_record_2") = categorical_record,
                      Named("mean_posterior") = mu,
                      Named("variance_posterior") = variance,
                      Named("class_prob_posterior") = class_probabilities_saved,
                      Named("context_similarity") = context_similarity_record,
                      Named("entropy") = entropy_cw,
                      Named("outlier") = outlier_probs_saved,
                      Named("likelihood") = mdi_recorded_likelihood);
}


///////////////////////////////////////////////////////////////////////////////

// MDI for different types

///////////////////////////////////////////////////////////////////////////////

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
                           std::string lab_file_1 = "allocation_1/labels_1_iter_",
                           std::string lab_file_2 = "allocation_2/labels_2_iter_",
                           std::string out_lab_file_1 = "outlier_allocation_1/outlier_1_labels_iter_",
                           std::string out_lab_file_2 = "outlier_allocation_2/outlier_2_labels_iter_"
){
  
  // As we tend to iterate to less than num_iter, add 1 to itself
  // Done here as many objects are sized based on this.
  num_iter++;
  
  // Declare the sample size and dimensionality of the continuous and 
  // categorical data
  
  // std::cout << "In function \n";
  
  arma::uword n = data_1.n_rows;
  arma::uword n_cols_1 = data_1.n_cols;
  
  // arma::uword n_cat = categorical_data.n_rows;
  arma::uword n_cols_2 = data_2.n_cols;
  
  // Frequently will compare clusters across contexts so this is a useful bound
  // to iterate to
  arma::uword min_n_clust = std::min(n_clust_1, n_clust_2);
  
  // std::cout << "Nomalising data\n";
  
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
  double v = 0.0;
  
  // Cluster weights for each dataset
  arma::vec clust_weights_1(n_clust_1);
  clust_weights_1.zeros();
  
  arma::vec clust_weights_2(n_clust_2);
  clust_weights_2.zeros();
  
  // These are the local cubes of posterior mean and variance overwritten each
  // iteration
  arma::field<arma::cube> loc_mu_variance_1(2);
  arma::field<arma::cube> loc_mu_variance_2(2);
  
  // std::cout << "Generic message. Parameter fields declared.\n";
  
  // Context similarity - sample prior
  double phi = arma::randg(arma::distr_param(a0, 1/b0) );
  
  // Declare normalising constant
  double Z = 1.0;
  
  // Used in each iteration
  arma::vec curr_prob_vec_1(n_clust_1);
  arma::vec curr_prob_vec_2(n_clust_2);
  
  // Various objects to record values for posterior distributions and clustering
  // the record for similarity in each clustering
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  
  // These are the lists recording the posterior mean and 
  // variance for each class for each recorded iteration
  ListMatrix variance_1(eff_count, n_clust_1);
  ListMatrix mu_1(eff_count, n_clust_1);
  
  ListMatrix variance_2(eff_count, n_clust_2);
  ListMatrix mu_2(eff_count, n_clust_2);
  
  // Record the class allocations
  arma::umat record_1(n, eff_count);
  record_1.zeros();
  
  arma::umat record_2(n, eff_count);
  record_2.zeros();
  
  // The vector to record the context similarity parameter
  arma::vec phi_record(eff_count);
  
  arma::vec labels_weights_phi(n + n_clust_2 + 1);
  
  arma::vec entropy_cw(num_iter);
  
  // Filenames for saving posteriors
  // std::string gauss_lab_file = "gaussian_allocation/gaussian_labels_iter_";
  // std::string cat_lab_file = "categorical_allocation/categorical_labels_iter_";
  // std::string out_lab_file = "outlier_allocation/outlier_labels_iter_";
  std::string i_str;
  
  std::string lab_file_1_loc;
  std::string lab_file_2_loc;
  std::string out_lab_file_1_loc;
  std::string out_lab_file_2_loc;
  
  arma::uword record_iter = 0;
  
  
  // std::cout << "Entropy vec declared.\n";
  
  arma::vec rate_0_1(n_clust_1);
  arma::vec rate_0_2(n_clust_2);
  
  // Placeholder prior
  rate_0_1.fill(rate_1_0);
  rate_0_2.fill(rate_2_0);
  
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
  arma::uvec b_k_1(n);
  double b_1 = 0.0; 
  
  arma::uvec b_k_2(n);
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
  
  // Initialise v based on the prior
  arma::uword v_a_0 = 1;
  arma::uword v_b_0 = 1;
  
  // v = arma::randg( arma::distr_param(v_a_0, 1.0/(v_b_0) ) );
  v = 1.0;
  
  
  for(arma::uword i = 0; i < num_iter; i++){
    
    // sample cluster weights for the two datasets
    // clust_weights_1 = gamma_posterior(clust_weight_priors_1,
    //                                   clust_labels_1,
    //                                   n_clust_1);
    
    
    // Consider only the labels of points not considered to be outliers
    relevant_labels_1 = clust_labels_1 % (1 - outlier_vec_1);
    relevant_labels_2 = clust_labels_2 % (1 - outlier_vec_2);
    
    // std::cout << "\nStrategic latent variable sampling\n";
    
    
    // // sample the strategic latent variable, v
    // v = arma::randg( arma::distr_param(n, 1/Z) );
    // 
    
    // Entropy for graphing convergence
    entropy_cw(i) = entropy(clust_weights_1);
    
    // Sample the posterior mean and variance for the first dataset
    loc_mu_variance_1 = mean_variance_sampling(data_1,
                                               relevant_labels_1,
                                               n_clust_1,
                                               df_0_1,
                                               n_cols_1,
                                               scale_0_1,
                                               lambda_0_1,
                                               mu_0_1);
    
    // std::cout << "Variance sampled\n";
    
    // Sample the posterior mean and variance for the second dataset
    loc_mu_variance_2 = mean_variance_sampling(data_2,
                                               relevant_labels_2,
                                               n_clust_2,
                                               df_0_2,
                                               n_cols_2,
                                               scale_0_2,
                                               lambda_0_2,
                                               mu_0_2);
    
    // std::cout << "Sampled parameters for both datasets\n";
    
    // Sample cluster weights within each context
    clust_weights_1 = mdi_cluster_weights(clust_weight_priors_1,
                                          rate_0_1,
                                          v,
                                          n_clust_1,
                                          n_clust_2,
                                          clust_weights_2,
                                          relevant_labels_1,
                                          clust_labels_2,
                                          phi);
    
    clust_weights_2 = mdi_cluster_weights(clust_weight_priors_2,
                                          rate_0_2,
                                          v,
                                          n_clust_2,
                                          n_clust_1,
                                          clust_weights_1,
                                          relevant_labels_2,
                                          clust_labels_1,
                                          phi);
    
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
    
    // std::cout << "Sampled phi\n";
    
    // Calculate the current normalising constant (consider being more clever 
    // about this) 
    Z = calculate_normalising_constant(clust_weights_1,
                                       clust_weights_2,
                                       phi,
                                       n_clust_1,
                                       n_clust_2);
    
    // std::cout << "Z calculated \n";
    
    // sample the strategic latent variable, v
    v = arma::randg( arma::distr_param(v_a_0 + n, 1.0/(v_b_0 + Z) ) );
    
    // sample 
    for(arma::uword j = 0; j < n; j++){
      
      // Calculate the log-likelihoods for each cluster in each context
      norm_likelihoods_1 = mdi_gauss_clust_probs(j,
                                                 data_1,
                                                 n_clust_1,
                                                 loc_mu_variance_1(1),
                                                 loc_mu_variance_1(0),
                                                 phi,
                                                 clust_weights_1,
                                                 relevant_labels_1,
                                                 relevant_labels_2);
      
      norm_likelihoods_2 = mdi_gauss_clust_probs(j,
                                                 data_2,
                                                 n_clust_2,
                                                 loc_mu_variance_2(1),
                                                 loc_mu_variance_2(0),
                                                 phi,
                                                 clust_weights_2,
                                                 relevant_labels_2,
                                                 relevant_labels_1);
      
      // std::cout << "Normal likelihoods calculated.\n";
      
      // Convert to likelihoods, handle overflow and normalise
      curr_prob_vec_1 = over_flow_handling(norm_likelihoods_1);
      curr_prob_vec_2 = over_flow_handling(norm_likelihoods_2);
      
      // Predict the component to which the obersation belongs
      predicted_class_1 = cluster_predictor(curr_prob_vec_1);
      predicted_class_2 = cluster_predictor(curr_prob_vec_2);
      
      // std::cout << "Classes predicted.\n";
      
      if(outlier_1){
        
        // std::cout << "Outlier 1.\n";
        
        // Calculate the likelihood associated with the outlier component
        outlier_likelihood_1 = sample_outlier(arma::trans(data_1.row(j)),
                                              data_1,
                                              outlier_weight_1,
                                              global_mean_1,
                                              global_variance_1,
                                              t_df_1,
                                              u_1 + b_1,
                                              v_1 + n - b_1);
        
        // std::cout << "Outlier likelihood calculated.\n";
        
        // Put it into a vector (for normalising and whatnot)
        curr_outlier_prob_1(1) = outlier_likelihood_1;
        
        predicted_norm_likelihood_1 = normal_likelihood(arma::trans(data_1.row(j)),
                                                        loc_mu_variance_1(1).slice(predicted_class_1 - 1),
                                                        loc_mu_variance_1(0).slice(predicted_class_1 - 1),
                                                        n_cols_1);
        
        predicted_norm_likelihood_1 += log(1 - outlier_weight_1);
        curr_outlier_prob_1(0) = predicted_norm_likelihood_1;
        
        
        // Overflow handling
        curr_outlier_prob_1 = exp(curr_outlier_prob_1 - max(curr_outlier_prob_1));
        
        // Normalise the vector
        curr_outlier_prob_1 = curr_outlier_prob_1 / sum(curr_outlier_prob_1);
        
        // Predict if the current observation is an outlier or not  
        predicted_outlier_1 = cluster_predictor(curr_outlier_prob_1) - 1; // as +1 to handle R
        
      }
      
      // Similarly for context 2
      if(outlier_2){
        
        // Calculate the likelihood associated with the outlier component
        outlier_likelihood_2 = sample_outlier(arma::trans(data_2.row(j)),
                                              data_2,
                                              outlier_weight_2,
                                              global_mean_2,
                                              global_variance_2,
                                              t_df_2,
                                              u_2 + b_2,
                                              v_2 + n - b_2);
        
        // std::cout << "Outlier likelihood calculated.\n";
        
        // Put it into a vector (for normalising and whatnot)
        curr_outlier_prob_2(1) = outlier_likelihood_2;
        
        predicted_norm_likelihood_2 = normal_likelihood(arma::trans(data_2.row(j)),
                                                        loc_mu_variance_2(1).slice(predicted_class_2 - 1),
                                                        loc_mu_variance_2(0).slice(predicted_class_2 - 1),
                                                        n_cols_2);
        
        predicted_norm_likelihood_2 += log(1 - outlier_weight_2);
        curr_outlier_prob_2(0) = predicted_norm_likelihood_2;
        
        
        // Overflow handling
        curr_outlier_prob_2 = exp(curr_outlier_prob_2 - max(curr_outlier_prob_2));
        
        // Normalise the vector
        curr_outlier_prob_2 = curr_outlier_prob_2 / sum(curr_outlier_prob_2);
        
        // Predict if the current observation is an outlier or not  
        predicted_outlier_2 = cluster_predictor(curr_outlier_prob_2) - 1; // as +1 to handle R
      }
      
      // std::cout << "\nContext 1 prediction: " << predicted_class_1
      //           << "\nContext 2 prediction: " << predicted_class_2 << "\n";
      
      // update labels if current point's label is not fixed
      if(fix_vec_1[j] == 0){
        // std::cout << "Gaussian 1\n";
        clust_labels_1(j) = predicted_class_1;
        if(outlier_1){
          // std::cout << "Gaussian 1 outlier\n";
          outlier_vec_1(j) = predicted_outlier_1;
        }
      }
      
      if(fix_vec_2[j] == 0){
        // std::cout << "Gaussian 2\n";
        clust_labels_2(j) = predicted_class_2;
        if(outlier_2){
          // std::cout << "Gaussian 2 outlier\n";
          outlier_vec_2(j) = predicted_outlier_2;
        }
      }
                         
    }
    
    // Update cluster labels in second dataset
    // Return the new labels, weights and similarity in a single vector
    
    // Do not allow label flipping if any of context 2 have fixed labels
    
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
    
    // Separate the output into the relevant components
    
    // std::cout << "Values calculated now sharing out\n";
    clust_labels_2 = arma::conv_to<arma::uvec>::from(labels_weights_phi.subvec(0, n - 1));
    clust_weights_2 = labels_weights_phi.subvec(n, n + n_clust_2 - 1);
    phi = arma::as_scalar(labels_weights_phi(n + n_clust_2));
    
    // }
    
    // if current iteration is a recorded iteration, save the labels
    if (i >= burn && (i - burn) % thinning == 0) {
      
      phi_record((i - burn) / thinning) = phi;
      
      // Variable to create correct iteration for filenames
      record_iter = (i - burn) / thinning;
      
      if(save_results){
      // Save to file
      i_str = std::to_string(record_iter);
      
      // Create current interation file name
      lab_file_1_loc = lab_file_1 + i_str;
      lab_file_2_loc = lab_file_2 + i_str;
      out_lab_file_1_loc = out_lab_file_1 + i_str;
      out_lab_file_2_loc = out_lab_file_2 + i_str;
      
      // Save
      clust_labels_1.save(lab_file_1_loc);
      clust_labels_2.save(lab_file_2_loc);
      
      outlier_vec_1.save(out_lab_file_1_loc);
      outlier_vec_2.save(out_lab_file_2_loc);
      } else {
        
        record_1.col((i - burn) / thinning) = clust_labels_1;
        record_2.col((i - burn) / thinning) = clust_labels_2;
        outlier_probs_saved_1.col((i - burn) / thinning) = outlier_vec_1;
        outlier_probs_saved_2.col((i - burn) / thinning) = outlier_vec_2;
        
      }
      
      // Record posteriors of parameters for Gaussian and Categorical
      // distributions
      if(record_posteriors){
        for(arma::uword j = 0; j < n_clust_1; j++){
          mu_1(record_iter, j) = loc_mu_variance_1(1).slice(j);
          variance_1(record_iter, j) = loc_mu_variance_1(0).slice(j);
        }
        for(arma::uword j = 0; j < n_clust_2; j++){
          mu_2(record_iter, j) = loc_mu_variance_2(1).slice(j);
          variance_2(record_iter, j) = loc_mu_variance_2(0).slice(j);
        }
      }
    }
    
    if(outlier_1){
      b_k_1 = arma::find(outlier_vec_1 == 0);
      b_1 = b_k_1.n_elem;
      
      // std::cout << "Outlier weight:\n";
      outlier_weight_1 = 1 - sample_beta(u_1 + b_1, v_1 + n - b_1);
      // std::cout << "Outlier weight success!\n";
    }
    
    if(outlier_2){
      b_k_2 = arma::find(outlier_vec_2 == 0);
      b_2 = b_k_2.n_elem;
      
      // std::cout << "Outlier weight:\n";
      outlier_weight_2 = 1 - sample_beta(u_2 + b_2, v_2 + n - b_2);
      // std::cout << "Outlier weight success!\n";
    }
  }
  
  // construct similarity matrix
  arma::mat sim_1(n, n); 
  arma::mat sim_2(n, n);
  sim_1 = similarity_mat(record_1);
  sim_2 = similarity_mat(record_2);
  
  // std::cout << "Context similarity: " << context_similarity << "\n";
  
  return List::create(Named("similarity_1") = sim_1,
                      Named("similarity_2") = sim_2,
                      Named("class_record_1") = record_1,
                      Named("class_record_2") = record_2,
                      Named("context_similarity") = phi_record,
                      Named("entropy") = entropy_cw,
                      Named("outliers_1") = outlier_probs_saved_1,
                      Named("outliers_2") = outlier_probs_saved_2);
                     
}
