# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

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
  global_variance = 0.5 * arma::cov(data); 
  
  arma::vec global_mean(num_cols);
  global_mean = arma::trans(arma::mean(data, 0));
  
  arma::vec entropy_cw(num_iter);
  
  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  arma::uword record_ind;
  
  arma::umat record(N, eff_count);
  
  // for recording the probabilities for being declared an outlier
  arma::umat outlier_probs_saved(N, eff_count);
  record.zeros();
  
  arma::mat sim(N, N);
  sim.zeros();
  
  arma::mat cluster_data;
  cluster_data.zeros();
  
  
  // These are the lists recording the posterior mean and 
  // variance for each class for each recorded iteration
  arma::vec point;
  
  arma::vec class_weights(k);
  class_weights.zeros();
  
  arma::vec outlier_weights(2);
  outlier_weights.zeros();
  
  // for the calculation of the a, b parameters in the posterior of the outlier 
  // class weight (epsilon in Crook et al 2018)
  double b = 0.0; // a will be N - b, no need to declare
  
  // Vector of 0 and 1 for points assigned to outlier group or not
  arma::uvec outlier_vec(N);
  outlier_vec.zeros();
  
  // Begin with all non-fixed points (i.e. unlabelled) to outlier component
  outlier_vec = 1 - fix_vec;
  
  // the current iterations probabilities (overwritten - saved to the above cube
  // after burn in and as defined by thinning)
  arma::vec cluster_prob(k);
  
  // Class labels of points not currently assigned as outliers
  arma::uvec rel_labels(N);
  rel_labels.zeros();
  
  // the current iterations probabilities of being an outlier (non-outlier prob
  // is 1 - outlier_prob)
  arma::vec outlier_prob(2);
  outlier_prob.zeros();
  arma::uword predicted_outlier = 0;
  double outlier_weight = 0.0;

  // the predicted class assuming the point is not an outlier
  arma::uword predicted_class = 0;
  
  double non_outlier_weight = 0.0;
  arma::mat mu_n(num_cols, k);
  arma::cube variance_n(num_cols, num_cols, k);
  
  for(arma::uword i = 0; i < num_iter; i++){
    
    // To see which points are relevant for defining component parameters
    // use pairwise multiplication between the current label and the outlier
    rel_labels = class_labels % (1 - outlier_vec);
    
    // Update class weights
    class_weights = dirichlet_posterior(concentration_0, class_labels, k);
    
    // Calculate the entropy of the current weights
    entropy_cw(i) = entropy(class_weights);
    
    // Sample cluster specific parameters
    variance_n = sample_cluster_variance(data,
                                         rel_labels,
                                         k,
                                         df_0,
                                         num_cols,
                                         scale_0,
                                         lambda_0,
                                         mu_0);
    
    mu_n = sample_cluster_means(data,
                                rel_labels,
                                k,
                                num_cols,
                                variance_n,
                                lambda_0,
                                mu_0);
    
    // If outliers are allowed calculate the outlier component weight
    if(outlier){
      b = (double) sum(outlier_vec);
      outlier_weight = sample_beta(b + u, N + v - b);
      non_outlier_weight = 1.0 - outlier_weight;
    }
    
    
    for (arma::uword jj = 0; jj < N; jj++){
      // sample class allocation
      point = arma::trans(data.row(jj));
      
      cluster_prob = sample_gaussian_cluster(point, 
                                             data,
                                             k, 
                                             class_weights, 
                                             mu_n,
                                             variance_n
      );
      
      // Predict the label (+1 due to R using a 1:n system rather than 0:(n-1))
      predicted_class = predict_ind(cluster_prob) + 1; 
      
      if(outlier){
        
        // Calculate the probabilities of being an outlier or not
        outlier_prob = calculate_outlier_prob(point, 
                                              global_mean,
                                              global_variance,
                                              num_cols,
                                              t_df,
                                              outlier_weight,
                                              mu_n.col(predicted_class - 1),
                                              variance_n.slice(predicted_class - 1));
        
        // Predict membership in the outlier set or not
        predicted_outlier = predict_ind(outlier_prob);
      }
      
      // If point is not known, update allocation
      if(fix_vec[jj] == 0){
        class_labels(jj) = predicted_class;
        outlier_vec(jj) = predicted_outlier;
      }
    }
    
    
    // Record output
    if (i >= burn && (i - burn) % thinning == 0) {
      record_ind = (i - burn) / thinning;
      record.col(record_ind) = class_labels;
      outlier_probs_saved.col(record_ind) = outlier_vec;
      
    }
  }
  
  // Calculate similarity matrix
  sim = similarity_mat(record);
  
  // Convert sum of recorded allocaiton probabilities to average probability
  alloc_prob_gauss = alloc_prob_gauss / eff_count; 
  
  return List::create(Named("similarity") = sim,
                      Named("class_record") = record,
                      Named("entropy") = entropy_cw,
                      Named("outliers") = outlier_probs_saved);
}
