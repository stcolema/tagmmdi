# include <RcppArmadillo.h>
// # include <iostream>
// # include <fstream>
# include "common_functions.h"

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
) {
  
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
  arma::uword phi_count = 0; // index for recording phis
  double v = 1.0; // strategic latent variable
  double b = 0.0; // a will be N - b, no need to declare
  double context_similarity = 0.0;
  double Z = 1.0; // Declare the normalising constant
  double curr_outlier_likelihood = 0.0;
  double outlier_weight = 1.0 - sample_beta(u_1, v_1);
  double predicted_norm_likelihood = 0.0; // likelihood of being non-outlier
  double gaussian_score = 0.0; // Likelihood of the current gaussian model
  double cat_score = 0.0; // Likelihood of the current categorical model
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
  arma::mat sim(n, n); 
  arma::mat cat_sim(n, n);
  arma::cube variance_n(n_cols_cont, n_cols_cont, n_clust_gauss);
  arma::field<arma::mat> class_probabilities(n_cols_cat);
  
  //  Normalise the continuous data
  if(normalise){
    cont_data = arma::normalise(cont_data);
  }
  
  // Declare global variance and mean - used in outlier t-distribution 
  // (if included)
  global_variance = 0.5 * arma::cov(cont_data); // Olly's rec
  global_mean = arma::trans(arma::mean(cont_data, 0) );
  
  // Cluster weights for each dataset
  clust_wgt_gauss.zeros();
  clust_wgt_cat.zeros();
  
  // Posterior mean and variance
  mu_n.zeros();
  variance_n.zeros();

  // Declare the field for the phi variable for the categorical data
  cat_count = cat_counter(cat_data);
  class_probabilities = DeclareClassProbsField(cat_count,
                                               n_cols_cat,
                                               n_clust_cat
  );
  
  // Initialise the context similarity parameter based on priors
  context_similarity = arma::randg(arma::distr_param(a_0, 1/b_0) );

  // Records of allocation
  gaussian_record.zeros();
  categorical_record.zeros();
  
  rate_0_gauss.fill(rate_gauss_0);
  rate_0_cat.fill(rate_cat_0);

  // OUTLIER COMPONENT
  // Begin with all non-fixed points (i.e. unlabelled) to outlier component
  outlier_vec.zeros();
  if(allow_outliers && arma::any(fix_vec_1) ) {
    outlier_vec = 1 - fix_vec_1;
  }

  // Set all entries to zero as we will be using += rather than = (i.e will 
  // never over-write the initial value held, only add to it)
  alloc_prob_gauss.zeros();
  alloc_prob_cat.zeros();
  
  for(arma::uword i = 0; i < n_iter; i++) {
    
    // Consider only the labels of points not considered to be outliers
    relevant_labels = clust_labels_gauss % (1 - outlier_vec);
    
    // Entropy for graphing convergence
    entropy_cw(i) = CalcEntropy(clust_wgt_gauss);
    
    // ## Sample the parameters for the two datasets ##
    variance_n = SampleClusterVariance(cont_data,
                                       relevant_labels,
                                       n_clust_gauss,
                                       nu_0,
                                       n_cols_cont,
                                       scale_0,
                                       lambda_0,
                                       mu_0
    );
    
    mu_n = SampleClusterMeans(cont_data,
                              relevant_labels,
                              n_clust_gauss,
                              n_cols_cont,
                              variance_n,
                              lambda_0,
                              mu_0
    );
    
    // For the categorical data, sample the probabilities for each class
    class_probabilities = SampleCategoryProbabilities(cat_data,
                                                      class_probabilities,
                                                      phi_prior,
                                                      clust_labels_cat,
                                                      cat_count,
                                                      n_clust_cat,
                                                      n_cols_cat
    );
    
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
                                              context_similarity
    );
    
    // Categorical weights
    clust_wgt_cat = SampleMDIClusterWeights(clust_wgt_priors_cat,
                                            rate_0_cat,
                                            v,
                                            n_clust_cat,
                                            n_clust_gauss,
                                            clust_wgt_gauss,
                                            clust_labels_cat,
                                            clust_labels_gauss,
                                            context_similarity
    );
    
    // Calculate the current normalising constant (consider being more clever
    // about this)
    Z = CalcNormalisingConst(clust_wgt_gauss,
                             clust_wgt_cat,
                             context_similarity,
                             n_clust_gauss,
                             n_clust_cat);
    
    // sample the strategic latent variable, v
    v = arma::randg( arma::distr_param(n, 1.0/Z) );
    
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
                                      b_0
    );

    // If using a t-augmented Gaussian mixture rather than a Gaussian mixture
    if(allow_outliers) {

      // Components of outlier weight
      b = sum(outlier_vec);
      outlier_weight = SampleBetaDistn(b + u_1, n + v_1 - b);
      
    }
    
    // sample class for each observation
    for(arma::uword j = 0; j < n; j++) {
      
      // for each point create the vector of probabilities associated with 
      // assignment to each cluster within the gaussian dataset
      curr_norm_likelihoods = SampleMDIGaussClustProbs(j,
                                                       cont_data,
                                                       n_clust_gauss,
                                                       mu_n,
                                                       variance_n,
                                                       context_similarity,
                                                       clust_wgt_gauss,
                                                       relevant_labels,
                                                       clust_labels_cat
      );
      
      // Normalise this vector
      // Predict membership
      predicted_class =  PredictIndex(curr_norm_likelihoods) + 1;
     
      // The various probabilities to determine if the observation is considered 
      // an outlier or not (if using TAGM rather than traditional mixtures)
      if(allow_outliers) {
        
        // The likelihood associated with the outlier t-distribution
        curr_outlier_likelihood = CalcTdistnLikelihood( (cont_data.row(j) ).t(), 
                                                       global_mean, 
                                                       global_variance, 
                                                       n_cols_cont,
                                                       t_df
        );
        
        curr_outlier_likelihood += log(outlier_weight);
        
        curr_outlier_prob(1) = curr_outlier_likelihood;
        
        // The normal likelihood for the current class allocation
        predicted_norm_likelihood = normal_likelihood( (cont_data.row(j) ).t(),
                                                      mu_n.col(predicted_class - 1),
                                                      variance_n.slice(predicted_class - 1),
                                                      n_cols_cont);
        
        predicted_norm_likelihood += log(1 - outlier_weight);
        curr_outlier_prob(0) = predicted_norm_likelihood;
        
        // Predict if point is considered an outlier in current cluster
        predicted_outlier =  PredictIndex(curr_outlier_prob);
    
      }
      
      curr_cat_prob_vec = SampleMDICatClustProb(j,
                                                cat_data,
                                                class_probabilities,
                                                n_clust_cat,
                                                n_cols_cat,
                                                context_similarity,
                                                clust_wgt_cat,
                                                clust_labels_gauss,
                                                clust_labels_cat
      );
      
      
      
      if (i >= burn && (i - burn) % thin == 0) {
        
        alloc_prob_gauss.row(j) += arma::trans(curr_gauss_prob_vec);
        alloc_prob_cat.row(j) += arma::trans(curr_cat_prob_vec);
        
        alloc_prob_gauss_curr.row(j) = arma::trans(curr_gauss_prob_vec);
        alloc_prob_cat_curr.row(j) = arma::trans(curr_cat_prob_vec);
        
      }
      
      // update labels - in gaussian data this is only if the current point is 
      // not fixed
      if(fix_vec_1[j] == 0) {
        clust_labels_gauss(j) = predicted_class; // cluster_predictor(curr_gauss_prob_vec);
        if(allow_outliers) {
          outlier_vec(j) = predicted_outlier;
        }
      }
      
      if (fix_vec_2[j] == 0) {
        clust_labels_cat(j) = PredictIndex(curr_cat_prob_vec) + 1;
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
    clust_labels_cat = arma::conv_to<arma::uvec>::from(labels_weights_phi.subvec(0, n - 1) );
    clust_wgt_cat = labels_weights_phi.subvec(n, n + n_clust_cat - 1);
    context_similarity = arma::as_scalar(labels_weights_phi(n + n_clust_cat) );
    
    // if current iteration is a recorded iteration, save the labels
    if (i >= burn && (i - burn) % thin == 0) {
      
      // Record the current phi
      context_similarity_record(phi_count) = context_similarity;
      phi_count++;
      
      gaussian_record.col(record_iter) = clust_labels_gauss;
      categorical_record.col(record_iter) = clust_labels_cat;
      outlier_probs_saved.col(record_iter) = outlier_vec;

      record_iter++;
    }
  }
  
  // Normalise the allocation probabilities
  alloc_prob_gauss = alloc_prob_gauss / eff_count;
  alloc_prob_cat = alloc_prob_cat / eff_count;
  
  // construct similarity matrices
  sim = CreateSimilarityMat(gaussian_record);
  cat_sim = CreateSimilarityMat(categorical_record);
  
  return List::create(Named("similarity_1") = sim,
                      Named("similarity_2") = cat_sim,
                      Named("allocation_mat_1") = alloc_prob_gauss,
                      Named("allocation_mat_2") = alloc_prob_cat,
                      Named("class_record_1") = gaussian_record,
                      Named("class_record_2") = categorical_record,
                      Named("context_similarity") = context_similarity_record,
                      Named("entropy") = entropy_cw,
                      Named("outlier") = outlier_probs_saved
  );
}
