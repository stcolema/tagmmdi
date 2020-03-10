# include <RcppArmadillo.h>
# include "CommonFunctions.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// --- Hyperparameters ---------------------------------------------------------

//' Returns a variable involved in updating the scale, it is similar to sample 
//' covariance
//' 
//' @param data Data in matrix format
//' @param sample_mean Sample mean for data
//' @param n The number of samples in data
//' @param n_col The number of columns in data
//' 
//' @return One of the parameters required to calculate the posterior of the
//'  Multivariate normal with uknown mean and covariance (the unnormalised 
//'  sample covariance).
arma::mat calcSampleCov(arma::mat data,
                        arma::vec sample_mean,
                        arma::uword n,
                        arma::uword n_col
) {
  
  arma::mat sample_covariance = arma::zeros<arma::mat>(n_col, n_col);
  // sample_covariance.zeros();
  
  // If n > 0 (as this would crash for empty clusters)
  if(n > 0){
    data.each_row() -= arma::trans(sample_mean);
    sample_covariance = arma::trans(data) * data;
  }
  return sample_covariance;
}

//' Returns a vector of the mean of a cluster of size n
//' 
//' @param lambda_0 The prior on the shrinkage of the variance estimated from 
//' the Inverse Wishart as the variance of the mean parameter.
//' @param mu_0 The prior on the mean.
//' @param n The number of samples in the pertinent dataset.
//' @param n_col The number of columns in the pertinent dataset.
//' @param sample_mean The sample mean of the pertinent dataset.
//' 
//' @return The updated values of the mean parameter for the Multivariate Normal.
arma::vec updateMu(double lambda_0,
                   arma::vec mu_0,
                   arma::uword n,
                   arma::uword n_col,
                   arma::vec sample_mean) {
    
  arma::vec mu_n(n_col);
  
  mu_n = (lambda_0 * mu_0 + n * sample_mean) / (double)(lambda_0 + n);
  return mu_n;
}

// Returns the matrix for the scale hyperparameter for n observations
//' @param scale_0 The prior on the scale parameter of the Inverse Wishart 
//' distribution.
//' @param sample_cov The sample covariance of the pertinent dataset.
//' @param lambda The prior on the shrinkage of the variance.
//' @param n The number of pobservations in the pertinent dataset.
//' @param sample_mean The sample mean from the pertinent dataset.
//' @param mu_0 The prior on the mean parameter.
//' @param d The number of columns in the dataset used.
//' 
//' @return The updated scale parameter.
arma::mat updateScale(arma::mat scale_0, 
                     arma::mat sample_cov, 
                     double lambda, 
                     int n, 
                     arma::vec sample_mean, 
                     arma::vec mu_0,
                     int d){
  
  arma::mat dist_from_prior(d, d);
  arma::mat scale_n(d, d);
  
  // Calculate the distance of the sample mean from the prior
  dist_from_prior = (sample_mean - mu_0) * arma::trans(sample_mean - mu_0);
  
  // Update the scale hyperparameter
  scale_n = scale_0 + sample_cov + ((lambda * n) / (double)(lambda + n)) * dist_from_prior;
  
  return scale_n;
}

// --- Variance and mean posteriors --------------------------------------------

//' Function to sample a nxn covariance matrix from the inverse wishart for given 
//' priors.
//' 
//' @param nu_0 The prior on the nu parameter of the inverse Wishart (must be
//' greater than or equal to the number of columns in the data).
//' @param scale_0 The prior on the scale parameter of the inverse-Wishart.
//' @param lambda The prior on the shrinkage on variance of the Normal 
//' distribution.
//' @param mu_0 The prior on the mean parameter.
//' @param data The dataset updating our belief of the variance.
//' 
//' @return A covariance matrix sampled from the posterior distribution.
arma::mat sampleVariancePosterior(int nu_0,
                                  arma::mat scale_0,
                                  double lambda,
                                  arma::vec mu_0,
                                  arma::mat data
) {
  
  
  int n = data.n_rows; // The sample size
  int d = data.n_cols; // The dimension of the data
  double df_n = nu_0; // the degrees of freedom
  arma::mat scale_n(d, d); // the scale
  arma::vec sample_mean(d);
  arma::mat sample_cov(d, d);
  arma::mat variance(d, d); // the output variable, the variance matrix

  variance.zeros();
  scale_n = scale_0;
  
  // if any points belong to this data update the priors
  if(n > 0){
    // The posterior degrees of freedom
    df_n += n;
    
    if(n > 1){
      // Calculate the sample mean (the 0 means the mean of the columns is calculated)
      sample_mean = arma::trans(arma::mean(data, 0));
      
    } else {
      // If only one point in cluster the mean is the entry
      sample_mean = arma::trans(data.row(0));
    }
    
    // Sample covariance without division by sample size
    sample_cov = calcSampleCov(data, sample_mean, n, d);
    
    // Calculate the scale variable posterior
    scale_n = updateScale(scale_0, sample_cov, lambda, n, sample_mean, mu_0, d);
  }
  
  // Draw the current variance from the inverse wishart distribution
  variance = arma::iwishrnd(scale_n, df_n);
  
  return variance;
}

//' Function sampling from the mean posterior.
//' 
//' @param mu_0 Prior on the mean parameter.
//' @param variance The covariance matrix for the Multivariate normal distribution (sampled from the inverse-Wishart).
//' @param lambda_0 Prior on the shrinkage applied to the variance parameter.
//' @param data The data used to update our prior.
//' 
//' @return A mean vector sampled from the posterior distribution.
arma::vec sampleMeanPosterior(arma::vec mu_0,
                              arma::mat variance,
                              double lambda_0,
                              arma::mat data
) {
  
  arma::uword n = data.n_rows;
  arma::uword n_col = data.n_cols;
  double lambda_n = lambda_0 + n;
  arma::vec sample_mean = arma::zeros<arma::vec>(n_col);
  arma::vec mu_out = arma::zeros<arma::vec>(n_col);
  arma::vec mu_n = arma::zeros<arma::vec>(n_col);
  arma::mat variance_n = arma::zeros<arma::mat>(n_col, n_col);
  
  
  if (n > 0) {
    sample_mean = trans(arma::mean(data, 0));
  }
  
  // Update mean hyperparameter
  mu_n = updateMu(lambda_0, mu_0, n, n_col, sample_mean);
  
  // Update variance hyperparameter
  variance_n = variance / lambda_n;
  
  // Draw from multivariate Normal distribution
  arma::mat x = arma::mvnrnd(mu_n, variance_n, 1);
  mu_out = x.col(0);
  
  return mu_out;
}


//' Returns a 2-D field of cubes (i.e. a 5D object) for the current iterations 
//' means and variances across each cluster (hence a field)
//' 
//' @param data The entrie dataset.
//' @param cluster_labels Binary vector of cluster labels with the ith entry 
//' corresponding to the membership of the ith observation in data.
//' @param k The maximum number of clusters allowed.
//' @param nu_0 The prior on the nu parameter for the inverse-Wishart.
//' @param num_cols The number of columns in data.
//' @param scale_0 The prior on the sclae parameter. A positive definite matrix.
//' @param mu_0 The prior on the mean parameter.
//' 
//' @return A cube of covariance matrices sampled from the posterior 
//' distribution for each cluster within the data.
arma::cube sampleClusterVariance(arma::mat data,
                                 arma::uvec cluster_labels,
                                 arma::uword k,
                                 int nu_0,
                                 arma::uword num_cols,
                                 arma::mat scale_0,
                                 double lambda_0,
                                 arma::vec mu_0
) {
  
  arma::mat cluster_data; // the data specific to a given cluster
  arma::cube variance(num_cols, num_cols, k); // holds the k cluster-specific variance matrices
  variance.zeros();
  
  // sample the variance for each cluster
  for (arma::uword j = 1; j < k + 1; j++) {

    cluster_data = data.rows(find(cluster_labels == j ));
    
    variance.slice(j - 1) = sampleVariancePosterior(
      nu_0,
      scale_0,
      lambda_0,
      mu_0,
      cluster_data
    );
    
  }
  return variance;
}

//' Returns a matrix of the sampled means for each cluster within the data. The
//' means are column vectors within the outputted matrix.
//' 
//' @param data Matrix of entire dataset.
//' @param cluster_labels Vector of membership labels with entries corresponding 
//' to the data (i.ew. the ith entry is for the observaiton in the ith row of 
//' data).
//' @param k The number of clusters present.
//' @param num_cols The number of columns in the dataset.
//' @param variance A cube of covaraince matrices with that associated with the
//' ith cluster in the ith entry.
//' @param lambda_0 The prior on the shrinkage parameter.
//' @param mu_0 The prior on the mean parameter.
//' 
//' @return A num_cols x k matrix consisting of column vectors of the sampled
//' means for each cluster.
arma::mat sampleClusterMeans(arma::mat data,
                             arma::uvec cluster_labels,
                             arma::uword k,
                             arma::uword num_cols,
                             arma::cube variance,
                             double lambda_0,
                             arma::vec mu_0
) {
  arma::mat cluster_data;
  arma::mat mu(num_cols, k);
  mu.zeros();
  
  // Sample the mean for each cluster
  for (arma::uword j = 1; j < k + 1; j++) {
    cluster_data = data.rows(find(cluster_labels == j ));
    mu.col(j - 1) = sampleMeanPosterior(mu_0, 
           variance.slice(j - 1), 
           lambda_0,
           cluster_data);
    
  }
  return mu;
}

// I think this function is not ever used.
//' Returns a 2-D field of cubes (i.e. a 5D object) for the current iterations 
//' means and variances across each cluster (hence a field)
//' 
//' @param data Matrix of the data.
//' @param cluster_labels Vector of membership labels with entries corresponding 
//' to the data (i.ew. the ith entry is for the observaiton in the ith row of 
//' data).
//' @param k The number of clusters present.
//' @param nu_0 The prior on the nu parameter for the inverse-Wishart.
//' @param num_cols The number of columns in the dataset.
//' @param @param scale_0 The prior on the sclae parameter. A positive definite matrix.
//' @param lambda_0 The prior on the shrinkage parameter.
//' @param mu_0 The prior on the mean parameter.
//' 
//' @return cube of the mean and covariance sampled from the posterior of each 
//' cluster.
arma::field<arma::cube> sampleGaussParams(arma::mat data,
                                          arma::uvec cluster_labels,
                                          arma::uword k,
                                          int nu_0,
                                          arma::uword num_cols,
                                          arma::mat scale_0,
                                          double lambda_0,
                                          arma::vec mu_0
) {
  arma::mat cluster_data;
  arma::mat variance(num_cols, num_cols);
  arma::vec mu(num_cols);
  arma::vec curr_mean(num_cols);
  arma::cube var_entry = arma::zeros<arma::cube>(num_cols, num_cols, k);
  arma::cube mu_entry = arma::zeros<arma::cube>(num_cols, 1, k);
  arma::field<arma::cube> mean_variance_field(2);
  
  mean_variance_field(0) = var_entry;
  mean_variance_field(1) = mu_entry;
  
  curr_mean.zeros();
  
  for (arma::uword j = 1; j < k + 1; j++) {
    cluster_data = data.rows(find(cluster_labels == j ) );
    
    mean_variance_field(0).slice(j - 1) = sampleVariancePosterior(
      nu_0,
      scale_0,
      lambda_0,
      mu_0,
      cluster_data
    );
    
    curr_mean = sampleMeanPosterior(mu_0, 
                                    mean_variance_field(0).slice(j - 1), 
                                    lambda_0,
                                    cluster_data
    );
    
    
    mean_variance_field(1).slice(j - 1) = curr_mean;
  }
  return mean_variance_field;
}