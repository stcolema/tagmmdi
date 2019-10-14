# include <RcppArmadillo.h>
# include <iostream>
# include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;




//' Calculates the t-distribution log likelihood for a given point.
//' 
//' @param point A sample from the dataset.
//' @param mu The mean.
//' @param variance Covariance matrix.
//' @param n_col The length of point.
//' @param df The degrees of freedom variable for a t-distribution.
//' 
//' @return The log-likelihood for a point within a t-distribution.
double CalcTdistnLikelihood(arma::vec point,
                            arma::vec mu,
                            arma::mat variance,
                            arma::uword n_col,
                            double df
){
  
  double log_det = 0.0;
  double exponent = 0.0;
  double log_likelihood = 0.0;
  
  exponent = arma::as_scalar(
    arma::trans(point - mu) 
    * arma::inv(variance)
    * (point - mu)
  );
  
  log_det = arma::log_det(variance).real();
  
  log_likelihood = lgamma((df + n_col) / 2.0) 
    - lgamma(df / 2.0) 
    - (double)n_col / 2.0 * log(df * M_PI) 
    - 0.5 * log_det 
    - ((df + (double) n_col) / 2.0) * log(1.0 + (1.0 / df) * exponent);
    
    return log_likelihood;
}

// See https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions.
//' Samples from a Beta distribution based on the idea that two independent gamma
//' functions can be used. If X ~ Gamma(a, 1) and Y ~ Gamma(b, 1) then
//' X / (X+Y) ~ Beta(a, b). Used in calculating the weights for allocation to 
//' outlier.
//' 
//' @param a Shape parameter of the first gamma distribution.
//' @param b Shape parameter of the second gamma distribution.
//' @param theta (Default = 1.0) Rate of the two Gamma distributions.
//' 
//' @return Sample from Beta(a, b).
double SampleBetaDistn(double a, double b, double theta = 1.0) {
  double X = arma::randg( arma::distr_param(a, 1/theta) );
  double Y = arma::randg( arma::distr_param(b, 1/theta) );
  double beta = X / (double)(X + Y);
  return beta;
}


//' Samples if the point is an outlier or not (comparing to assigned class).
//' 
//' @param point A row from the dataset corresponding to the measurements for a 
//' given sample.
//' @param data Matrix of data.
//' @param outlier_weight Weight for outlier assignment.
//' @param gloabl_mean Sample mean of data.
//' @param global_varinace Sample covariance matrix from data.
//' @param t_df (Default = 4.0) Degrees of freedom parameter for the 
//' t-distribution.
//' 
//' @return Probability that current point is an outlier within data.
double SampleOutlier(arma::vec point,
                     arma::mat data,
                     double outlier_weight,
                     arma::vec global_mean,
                     arma::mat global_variance,
                     double t_df = 4.0
) {
  
  arma::uword n = data.n_rows;
  arma::uword d = data.n_cols;
  double log_likelihood = 0.0;
  double prob = 0.0;
  
  // The probability of belonging to the outlier class
  log_likelihood = t_likelihood(point, global_mean, global_variance, d, t_df);
  
  prob = log(outlier_weight) + log_likelihood;
  
  return prob;
}

//' Calculate the probability of being an outlier or not.
//' @param point A row from the dataset corresponding to the measurements for a 
//' given sample.
//' @param gloabl_mean Sample mean of data.
//' @param global_varinace Sample covariance matrix from data.
//' @param n_col The number of columns present.
//' @param t_df (Default = 4.0) Degrees of freedom parameter for the 
//' t-distribution.
//' @param outlier_weight The weight associated with outlier assignment.
//' @param mu The mean.
//' @param variance Covariance matrix.
//' 
//' @return Vector of non-outlier and outlier assignment probabilities.
arma::vec CalculateOutlierProb(arma::vec point,
                               arma::vec global_mean,
                               arma::mat global_variance,
                               arma::uword n_col,
                               double t_df,
                               double outlier_weight,
                               arma::col mu,
                               arma::mat var
) {
  double out_likelihood = 0.0;
  double in_likelihood = 0.0;
  double non_outlier_weight = 1.0 - outlier_weight;
  arma::vec outlier_prob(2);
  outlier_prob.zeros();
  
  // Calculate outlier likelihood
  out_likelihood = t_likelihood(point, global_mean,  global_variance, n_col, t_df);
  out_likelihood += log(outlier_weight);
  
  outlier_prob(1) = out_likelihood;
  
  // Calculate likelihood of belonging to current cluster
  in_likelihood = normal_likelihood(point, mu, var, n_col);
  
  // Update with the non-outlier weight
  in_likelihood += log(non_outlier_weight);
  
  outlier_prob(0) = in_likelihood;
  
  return outlier_prob;
}