#ifndef MODSTRING_H
#define MODSTRING_H

#include <RcppArmadillo.h>
double CalcEntropy(arma::vec class_weights);
arma::mat CreateSimilarityMat(arma::umat cluster_record);

arma::vec CalcConcentrationn(arma::vec concentration_0,
                             arma::uvec cluster_labels,
                             arma::uword n_cat,
                             arma::uword count_from
);

arma::vec SampleDirichletPosterior(arma::vec concentration_0,
                         arma::uvec cluster_labels,
                         arma::uword n_clusters,
                         arma::uword count_from
);

arma::cube SampleClusterVariance(arma::mat data,
                                 arma::uvec cluster_labels,
                                 arma::uword k,
                                 int nu_0,
                                 arma::uword num_cols,
                                 arma::mat scale_0,
                                 double lambda_0,
                                 arma::vec mu_0
);

arma::mat SampleClusterMeans(arma::mat data,
                             arma::uvec cluster_labels,
                             arma::uword k,
                             arma::uword num_cols,
                             arma::cube variance,
                             double lambda_0,
                             arma::vec mu_0
);

double SampleBetaDistn(double a, double b, double theta = 1.0);

SampleGaussianMembership(arma::vec point,
                         arma::mat data,
                         arma::uword k,
                         arma::vec class_weights,
                         arma::mat mu,
                         arma::cube variance
);

arma::uword PredictIndex(arma::vec my_vec);

arma::vec CalculateOutlierProb(arma::vec point,
                               arma::vec global_mean,
                               arma::mat global_variance,
                               arma::uword n_col,
                               double t_df,
                               double outlier_weight,
                               arma::col mu,
                               arma::mat var
);


#endif