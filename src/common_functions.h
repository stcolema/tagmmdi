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

arma::uvec CountCatergories(arma::umat data);

arma::field<arma::mat> DeclareClassProbsField(arma::uvec cat_count,
                                              arma::uword n_col,
                                              arma::uword n_clust
);

arma::field<arma::mat> SampleCategoryProbabilities(arma::umat data,
                                                   arma::field<arma::mat> class_probs,
                                                   arma::field<arma::vec> phi_prior,
                                                   arma::uvec cluster_labels,
                                                   arma::uvec cat_count,
                                                   arma::uword n_clust,
                                                   arma::uword n_col,
                                                   arma::uword class_start
);

arma::vec SampleCategoricalDistn(arma::urowvec point,
                                 arma::umat data,
                                 arma::field<arma::mat> class_probabilities,
                                 arma::vec cluster_weights,
                                 arma::uword n_clust,
                                 arma::uword n_col
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

double CalcNormalLikelihood(arma::vec point,
                            arma::vec mu,
                            arma::mat variance,
                            arma::uword d
);                              


arma::vec SampleGaussianMembership(arma::vec point,
                                   arma::mat data,
                                   arma::uword k,
                                   arma::vec class_weights,
                                   arma::mat mu,
                                   arma::cube variance
);

arma::uword PredictIndex(arma::vec my_vec);

double CalcTdistnLikelihood(arma::vec point,
                            arma::vec mu,
                            arma::mat variance,
                            arma::uword n_col,
                            double df
);

double SampleBetaDistn(double a, double b, double theta = 1.0);

arma::vec CalculateOutlierProb(arma::vec point,
                               arma::vec global_mean,
                               arma::mat global_variance,
                               arma::uword n_col,
                               double t_df,
                               double outlier_weight,
                               arma::vec mu,
                               arma::mat var
);


arma::vec SampleMDIClusterWeights(arma::vec shape_0,
                                  arma::vec rate_0,
                                  double v,
                                  arma::uword n_clust,
                                  arma::uword n_clust_comp,
                                  arma::vec cluster_weights_comp,
                                  arma::uvec cluster_labels,
                                  arma::uvec cluster_labels_comp,
                                  double phi
);

double CalcNormalisingConst(arma::vec cl_wgts_1,
                            arma::vec cl_wgts_2,
                            double phi,
                            arma::uword n_clust_1,
                            arma::uword n_clust_2
);

double SampleMDIPhi(arma::uvec cl_1,
                    arma::uvec cl_2,
                    arma::vec cl_wgts_1,
                    arma::vec cl_wgts_2,
                    double v,
                    arma::uword n,
                    arma::uword min_n_clust,
                    double a_0,
                    double b_0
);

arma::vec SampleMDICatClustProb(arma::uword row_index,
                                arma::umat data,
                                arma::field<arma::mat> class_probs,
                                arma::uword num_clusters,
                                arma::uword num_cols_cat,
                                double phi,
                                arma::vec cluster_weights_categorical,
                                arma::uvec clust_labels,
                                arma::uvec clust_labels_comp
);

arma::vec SampleMDIGaussClustProbs(arma::uword row_index,
                                   arma::mat data,
                                   arma::uword k,
                                   arma::mat mu,
                                   arma::cube variance,
                                   double context_similarity,
                                   arma::vec cluster_weights,
                                   arma::uvec cluster_labels,
                                   arma::uvec cluster_labels_comp
);

arma::vec UpdateClusterLabels(arma::uvec cluster_labels_1,
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
                              double Z
);



#endif