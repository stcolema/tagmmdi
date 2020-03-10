// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// calcEntropy
double calcEntropy(arma::vec class_weights);
RcppExport SEXP _tagmmdi_calcEntropy(SEXP class_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type class_weights(class_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(calcEntropy(class_weights));
    return rcpp_result_gen;
END_RCPP
}
// countCategories
arma::uvec countCategories(arma::umat data);
RcppExport SEXP _tagmmdi_countCategories(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(countCategories(data));
    return rcpp_result_gen;
END_RCPP
}
// declareClassProbsField
arma::field<arma::mat> declareClassProbsField(arma::uvec cat_count, arma::uword n_col, arma::uword n_clust);
RcppExport SEXP _tagmmdi_declareClassProbsField(SEXP cat_countSEXP, SEXP n_colSEXP, SEXP n_clustSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type cat_count(cat_countSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_col(n_colSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_clust(n_clustSEXP);
    rcpp_result_gen = Rcpp::wrap(declareClassProbsField(cat_count, n_col, n_clust));
    return rcpp_result_gen;
END_RCPP
}
// sampleDirichletPosterior
arma::vec sampleDirichletPosterior(arma::vec concentration_0, arma::uvec cluster_labels, arma::uword n_clusters, arma::uword count_from);
RcppExport SEXP _tagmmdi_sampleDirichletPosterior(SEXP concentration_0SEXP, SEXP cluster_labelsSEXP, SEXP n_clustersSEXP, SEXP count_fromSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type concentration_0(concentration_0SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cluster_labels(cluster_labelsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_clusters(n_clustersSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type count_from(count_fromSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleDirichletPosterior(concentration_0, cluster_labels, n_clusters, count_from));
    return rcpp_result_gen;
END_RCPP
}
// sampleCategoryProbabilities
arma::field<arma::mat> sampleCategoryProbabilities(arma::umat data, arma::field<arma::mat> class_probs, arma::field<arma::vec> phi_prior, arma::uvec cluster_labels, arma::uvec cat_count, arma::uword n_clust, arma::uword n_col, arma::uword class_start);
RcppExport SEXP _tagmmdi_sampleCategoryProbabilities(SEXP dataSEXP, SEXP class_probsSEXP, SEXP phi_priorSEXP, SEXP cluster_labelsSEXP, SEXP cat_countSEXP, SEXP n_clustSEXP, SEXP n_colSEXP, SEXP class_startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type class_probs(class_probsSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::vec> >::type phi_prior(phi_priorSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cluster_labels(cluster_labelsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cat_count(cat_countSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_clust(n_clustSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_col(n_colSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type class_start(class_startSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleCategoryProbabilities(data, class_probs, phi_prior, cluster_labels, cat_count, n_clust, n_col, class_start));
    return rcpp_result_gen;
END_RCPP
}
// sampleCategoricalDistn
arma::vec sampleCategoricalDistn(arma::urowvec point, arma::umat data, arma::field<arma::mat> class_probabilities, arma::vec cluster_weights, arma::uword n_clust, arma::uword n_col);
RcppExport SEXP _tagmmdi_sampleCategoricalDistn(SEXP pointSEXP, SEXP dataSEXP, SEXP class_probabilitiesSEXP, SEXP cluster_weightsSEXP, SEXP n_clustSEXP, SEXP n_colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::urowvec >::type point(pointSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type class_probabilities(class_probabilitiesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cluster_weights(cluster_weightsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_clust(n_clustSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_col(n_colSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleCategoricalDistn(point, data, class_probabilities, cluster_weights, n_clust, n_col));
    return rcpp_result_gen;
END_RCPP
}
// categoricalClustering
Rcpp::List categoricalClustering(arma::umat data, arma::field<arma::vec> phi_prior, arma::uvec cluster_labels, arma::uvec fix_vec, arma::vec cluster_weight_priors, arma::uword num_clusters, arma::uword num_iter, arma::uword burn, arma::uword thinning);
RcppExport SEXP _tagmmdi_categoricalClustering(SEXP dataSEXP, SEXP phi_priorSEXP, SEXP cluster_labelsSEXP, SEXP fix_vecSEXP, SEXP cluster_weight_priorsSEXP, SEXP num_clustersSEXP, SEXP num_iterSEXP, SEXP burnSEXP, SEXP thinningSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::vec> >::type phi_prior(phi_priorSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cluster_labels(cluster_labelsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fix_vec(fix_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cluster_weight_priors(cluster_weight_priorsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type num_clusters(num_clustersSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type thinning(thinningSEXP);
    rcpp_result_gen = Rcpp::wrap(categoricalClustering(data, phi_prior, cluster_labels, fix_vec, cluster_weight_priors, num_clusters, num_iter, burn, thinning));
    return rcpp_result_gen;
END_RCPP
}
// gaussianClustering
Rcpp::List gaussianClustering(arma::uword num_iter, arma::vec concentration_0, arma::mat scale_0, arma::uvec class_labels, arma::uvec fix_vec, arma::vec mu_0, double lambda_0, arma::mat data, int nu_0, arma::uword k, arma::uword burn, arma::uword thinning, bool outlier, double t_df, bool normalise, double u, double v);
RcppExport SEXP _tagmmdi_gaussianClustering(SEXP num_iterSEXP, SEXP concentration_0SEXP, SEXP scale_0SEXP, SEXP class_labelsSEXP, SEXP fix_vecSEXP, SEXP mu_0SEXP, SEXP lambda_0SEXP, SEXP dataSEXP, SEXP nu_0SEXP, SEXP kSEXP, SEXP burnSEXP, SEXP thinningSEXP, SEXP outlierSEXP, SEXP t_dfSEXP, SEXP normaliseSEXP, SEXP uSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type concentration_0(concentration_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type scale_0(scale_0SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type class_labels(class_labelsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fix_vec(fix_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_0(mu_0SEXP);
    Rcpp::traits::input_parameter< double >::type lambda_0(lambda_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type nu_0(nu_0SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type thinning(thinningSEXP);
    Rcpp::traits::input_parameter< bool >::type outlier(outlierSEXP);
    Rcpp::traits::input_parameter< double >::type t_df(t_dfSEXP);
    Rcpp::traits::input_parameter< bool >::type normalise(normaliseSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussianClustering(num_iter, concentration_0, scale_0, class_labels, fix_vec, mu_0, lambda_0, data, nu_0, k, burn, thinning, outlier, t_df, normalise, u, v));
    return rcpp_result_gen;
END_RCPP
}
// mdiCatCat
Rcpp::List mdiCatCat(arma::umat data_1, arma::umat data_2, arma::field<arma::vec> class_dist_prior_1, arma::field<arma::vec> class_dist_prior_2, arma::vec clust_weight_priors_1, arma::vec clust_weight_priors_2, arma::uvec clust_labels_1, arma::uvec clust_labels_2, arma::uword n_clust_1, arma::uword n_clust_2, arma::uvec fix_vec_1, arma::uvec fix_vec_2, double a_0, double b_0, arma::uword num_iter, arma::uword burn, arma::uword thinning);
RcppExport SEXP _tagmmdi_mdiCatCat(SEXP data_1SEXP, SEXP data_2SEXP, SEXP class_dist_prior_1SEXP, SEXP class_dist_prior_2SEXP, SEXP clust_weight_priors_1SEXP, SEXP clust_weight_priors_2SEXP, SEXP clust_labels_1SEXP, SEXP clust_labels_2SEXP, SEXP n_clust_1SEXP, SEXP n_clust_2SEXP, SEXP fix_vec_1SEXP, SEXP fix_vec_2SEXP, SEXP a_0SEXP, SEXP b_0SEXP, SEXP num_iterSEXP, SEXP burnSEXP, SEXP thinningSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type data_1(data_1SEXP);
    Rcpp::traits::input_parameter< arma::umat >::type data_2(data_2SEXP);
    Rcpp::traits::input_parameter< arma::field<arma::vec> >::type class_dist_prior_1(class_dist_prior_1SEXP);
    Rcpp::traits::input_parameter< arma::field<arma::vec> >::type class_dist_prior_2(class_dist_prior_2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type clust_weight_priors_1(clust_weight_priors_1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type clust_weight_priors_2(clust_weight_priors_2SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type clust_labels_1(clust_labels_1SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type clust_labels_2(clust_labels_2SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_clust_1(n_clust_1SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_clust_2(n_clust_2SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fix_vec_1(fix_vec_1SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fix_vec_2(fix_vec_2SEXP);
    Rcpp::traits::input_parameter< double >::type a_0(a_0SEXP);
    Rcpp::traits::input_parameter< double >::type b_0(b_0SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type thinning(thinningSEXP);
    rcpp_result_gen = Rcpp::wrap(mdiCatCat(data_1, data_2, class_dist_prior_1, class_dist_prior_2, clust_weight_priors_1, clust_weight_priors_2, clust_labels_1, clust_labels_2, n_clust_1, n_clust_2, fix_vec_1, fix_vec_2, a_0, b_0, num_iter, burn, thinning));
    return rcpp_result_gen;
END_RCPP
}
// mdiGaussCat
Rcpp::List mdiGaussCat(arma::mat cont_data, arma::umat cat_data, arma::vec mu_0, double lambda_0, arma::mat scale_0, int nu_0, double a_0, double b_0, arma::vec clust_wgt_priors_gauss, arma::vec clust_wgt_priors_cat, arma::field<arma::vec> phi_prior, arma::uvec clust_labels_gauss, arma::uvec clust_labels_cat, arma::uword n_clust_gauss, arma::uword n_clust_cat, arma::uvec fix_vec_1, arma::uvec fix_vec_2, arma::uword n_iter, arma::uword burn, arma::uword thin, bool allow_outliers, double t_df, bool normalise, double u_1, double v_1, arma::uword rate_gauss_0, arma::uword rate_cat_0);
RcppExport SEXP _tagmmdi_mdiGaussCat(SEXP cont_dataSEXP, SEXP cat_dataSEXP, SEXP mu_0SEXP, SEXP lambda_0SEXP, SEXP scale_0SEXP, SEXP nu_0SEXP, SEXP a_0SEXP, SEXP b_0SEXP, SEXP clust_wgt_priors_gaussSEXP, SEXP clust_wgt_priors_catSEXP, SEXP phi_priorSEXP, SEXP clust_labels_gaussSEXP, SEXP clust_labels_catSEXP, SEXP n_clust_gaussSEXP, SEXP n_clust_catSEXP, SEXP fix_vec_1SEXP, SEXP fix_vec_2SEXP, SEXP n_iterSEXP, SEXP burnSEXP, SEXP thinSEXP, SEXP allow_outliersSEXP, SEXP t_dfSEXP, SEXP normaliseSEXP, SEXP u_1SEXP, SEXP v_1SEXP, SEXP rate_gauss_0SEXP, SEXP rate_cat_0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type cont_data(cont_dataSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type cat_data(cat_dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_0(mu_0SEXP);
    Rcpp::traits::input_parameter< double >::type lambda_0(lambda_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type scale_0(scale_0SEXP);
    Rcpp::traits::input_parameter< int >::type nu_0(nu_0SEXP);
    Rcpp::traits::input_parameter< double >::type a_0(a_0SEXP);
    Rcpp::traits::input_parameter< double >::type b_0(b_0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type clust_wgt_priors_gauss(clust_wgt_priors_gaussSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type clust_wgt_priors_cat(clust_wgt_priors_catSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::vec> >::type phi_prior(phi_priorSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type clust_labels_gauss(clust_labels_gaussSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type clust_labels_cat(clust_labels_catSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_clust_gauss(n_clust_gaussSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_clust_cat(n_clust_catSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fix_vec_1(fix_vec_1SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fix_vec_2(fix_vec_2SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_iter(n_iterSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< bool >::type allow_outliers(allow_outliersSEXP);
    Rcpp::traits::input_parameter< double >::type t_df(t_dfSEXP);
    Rcpp::traits::input_parameter< bool >::type normalise(normaliseSEXP);
    Rcpp::traits::input_parameter< double >::type u_1(u_1SEXP);
    Rcpp::traits::input_parameter< double >::type v_1(v_1SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type rate_gauss_0(rate_gauss_0SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type rate_cat_0(rate_cat_0SEXP);
    rcpp_result_gen = Rcpp::wrap(mdiGaussCat(cont_data, cat_data, mu_0, lambda_0, scale_0, nu_0, a_0, b_0, clust_wgt_priors_gauss, clust_wgt_priors_cat, phi_prior, clust_labels_gauss, clust_labels_cat, n_clust_gauss, n_clust_cat, fix_vec_1, fix_vec_2, n_iter, burn, thin, allow_outliers, t_df, normalise, u_1, v_1, rate_gauss_0, rate_cat_0));
    return rcpp_result_gen;
END_RCPP
}
// mdiGaussGauss
Rcpp::List mdiGaussGauss(arma::mat data_1, arma::mat data_2, arma::vec mu_0_1, double lambda_0_1, arma::mat scale_0_1, int df_0_1, arma::vec mu_0_2, double lambda_0_2, arma::mat scale_0_2, int df_0_2, arma::vec clust_weight_priors_1, arma::vec clust_weight_priors_2, arma::uvec clust_labels_1, arma::uvec clust_labels_2, arma::uword n_clust_1, arma::uword n_clust_2, arma::uvec fix_vec_1, arma::uvec fix_vec_2, double a_0, double b_0, arma::uword num_iter, arma::uword burn, arma::uword thinning, bool outlier_1, double t_df_1, bool outlier_2, double t_df_2, bool normalise_1, bool normalise_2, double u_1, double v_1, double u_2, double v_2, arma::uword rate_1_0, arma::uword rate_2_0);
RcppExport SEXP _tagmmdi_mdiGaussGauss(SEXP data_1SEXP, SEXP data_2SEXP, SEXP mu_0_1SEXP, SEXP lambda_0_1SEXP, SEXP scale_0_1SEXP, SEXP df_0_1SEXP, SEXP mu_0_2SEXP, SEXP lambda_0_2SEXP, SEXP scale_0_2SEXP, SEXP df_0_2SEXP, SEXP clust_weight_priors_1SEXP, SEXP clust_weight_priors_2SEXP, SEXP clust_labels_1SEXP, SEXP clust_labels_2SEXP, SEXP n_clust_1SEXP, SEXP n_clust_2SEXP, SEXP fix_vec_1SEXP, SEXP fix_vec_2SEXP, SEXP a_0SEXP, SEXP b_0SEXP, SEXP num_iterSEXP, SEXP burnSEXP, SEXP thinningSEXP, SEXP outlier_1SEXP, SEXP t_df_1SEXP, SEXP outlier_2SEXP, SEXP t_df_2SEXP, SEXP normalise_1SEXP, SEXP normalise_2SEXP, SEXP u_1SEXP, SEXP v_1SEXP, SEXP u_2SEXP, SEXP v_2SEXP, SEXP rate_1_0SEXP, SEXP rate_2_0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_1(data_1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data_2(data_2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_0_1(mu_0_1SEXP);
    Rcpp::traits::input_parameter< double >::type lambda_0_1(lambda_0_1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type scale_0_1(scale_0_1SEXP);
    Rcpp::traits::input_parameter< int >::type df_0_1(df_0_1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_0_2(mu_0_2SEXP);
    Rcpp::traits::input_parameter< double >::type lambda_0_2(lambda_0_2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type scale_0_2(scale_0_2SEXP);
    Rcpp::traits::input_parameter< int >::type df_0_2(df_0_2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type clust_weight_priors_1(clust_weight_priors_1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type clust_weight_priors_2(clust_weight_priors_2SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type clust_labels_1(clust_labels_1SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type clust_labels_2(clust_labels_2SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_clust_1(n_clust_1SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n_clust_2(n_clust_2SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fix_vec_1(fix_vec_1SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fix_vec_2(fix_vec_2SEXP);
    Rcpp::traits::input_parameter< double >::type a_0(a_0SEXP);
    Rcpp::traits::input_parameter< double >::type b_0(b_0SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type thinning(thinningSEXP);
    Rcpp::traits::input_parameter< bool >::type outlier_1(outlier_1SEXP);
    Rcpp::traits::input_parameter< double >::type t_df_1(t_df_1SEXP);
    Rcpp::traits::input_parameter< bool >::type outlier_2(outlier_2SEXP);
    Rcpp::traits::input_parameter< double >::type t_df_2(t_df_2SEXP);
    Rcpp::traits::input_parameter< bool >::type normalise_1(normalise_1SEXP);
    Rcpp::traits::input_parameter< bool >::type normalise_2(normalise_2SEXP);
    Rcpp::traits::input_parameter< double >::type u_1(u_1SEXP);
    Rcpp::traits::input_parameter< double >::type v_1(v_1SEXP);
    Rcpp::traits::input_parameter< double >::type u_2(u_2SEXP);
    Rcpp::traits::input_parameter< double >::type v_2(v_2SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type rate_1_0(rate_1_0SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type rate_2_0(rate_2_0SEXP);
    rcpp_result_gen = Rcpp::wrap(mdiGaussGauss(data_1, data_2, mu_0_1, lambda_0_1, scale_0_1, df_0_1, mu_0_2, lambda_0_2, scale_0_2, df_0_2, clust_weight_priors_1, clust_weight_priors_2, clust_labels_1, clust_labels_2, n_clust_1, n_clust_2, fix_vec_1, fix_vec_2, a_0, b_0, num_iter, burn, thinning, outlier_1, t_df_1, outlier_2, t_df_2, normalise_1, normalise_2, u_1, v_1, u_2, v_2, rate_1_0, rate_2_0));
    return rcpp_result_gen;
END_RCPP
}
// createSimilarityMat
arma::mat createSimilarityMat(arma::umat cluster_record);
RcppExport SEXP _tagmmdi_createSimilarityMat(SEXP cluster_recordSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type cluster_record(cluster_recordSEXP);
    rcpp_result_gen = Rcpp::wrap(createSimilarityMat(cluster_record));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _tagmmdi_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _tagmmdi_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _tagmmdi_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _tagmmdi_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tagmmdi_calcEntropy", (DL_FUNC) &_tagmmdi_calcEntropy, 1},
    {"_tagmmdi_countCategories", (DL_FUNC) &_tagmmdi_countCategories, 1},
    {"_tagmmdi_declareClassProbsField", (DL_FUNC) &_tagmmdi_declareClassProbsField, 3},
    {"_tagmmdi_sampleDirichletPosterior", (DL_FUNC) &_tagmmdi_sampleDirichletPosterior, 4},
    {"_tagmmdi_sampleCategoryProbabilities", (DL_FUNC) &_tagmmdi_sampleCategoryProbabilities, 8},
    {"_tagmmdi_sampleCategoricalDistn", (DL_FUNC) &_tagmmdi_sampleCategoricalDistn, 6},
    {"_tagmmdi_categoricalClustering", (DL_FUNC) &_tagmmdi_categoricalClustering, 9},
    {"_tagmmdi_gaussianClustering", (DL_FUNC) &_tagmmdi_gaussianClustering, 17},
    {"_tagmmdi_mdiCatCat", (DL_FUNC) &_tagmmdi_mdiCatCat, 17},
    {"_tagmmdi_mdiGaussCat", (DL_FUNC) &_tagmmdi_mdiGaussCat, 27},
    {"_tagmmdi_mdiGaussGauss", (DL_FUNC) &_tagmmdi_mdiGaussGauss, 35},
    {"_tagmmdi_createSimilarityMat", (DL_FUNC) &_tagmmdi_createSimilarityMat, 1},
    {"_tagmmdi_rcpparma_hello_world", (DL_FUNC) &_tagmmdi_rcpparma_hello_world, 0},
    {"_tagmmdi_rcpparma_outerproduct", (DL_FUNC) &_tagmmdi_rcpparma_outerproduct, 1},
    {"_tagmmdi_rcpparma_innerproduct", (DL_FUNC) &_tagmmdi_rcpparma_innerproduct, 1},
    {"_tagmmdi_rcpparma_bothproducts", (DL_FUNC) &_tagmmdi_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_tagmmdi(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
