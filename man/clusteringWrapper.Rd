% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper.R
\name{clusteringWrapper}
\alias{clusteringWrapper}
\title{Clustering Wrapper}
\usage{
clusteringWrapper(
  MS_object,
  fcol = "markers",
  data_2 = NULL,
  labels_0_1 = NULL,
  args_1 = NULL,
  cluster_weight_0_1 = NULL,
  type_1 = "Gaussian",
  type_2 = "Categorical",
  cluster_weight_0_2 = 1,
  args_2 = NULL,
  labels_0_2 = NULL,
  n_clust_2 = 50,
  fix_vec_1 = NULL,
  fix_vec_2 = NULL,
  a_0 = 1,
  b_0 = 0.2,
  outlier_1 = FALSE,
  t_df_1 = 4,
  normalise_1 = FALSE,
  outlier_2 = FALSE,
  t_df_2 = 4,
  normalise_2 = FALSE,
  train = NULL,
  num_iter = NULL,
  burn = floor(num_iter/10),
  thinning = 25,
  heat_plot = TRUE,
  heat_plot_2 = FALSE,
  main = "heatmap_for_similarity",
  cluster_row = T,
  cluster_cols = T,
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 6,
  entropy_plot = TRUE,
  window_length = min(25, num_iter/5),
  mean_tolerance = 5e-04,
  sd_tolerance = 5e-04,
  sense_check_map = TRUE,
  sense_check_map_2 = FALSE,
  sense_check_main = "component_level_clustering",
  prediction_threshold = 0.5,
  verbose = FALSE
)
}
\arguments{
\item{MS_object}{A dataset in the format used by pRolocdata.}

\item{data_2}{Second dataset for MDI. Either a matrix of continuous data,
categorical data or NULL (in which case MDI is not used).
Default is NULL in which case a generic mixture of Gaussians is used.}

\item{labels_0_1}{An optional prior for clusters in MS_object. If NULL
defaults to a randomly generated set using the proportion in the labelled
data.}

\item{args_1}{A named list containing the relevant objects for the data type
of MS_object. If NULL this is generated based on the input type_1 and using
the function gaussian_arguments or categorical arguments.}

\item{cluster_weight_0_1}{The prior for dirichlet distribution of cluster
weights.}

\item{type_1}{String denoting the type of data in MS_object. One of
"Gaussian" or "Categorical" (shorthand "G" and "C" also suffice). Default is
"Gaussian".}

\item{type_2}{String denoting the type of data in MS_object. One of
"Gaussian" or "Categorical" (shorthand "G" and "C" also suffice). Default is
"Categorical".}

\item{cluster_weight_0_2}{Vector of the prior on cluster
weights in the categorical data. If a single number is given a uniform vector
is generated.}

\item{args_2}{A named list containing the relevant objects for the data type
of data_2. If NULL this is generated based on the input type_1 and using
the function gaussian_arguments or categorical arguments.}

\item{labels_0_2}{Vector of labels for the prior clustering of the
categorical data.}

\item{n_clust_2}{Integer of the number of clusters to have as a
maximum in the categorical dataset. Default is 100.}

\item{fix_vec_1}{A vector of 0's and 1's indicating which points are fixed in
their initial allocation in dataset 1. Default of NULL means that the
information in MS_object is used to define this (any labelled points in the
dataset are held in this allocation).}

\item{fix_vec_2}{A vector of 0's and 1's indicating which points are fixed in
their initial allocation for dataset 2. Default of all points unfixed.}

\item{a_0}{Prior shape parameter for gamma distribution of strategic latent
variable, v, in MDI. Default is 1.}

\item{b_0}{Prior rate parameter for gamma distribution of strategic latent
variable, v, in MDI. Default is 0.2.}

\item{outlier_1}{A bool instructing the sampler to consider an outlier
cluster following a t-distribution (i.e. implementing TAGM). Only
applicable to gaussian datasets.}

\item{t_df_1}{The degrees of freedom for the outlier t-distribution (default
is 4).}

\item{normalise_1}{Bool instructing normalisation of continuous data in
dataset 1 (default is false).}

\item{outlier_2}{A bool instructing the sampler to consider an outlier
cluster following a t-distribution (i.e. implementing TAGM). Only
applicable to gaussian datasets.}

\item{t_df_2}{The degrees of freedom for the outlier t-distribution (default
is 4).}

\item{normalise_2}{Bool instructing normalisation of continuous data in
dataset 2 (default is false).}

\item{train}{instruction to include all data (NULL), labelled data (TRUE) or
unlabelled data (FALSE). Default is NULL.}

\item{num_iter}{The number of iterations to sample over.}

\item{burn}{The number of iterations to record after (i.e. the burn-in).}

\item{thinning}{The step between iterations for which results are recorded in
the mcmc output.}

\item{heat_plot}{A bool. Instructs saving and printing of heatmap of
similarity matrix from Gibbs sampling. Default is TRUE.}

\item{main}{String. The title for heatmap, default is "heatmap_for_similarity".}

\item{cluster_row}{A bool. Instructs pheatmap to cluster rows using a tree.}

\item{cluster_cols}{A bool. instructs pheatmap to cluster columns using a
tree.}

\item{fontsize}{The size of font in pheatmap.}

\item{fontsize_row}{The fontsize in rows in pheatmap.}

\item{fontsize_col}{The fontsize in columns in pheatmap.}

\item{entropy_plot}{A bool instructing function to save a plot of the entropy
across all iterations. Default is TRUE.}

\item{window_length}{A number. Input to entropy ploy function. Default is
min(25, num_iter / 5).}

\item{mean_tolerance}{Input to entropy ploy function. Default is 0.0005.}

\item{sd_tolerance}{Input to entropy ploy function. Default is 0.0005.}

\item{sense_check_map}{A bool instructing function to save a heatmap of the
componenet level clusterings.}

\item{sense_check_main}{String. Title for sense_check_map. Default is
"component_level_clustering".}

\item{prediction_threshold}{The minimum proportion of recorded iterations
for which a point is in its most common cluster for which a prediction is
returned. If below this predicted class is NA.}
}
\value{
A named list including at least the output from the gibbs sampler,
but can include two pheatmaps and a scatter plot of the entropy over
iterations.
}
\description{
Returns mean, variance and similarity posteriors from Gibbs
sampling with option of pheatmap
}
\examples{
# Single dataset clustering
data("hyperLOPIT2015") # MS object from pRolocData
mcmc_object <- clusteringWrapper(hyperLOPIT2015,
  num_iter = 10000,
  burn = 1000,
  thinning = 50,
  outlier_1 = TRUE,
  heat_plot = TRUE,
  main = "Gene clustering by organelle",
  prediction_threshold = 0.5
)

# Use some categorical data
cat_data <- as.matrix(exprs(tan2009r1goCC))

# Implement MDI
stuff <- clusteringWrapper(tan2009r1,
  data_2 = cat_data,
  num_iter = 10,
  burn = 1,
  thinning = 1,
  outlier_1 = TRUE,
  heat_plot = T,
  main = "Gene clustering by organelle",
  prediction_threshold = 0.4,
  sense_check_map = F
)
}
