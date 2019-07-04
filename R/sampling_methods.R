#!/usr/bin/env Rscript

# Gaussian, categorical and mdi mixture models

# === Sampling methods =========================================================

#' @title Gibbs sampling
#' @description Carries out gibbs sampling of data and returns a similarity matrix for points
#'
#' @param data A matrix of the data being analysed.
#' @param k The number of clusters.
#' @param class_labels A vector of unsigned integers representing the initial
#' cluster of the corresponding point in data
#' @param num_iter The number of iterations to sample over.
#' @param burn The number of iterations to record after (i.e. the burn-in).
#' @param mu_0 A d-vector; prior on mean. If NULL defaults to mean of data.
#' @param df_0 The prior on the degrees of freedom. if NULL defaults to d + 2.
#' @param scale_0 The prior on the scale for the Inverse Wishart. If NULL
#' generated using an empirical method.
#' @param lambda_0 The prior of shrinkage for mean distribution.
#' @param concentration_0 The prior for dirichlet distribution of cluster
#' weights.
#' @param cat_data Matrix of 1's and 0's used for multiple dataset integration,
#' Default is NULL in which case a generic mixture of Gaussians is used.
#' @param cluster_weight_priors_categorical Vector of the prior on cluster
#' weights in the categorical data.
#' @param phi_0 List of vectors, the prior on the distribution of the classes
#' over clusters.
#' @param c_clusters_label_0 Vector of labels for the prior clustering of the
#' categorical data.
#' @param num_clusters_cat Integer of the number of clusters to have as a
#' maximum in the categorical dataset. Default is 100.
#' @param thinning The step between iterations for which results are recorded in
#' the mcmc output.
#' @param outlier A bool instructing the sampler to consider an additional
#' outlier cluster following a t-distribution
#' @param t_df The degrees of freedom for the outlier t-distribution (default
#' is 4)
#' @param record_posteriors A bool instructing the mcmc function to record the
#' posterior distributions of the mean and variance for each cluster
#' (default is FALSE)
gibbs_sampling <- function(data, k, class_labels,
                           fix_vec = rep(F, nrow(data)),
                           d = ncol(data),
                           N = nrow(data),
                           num_iter = NULL,
                           burn = NULL,
                           thinning = 25,
                           mu_0 = NULL,
                           df_0 = NULL,
                           scale_0 = NULL,
                           lambda_0 = NULL,
                           concentration_0 = NULL,
                           a0 = 1.0,
                           b0 = 0.2,
                           cat_data = NULL,
                           cluster_weight_priors_categorical = 1,
                           phi_0 = NULL,
                           c_clusters_label_0 = NULL,
                           num_clusters_cat = 100,
                           outlier = FALSE,
                           t_df = 4.0,
                           record_posteriors = FALSE,
                           normalise = FALSE) {
  
  # REDUNDANT - DONE IN WRAPPER FUNCTION
  if (is.null(num_iter)) {
    num_iter <- min((d^2) * 1000 / sqrt(N), 10000)
  }

  if (is.null(burn)) {
    burn <- floor(num_iter / 10)
  }

  if (burn > num_iter) {
    stop("Burn in exceeds total iterations. None will be recorded.\nStopping.")
  }

  thinning_warning(thinning, num_iter, burn)

  data <- as.matrix(data)

  # Empirical Bayes
  parameters_0 <- gaussian_arguments(data, k,
    mu_0 = mu_0,
    scale_0 = scale_0,
    lambda_0 = lambda_0,
    df_0 = df_0
  )

  # Extract the arguments from the list
  mu_0 <- parameters_0$mu_0
  df_0 <- parameters_0$df_0
  scale_0 <- parameters_0$scale_0
  lambda_0 <- parameters_0$lambda_0


  # Prior on mass parameter for cluster  weights
  if (!is.null(cat_data)) {
    if (is.null(cluster_weight_priors_categorical)) {
      cluster_weight_priors_categorical <- rep(1, num_clusters_cat)
    } else if (length(cluster_weight_priors_categorical) < num_clusters_cat) {
      cluster_weight_priors_categorical <- rep(
        cluster_weight_priors_categorical,
        num_clusters_cat
      )
    }
  }

  # Gaussian clustering
  sim <- gaussian_clustering(
    num_iter,
    concentration_0,
    scale_0,
    class_labels,
    fix_vec,
    mu_0,
    lambda_0,
    data,
    df_0,
    k,
    burn,
    thinning,
    outlier,
    t_df,
    record_posteriors,
    normalise,
    2,
    10
  )

  sim
}

#' @title Categorical gibbs sampling
#' @description Carries out gibbs sampling of data and returns a similarity matrix for points
#'
#' @param data A matrix of the data being analysed.
#' @param cluster_weight_priors_categorical Vector of the prior on cluster
#' weights in the categorical data.
#' @param phi_0 List of vectors, the prior on the distribution of the classes
#' over clusters.
#' @param c_clusters_label_0 Vector of labels for the prior clustering of the
#' categorical data.
#' @param num_clusters_cat Integer of the number of clusters to have as a
#' maximum in the categorical dataset. Default is 100.
#' @param num_iter The number of iterations to sample over.
#' @param burn The number of iterations to record after (i.e. the burn-in).
#' @param thinning The step between iterations for which results are recorded in
#' the mcmc output.
categorical_gibbs_sampling <- function(data,
                                       fix_vec = rep(F, nrow(data)),
                                       d = ncol(data),
                                       N = nrow(data),
                                       cluster_weight_priors_categorical = 1,
                                       phi_0 = NULL,
                                       c_clusters_label_0 = NULL,
                                       num_clusters_cat = NULL,
                                       num_iter = 10000,
                                       burn = floor(num_iter / 10),
                                       thinning = 25) {
  if (burn > num_iter) {
    stop("Burn in exceeds total iterations. None will be recorded.\nStopping.")
  }

  if (thinning > (num_iter - burn)) {
    if (thinning > (num_iter - burn) & thinning < 5 * (num_iter - burn)) {
      stop("Thinning factor exceeds iterations feasibly recorded. Stopping.")
    } else if (thinning > 5 * (num_iter - burn) & thinning < 10 * (num_iter - burn)) {
      stop("Thinning factor relatively large to effective iterations. Stopping algorithm.")
    } else {
      warning(paste0(
        "Thinning factor relatively large to effective iterations.",
        "\nSome samples recorded. Continuing but please check input"
      ))
    }
  }

  if (is.null(num_clusters_cat)) {
    num_clusters_cat <- min(100, ceiling(nrow(data) / 4))
  }

  data <- as.matrix(data)

  # Empirical Bayes
  if (is.null(cluster_weight_priors_categorical)) {
    cluster_weight_priors_categorical <- rep(1, num_clusters_cat)
  } else if (length(cluster_weight_priors_categorical) < num_clusters_cat) {
    # print(paste0(
    #   "Creating vector of ",
    #   num_clusters_cat,
    #   " repetitions of ",
    #   cluster_weight_priors_categorical,
    #   " for categorical cluster weights prior."
    # ))
    cluster_weight_priors_categorical <- rep(
      cluster_weight_priors_categorical,
      num_clusters_cat
    )
  }

  if (is.null(phi_0)) {
    phi_0 <- phi_prior(cat_data)
  }

  if (is.null(c_clusters_label_0)) {
    c_clusters_label_0 <- sample(1:num_clusters_cat,
      size = N,
      replace = T,
      prob = cluster_weight_priors_categorical
    )
  }

  sim <- categorical_clustering(
    data,
    phi_0,
    c_clusters_label_0,
    fix_vec,
    cluster_weight_priors_categorical,
    num_clusters_cat,
    num_iter,
    burn,
    thinning
  )

  sim
}

#' @title MDI gibbs sampling
#' @description Carries out gibbs sampling of data and returns a similarity
#' matrix for points.
#'
#' @param data_1 A matrix of the data being analysed (context 1).
#' @param data_2 A matrix of the data being analysed (context 2).
#' @param args_1 A named list. The output of either Gaussian arguments or
#' Categorical arguments dpending on the type of data_1, containing the relevant
#' priors.
#' @param args_2 A named list similar to args_1 corrected for the relevant type
#' of data_2.
#' @param type_1 A string, one of "Gaussian" or "Categorical" instructing the
#' function which MDI comparison to use for data_1.
#' @param type_2 A string similar to type_1, one of "Gaussian" or "Categorical"
#' instructing the function which MDI comparison to use for data_2.
#' @param n_clust_1 The number of clusters to use for data_1.
#' @param n_clust_2 The number of clusters to use for data_2.
#' @param labels_0_1 A vector of unsigned integers representing the initial
#' cluster of the corresponding point in data_1.
#' @param labels_0_2 A vector of unsigned integers representing the initial
#' cluster of the corresponding point in data_2.
#' @param d The number of columns in the data (if not input calculates this).
#' @param N The number of observations in the data (defaults to the number of
#' rows in the data).
#' @param a_0 The prior shape for the context similarity parameter.
#' @param b_0 The prior rate for the context similarity parameter.
#' @param fix_vec_1 A vector of 1's and 0's or else bools used for
#' semi-supervised data (use all FALSE or 0 if using an unsupervised case).
#' @param fix_vec_2 A vector of the same type and length as fix_vec but used
#' for the categorical data - defaults to a vector of FALSEs of length N.
#' @param cluster_weight_0_1 The prior for dirichlet distribution of cluster
#' weights.
#' @param cluster_weight_0_2 Vector of the prior on cluster
#' weights in the categorical data.
#' @param num_iter The number of iterations to sample over.
#' @param burn The number of iterations to record after (i.e. the burn-in).
#' @param thinning The step between iterations for which results are recorded in
#' the mcmc output.
#' @param outlier_1 A bool instructing the sampler to consider an additional
#' outlier cluster following a t-distribution for data_1 (only relevant if
#' type_1 is "Gaussian").
#' @param t_df_1 The degrees of freedom for the outlier t-distribution (default
#' is 4).
#' @param outlier_2 A bool instructing the sampler to consider an additional
#' outlier cluster following a t-distribution for data_2 (only relevant if
#' type_2 is "Gaussian").
#' @param t_df_2 The degrees of freedom for the outlier t-distribution (default
#' is 4).
#' @param record_posteriors A bool instructing the mcmc function to record the
#' posterior distributions of the mean and variance for each cluster
#' (default is FALSE)
#' @return A named list of the relevant outputs from the MDI MCMC.
mdi <- function(data_1, data_2,
                args_1 = NULL,
                args_2 = NULL,
                type_1 = "Gaussian",
                type_2 = "Categorical",
                n_clust_1 = 50,
                n_clust_2 = 50,
                labels_0_1 = NULL,
                labels_0_2 = NULL,
                d = ncol(data_1),
                N = nrow(data_1),
                a_0 = 1,
                b_0 = 0.2,
                fix_vec_1 = rep(0, nrow(data_1)),
                fix_vec_2 = rep(0, nrow(data_2)),
                cluster_weight_0_1 = NULL,
                cluster_weight_0_2 = NULL,
                num_iter = NULL,
                burn = floor(num_iter / 10),
                thinning = 25,
                outlier_1 = FALSE,
                t_df_1 = 4.0,
                normalise_1 = FALSE,
                outlier_2 = FALSE,
                t_df_2 = 4.0,
                normalise_2 = FALSE,
                record_posteriors = FALSE,
                save_results = FALSE,
                load_results = FALSE,
                num_load = 0) {

  # Calculate all the relevant parameters
  
  # Check the datasets are of equal height
  if (N != nrow(data_2)) {
    stop("Unequal number of observations in datasets. Incomparable.\nStopping.")
  }

  # Arguments determining the number of MCMC samples generated and recorded
  if (is.null(num_iter)) {
    num_iter <- min((d^2) * 1000 / sqrt(N), 10000)
  }

  # if (is.null(burn)) {
  #   burn <- floor(num_iter / 10)
  # }

  if (burn > num_iter) {
    stop("Burn in exceeds total iterations. None will be recorded.\nStopping.")
  }

  thinning_warning(thinning, num_iter, burn)

  # Declare arguments if not given already
  if (is.null(args_1) & is.null(args_2)) {
    if ((type_1 == "Gaussian" | type_1 == "G")
    & (type_2 == "Categorical" | type_2 == "C")) {
      args_1 <- gaussian_arguments(data_1, n_clust_1)
      args_2 <- categorical_arguments(data_2, n_clust_2)
    } else if ((type_1 == "Gaussian" | type_1 == "G")
    & (type_2 == "Gaussian" | type_2 == "G")) {
      args_1 <- gaussian_arguments(data_1, n_clust_1)
      args_2 <- gaussian_arguments(data_2, n_clust_2)
    } else if ((type_1 == "Categorical" | type_1 == "C")
    & (type_2 == "Categorical" | type_2 == "C")) {
      args_1 <- categorical_arguments(data_1, n_clust_1)
      args_2 <- categorical_arguments(data_2, n_clust_2)
    } else if ((type_1 == "Categorical" | type_1 == "C")
    & (type_2 == "Gaussian" | type_2 == "G")) {
      args_1 <- categorical_arguments(data_1, n_clust_1)
      args_2 <- gaussian_arguments(data_2, n_clust_2)
    }
  }

  # Convert data to matrix format (maybe should just let it hit an error if this
  # isn't done in advance)
  data_1 <- as.matrix(data_1)
  data_2 <- as.matrix(data_2)

  # Declare the cluster weights if not declared in advance
  cluster_weight_0_1 <- declare_cluster_weights(c(0, 0), labels_0_1, n_clust_1,
    weight_0 = cluster_weight_0_1
  )

  cluster_weight_0_2 <- declare_cluster_weights(c(0, 0), labels_0_2, n_clust_2,
    weight_0 = cluster_weight_0_2
  )

  # Generate random labels if not declared in advance
  if (is.null(labels_0_1)) {
    labels_0_1 <- sample(1:n_clust_1,
      size = N,
      replace = T
    )
  }

  if (is.null(labels_0_2)) {
    labels_0_2 <- sample(1:n_clust_2,
      size = N,
      replace = T
    )
  }

  # Call the relevant MDI function
  if ((type_1 == "Gaussian" | type_1 == "G")
  & (type_2 == "Categorical" | type_2 == "C")) {
    mu_0 <- args_1$mu_0
    df_0 <- args_1$df_0
    scale_0 <- args_1$scale_0
    lambda_0 <- args_1$lambda_0

    phi_0 <- args_2$phi

    sim <- mdi_gauss_cat(
      data_1,
      data_2,
      mu_0,
      lambda_0,
      scale_0,
      df_0,
      a_0,
      b_0,
      cluster_weight_0_1,
      cluster_weight_0_2,
      phi_0,
      labels_0_1,
      labels_0_2,
      n_clust_1,
      n_clust_2,
      fix_vec_1,
      fix_vec_2,
      num_iter,
      burn,
      thinning,
      outlier_1,
      t_df_1,
      record_posteriors,
      normalise_1,
      2, # u_1
      10, # v_1
      1, # rate_gauss_0
      1, # rate_cat_0
      save_results,
      load_results,
      num_load
    )
  } else if ((type_1 == "Categorical" | type_1 == "C")
  & (type_2 == "Categorical" | type_2 == "C")) {
    phi_0_1 <- args_1$phi
    phi_0_2 <- args_2$phi

    sim <- mdi_cat_cat(
      data_1,
      data_2,
      phi_0_1,
      phi_0_2,
      cluster_weight_0_1,
      cluster_weight_0_2,
      labels_0_1,
      labels_0_2,
      n_clust_1,
      n_clust_2,
      fix_vec_1,
      fix_vec_2,
      a_0,
      b_0,
      num_iter,
      burn,
      thinning
    )
  } else if ((type_1 == "Gaussian" | type_1 == "G")
  & (type_2 == "Gaussian" | type_2 == "G")) {
    mu_0_1 <- args_1$mu_0
    df_0_1 <- args_1$df_0
    scale_0_1 <- args_1$scale_0
    lambda_0_1 <- args_1$lambda_0

    mu_0_2 <- args_2$mu_0
    df_0_2 <- args_2$df_0
    scale_0_2 <- args_2$scale_0
    lambda_0_2 <- args_2$lambda_0

    sim <- mdi_gauss_gauss(
      data_1,
      data_2,
      mu_0_1,
      lambda_0_1,
      scale_0_1,
      df_0_1,
      mu_0_2,
      lambda_0_2,
      scale_0_2,
      df_0_2,
      cluster_weight_0_1,
      cluster_weight_0_2,
      labels_0_1,
      labels_0_2,
      n_clust_1,
      n_clust_2,
      fix_vec_1,
      fix_vec_2,
      a_0,
      b_0,
      num_iter,
      burn,
      thinning,
      outlier_1,
      t_df_1,
      outlier_2,
      t_df_2,
      record_posteriors,
      normalise_1,
      normalise_2
    )
  }
  sim
}
