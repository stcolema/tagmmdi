#!/usr/bin/env Rscript

# === Functions ================================================================

# --- MCMC analysis ------------------------------------------------------------
#' @title entropy_window
#' @description  Find the point at which entropy stabilises in the iterations.
#'
#' @param entropy_vec A vector of numbers corresponding to entropy of each
#' iteration.
#' @param start An integer instructing which iteration to start with (default is
#' 1).
#' @param window_length The number of iterations to consider when considering
#' convergence (default is 25).
#' @param mean_tolerance A number. The threshold for how close the mean of the
#' two windows must be to be considered converged (default is 0.001).
#' @param sd_tolerance: A number. The threshold for how close the standard
#' deviation of the two windows must be to be considered converged (default is
#' 0.001).
#' @return The iteration at which convergence occurs in the clustering
entropy_window <- function(entropy_vec,
                           start = 1,
                           window_length = 25,
                           mean_tolerance = 0.001,
                           sd_tolerance = 0.001) {
  n <- length(entropy_vec)

  search_range <- seq(
    from = start,
    to = n - window_length,
    by = window_length
  )

  for (i in search_range) {

    # Create two windows looking forward from the current iteration and compare
    # their means and standard deviations
    win_1 <- entropy_vec[i:(i + window_length - 1)]
    win_2 <- entropy_vec[(i + window_length)
    :min((i + 2 * window_length - 1), n)]

    mean_1 <- mean(win_1)
    mean_2 <- mean(win_2)

    sd_1 <- sd(win_1)
    sd_2 <- sd(win_2)

    # If the differences are less than the predefined tolerances, return this
    # iteration as the point to burn up to
    if ((abs(mean_1 - mean_2) < mean_tolerance)
    & (abs(sd_1 - sd_2) < sd_tolerance)
    ) {
      return(i)
    }
  }
}

# === PRE-PROCESSING ===========================================================

# Functions to prepare the data and various parameters for the sampling

# --- Data preparation ---------------------------------------------------------

#' @title MS dataset
#' @description Given an MS object from pRolocData, returns the numerical data,
#' the count of each class and the vector recording which labels are fixed (i.e
#' from the training set) and unfixed.
#'
#' @param MS_object A dataset in the format used by pRolocdata.
#' @param train: instruction to include all data (NULL), labelled data (TRUE) or
#' unlabelled data (FALSE). Default is NULL.
#' @examples
#' data("hyperLOPIT2015") # MS object from pRolocData
#' MS_data <- MS_dataset(hyperLOPIT2015)
MS_dataset <- function(MS_object, train = NULL) {
  # Data with labels
  mydata_labels <- pRoloc:::subsetAsDataFrame(
    object = MS_object,
    fcol = "markers",
    train = TRUE
  )

  # fixed <- rep(TRUE, nrow(mydata_labels))
  fixed <- rep(1, nrow(mydata_labels))

  mydata_no_labels <- pRoloc:::subsetAsDataFrame(
    object = MS_object,
    fcol = "markers",
    train = FALSE
  )

  # not_fixed <- rep(FALSE, nrow(mydata_no_labels))
  not_fixed <- rep(0, nrow(mydata_no_labels))

  nk <- tabulate(fData(markerMSnSet(MS_object))[, "markers"])

  mydata_no_labels$markers <- NA

  if (is.null(train)) {
    row_names <- c(rownames(mydata_labels), rownames(mydata_no_labels))
    mydata <- suppressWarnings(bind_rows(mydata_labels, mydata_no_labels))
    fix_vec <- c(fixed, not_fixed)
  } else if (isTRUE(train)) {
    row_names <- c(rownames(mydata_labels))
    mydata <- mydata_labels
    fix_vec <- fixed
  } else {
    row_names <- c(rownames(mydata_no_labels))
    mydata <- mydata_no_labels
    fix_vec <- not_fixed
  }
  return(list(data = mydata, fix_vec = fix_vec, row_names = row_names, nk = nk))
}

#' @title Prepare cat data
#' @description Converts categorical data to numerical format appropriate for
#' analysis (i.e. integers with a null class of 0)
#'
#' @param data  A data frame or matrix of categorical data
#' @return A matrix of integers where the lowest class fval
prepare_cat_data <- function(data) {
  data <- data %>%
    as.data.frame() %>%
    lapply(as.character) %>%
    lapply(as.factor) %>%
    lapply(as.numeric) %>%
    do.call(rbind, .) %>%
    "-"(1) %>%
    t()
}

# --- Parameters  --------------------------------------------------------------

#' @title Thinning warning
#' @description Generates an error or warning if the relative values of
#' thinning, num_iter and burn (the arguments controlling the number of recorded
#' MCMC samples) are strange.
#' @return None.
thinning_warning <- function(thinning, num_iter, burn) {
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
}

#' Generates priors for the mean, degrees of freedom and scale parameters if not
#' set.
#'
#' @param data A matrix of the data being analysed.
#' @param mu_0 A d-vector. If NULL defaults to a vector of column means of data.
#' @param df_0 An integer. If NULL defaults to d + 2.
#' @param scale_0 A positive definite matrix. The prior for the scale parameter
#' of the inverse wishart distribution; if NULL defaults to a diagonal matrix.
#' @param N The number of entries in data.
#' @param k The number of clusters used.
#' @param d The number of columns in data.
#' @param lambda_0 A positive real number; the shrinkage prior for the mean.
#' @return A named list of the three hyperparameters, mean, scale and degrees of
#'  freedom
empirical_bayes_gaussian <- function(data, mu_0, df_0, scale_0, N, k, d,
                                     lambda_0 = 0.01) {
  parameters <- list()
  if (is.null(mu_0)) {
    mu_0 <- colMeans(data)
  }

  if (is.null(df_0)) {
    df_0 <- d + 2
  }

  if (is.null(scale_0)) {
    # scale_0 <- diag(colSums((data - mean(data))^2) / N) / (k^(1 / d))
    # if (any(is.na(scale_0))) {
    scale_0 <- diag(d) / (k^(1 / d))
    # }
  }
  parameters$mu_0 <- mu_0
  parameters$df_0 <- df_0
  parameters$scale_0 <- scale_0
  parameters$lambda_0 <- lambda_0

  return(parameters)
}

#' @title Phi prior
#' @description Generates a prior for the phi vector for each variable for the Dirichlet
#' distribution
#' @param matrix_data A matrix of data.
#' @return A list of vectors of the proportion of each level across all of
#' matrix_data.
phi_prior <- function(matrix_data) {

  # lambda function applies ``table'' to each column of matr_data before
  # removing names. This ensures the same output type for the case when all
  # variables have the same number of levels and when they do not (rather than
  # merely using one or two lapply's of table)
  unnamed_list_prop <- lapply(
    1:ncol(matrix_data),
    function(i) {
      table(matrix_data[, i])
    }
  ) %>%
    unname() %>%
    lapply("/", nrow(matrix_data))
}


#' @title Gaussian arguments
#' @description Creates a named list of the priors required for a mixture of
#' Gaussians. If the user provides values the function will wrap them in a
#' single object appropriate for input into the MDI function. If values are not
#' given the priors are generated using an empirical method (hence the need for
#' inputting the data and the number of clusters).
#' @param data A matrix of continuous values.
#' @param n_clust A unsigned integer; the number of clusters to be used in the
#' clustering method.
#' @param mu_0 A d-vector where d is the number of columns in data representing
#' the prior beliefs about the mean for data. Defaults to the mean of data if
#' given value is NULL.
#' @param scale_0 A positive definite matrix. The prior for the scale parameter
#' of the inverse wishart distribution; if NULL defaults to a diagonal matrix.
#' @param lambda_0 A positive real number; the shrinkage prior for the mean.
#' @param df_0 An integer. The degrees of freedom used in the inverse Wishart
#' distribution from which the variance is sampled. If NULL defaults to d + 2.
#' @return A named list ready to be used as input into the MDI function.
#' @examples
#' args <- gaussian_arguments(data, 5)
gaussian_arguments <- function(data, n_clust,
                               mu_0 = NULL,
                               scale_0 = NULL,
                               lambda_0 = NULL,
                               df_0 = NULL) {
  N <- nrow(data)
  d <- ncol(data)

  if (is.null(lambda_0)) {
    lambda_0 <- 0.01
  }

  args <- empirical_bayes_gaussian(data, mu_0, df_0, scale_0, N, n_clust, d,
    lambda_0 = lambda_0
  )

  args
}

#' @title Categorical arguments
#' @description Creates a named list of the priors required for a mixture of
#' Dirichlets. The priors are generated using an empirical method (hence the
#' need for inputting the data and the number of clusters).
#' @param data A matrix of continuous values.
#' @param n_clust A unsigned integer; the number of clusters to be used in the
#' clustering method.
#' @return A named list ready to be used as input into the MDI function.
categorical_arguments <- function(data, n_clust) {
  args <- list(phi = phi_prior(data))
  args
}

#' @title Cluster weight prior
#' @description Produces a vector of prior weights for clusters
#'
cluster_weight_prior <- function(train_data, outlier = FALSE) {
  weights <- train_data %>%
    table() %>%
    unname() %>%
    lapply("/", nrow(.)) %>%
    unlist()

  if (outlier) {
    weights <- c(weights, mean(weights))
  }

  weights
}

#' @title Declare cluster weights
#' @description Creates a vector for the mass parameter for cluster weights.
#' @param fix_vec Vector of 1's and 0's indicating if the corresponding entry in
#' the data is a fixed label
#' @param clust_labels A vector of unsigned integers representing the current
#' labelling of points in the data
#' @param n_clust A unsigned integer; the number of clusters to be used in the
#' clustering method.
#' @param weight_0 One of an int or a p-vector where p = n_clust. If an int
#' defaults to a p-vector of p repetitions of said int.
#' @return A vector of values to be used as the mass parameters for the weights
#' for the clusters.
declare_cluster_weights <- function(fix_vec, clust_labels, n_clust,
                                    weight_0 = NULL) {
  if (any(fix_vec == 1) & is.null(weight_0)) {
    relevant_labels <- clust_labels[fix_vec == 1]
    weight_0 <- unname(table(relevant_labels) / sum(table(relevant_labels)))
  } else {
    if (is.null(weight_0)) {
      weight_0 <- 1
    }
    if (length(weight_0) < n_clust) {
      print(paste0(
        "Creating vector of ", n_clust, " repetitions of ", weight_0,
        " for mass parameter prior."
      ))
      weight_0 <- rep(weight_0, n_clust)
    }
  }
  weight_0
}

#' @title Cluster label prior
#' @description Generates a vector of labels if required
#' @param class_labels_0 An optional prior for clusters in MS_object. If NULL
#' defaults to a randomly generated set using the proportion in the labelled
#' data.
#' @param nk Table of the frequency of the classes in the datset
#' @param train A NULL, TRUE or FALSE value informing if supervised (TRUE),
#' semi-supervised (NULL) or unsupervised (FALSE)
#' @param MS_object The MS object from pRolocData being analysed
#' @param k The number of clusters in the dataset
#' @param N The number of observations
#' @return A vector of integers corresponding to the cluster allocation of the N
#' observations
cluster_label_prior <- function(class_labels_0,
                                nk,
                                train,
                                MS_object,
                                k,
                                N) {
  # Generate class labels
  if (is.null(class_labels_0)) {
    class_weights <- nk / sum(nk)
    if (is.null(train)) {
      fixed_labels <- as.numeric(fData(markerMSnSet(MS_object))[, "markers"])
      class_labels_0 <- c(fixed_labels, sample(1:k, N - length(fixed_labels),
        replace = T # ,
        # prob = class_weights # Remove this as not a good assumption
      ))
    } else if (isTRUE(train)) {
      class_labels_0 <- as.numeric(fData(markerMSnSet(MS_object))[, "markers"])
    } else {
      class_labels_0 <- sample(1:k, N, replace = T) # , prob = class_weights)
    }
  }
  return(class_labels_0)
}

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
                     df_0 = df_0)
  
  # parameters_0 <- empirical_bayes_gaussian(data, mu_0, df_0, scale_0, N, k, d)

  mu_0 <- parameters_0$mu_0
  df_0 <- parameters_0$df_0
  scale_0 <- parameters_0$scale_0
  lambda_0  <- parameters_0$lambda_0


  # Prior on mass parameter for cluster  weights
  if (!is.null(cat_data)) {
    if (is.null(cluster_weight_priors_categorical)) {
      cluster_weight_priors_categorical <- rep(1, num_clusters_cat)
    } else if (length(cluster_weight_priors_categorical) < num_clusters_cat) {
      print(paste0(
        "Creating vector of ", num_clusters_cat, " repetitions of ", cluster_weight_priors_categorical,
        " for categorical cluster weights prior."
      ))
      cluster_weight_priors_categorical <- rep(cluster_weight_priors_categorical, num_clusters_cat)
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
    print(paste0(
      "Creating vector of ", num_clusters_cat, " repetitions of ", cluster_weight_priors_categorical,
      " for categorical cluster weights prior."
    ))
    cluster_weight_priors_categorical <- rep(cluster_weight_priors_categorical, num_clusters_cat)
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
#' @param data A matrix of the data being analysed.
#' @param cat_data Matrix of 1's and 0's used for multiple dataset integration.
#' @param k The number of clusters.
#' @param class_labels A vector of unsigned integers representing the initial
#' cluster of the corresponding point in data.
#' @param fix_vec A vector of 1's and 0's or else bools used for semi-supervised
#' data (use all FALSE or 0 if using an unsupervised case).
#' @param d The number of columns in the data (if not input calculates this).
#' @param N The number of observations in the data (defaults to the number of
#' rows in the data).
#' @param fix_vec_cat A vector of the same type and length as fix_vec but used
#' for the categorical data - defaults to a vector of FALSEs of length N.
#' @param num_iter The number of iterations to sample over.
#' @param burn The number of iterations to record after (i.e. the burn-in).
#' @param mu_0 A d-vector; prior on mean. If NULL defaults to mean of data.
#' @param df_0 The prior on the degrees of freedom. if NULL defaults to d + 2.
#' @param scale_0 The prior on the scale for the Inverse Wishart. If NULL
#' generated using an empirical method.
#' @param lambda_0 The prior of shrinkage for mean distribution.
#' @param concentration_0 The prior for dirichlet distribution of cluster
#' weights.
#' @param a_0 The prior shape for the context similarity parameter
#' @param b_0 The prior rate for the context similarity parameter
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
mdi_gauss_cat_clustering <- function(data, cat_data, k, class_labels, fix_vec,
                                     d = NULL,
                                     N = NULL,
                                     fix_vec_cat = NULL,
                                     num_iter = NULL,
                                     burn = NULL,
                                     thinning = 25,
                                     mu_0 = NULL,
                                     df_0 = NULL,
                                     scale_0 = NULL,
                                     lambda_0 = 0.01,
                                     concentration_0 = 0.1,
                                     a_0 = 1,
                                     b_0 = 1,
                                     cluster_weight_priors_categorical = 1,
                                     phi_0 = NULL,
                                     c_clusters_label_0 = NULL,
                                     num_clusters_cat = NULL,
                                     outlier = FALSE,
                                     t_df = 4.0,
                                     record_posteriors = FALSE,
                                     normalise = FALSE) {
  if (is.null(d)) {
    d <- ncol(data)
  }

  if (is.null(N)) {
    N <- nrow(data)
  }

  if (is.null(fix_vec_cat)) {
    fix_vec_cat <- rep(F, N)
  }


  if (is.null(num_iter)) {
    num_iter <- min((d^2) * 1000 / sqrt(N), 10000)
  }

  if (is.null(burn)) {
    burn <- floor(num_iter / 10)
  }

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
  parameters_0 <- empirical_bayes_gaussian(data, mu_0, df_0, scale_0, N, k, d)

  mu_0 <- parameters_0$mu_0
  df_0 <- parameters_0$df_0
  scale_0 <- parameters_0$scale_0

  if (is.null(concentration_0)) {
    concentration_0 <- rep(0.1, (k + outlier))
  } else if (length(concentration_0) < (k + outlier)) {
    print(paste0(
      "Creating vector of ", k + outlier, " repetitions of ", concentration_0,
      " for concentration prior."
    ))
    concentration_0 <- rep(concentration_0, k + outlier)
  }

  if (is.null(cluster_weight_priors_categorical)) {
    cluster_weight_priors_categorical <- rep(1, num_clusters_cat)
  } else if (length(cluster_weight_priors_categorical) < num_clusters_cat) {
    print(paste0(
      "Creating vector of ", num_clusters_cat, " repetitions of ", cluster_weight_priors_categorical,
      " for categorical cluster weights prior."
    ))
    cluster_weight_priors_categorical <- rep(cluster_weight_priors_categorical, num_clusters_cat)
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

  sim <- mdi_gauss_cat(
    data,
    cat_data,
    mu_0,
    lambda_0,
    scale_0,
    df_0,
    a_0,
    b_0,
    concentration_0,
    cluster_weight_priors_categorical,
    phi_0,
    class_labels,
    c_clusters_label_0,
    k,
    num_clusters_cat,
    fix_vec,
    fix_vec_cat,
    num_iter,
    burn,
    thinning,
    outlier,
    t_df,
    record_posteriors,
    normalise
  )
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
                record_posteriors = FALSE) {

  # Calculate all the relevant parameters
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
  # cluster_weight_0_1 <- declare_cluster_weights(fix_vec_1, labels_0_1, n_clust_1,
  #   weight_0 = cluster_weight_0_1
  # )

  cluster_weight_0_1 <- declare_cluster_weights(c(0, 0), labels_0_1, n_clust_1,
    weight_0 = cluster_weight_0_1
  )

  # cluster_weight_0_2 <- declare_cluster_weights(fix_vec_2, labels_0_2, n_clust_2,
  #   weight_0 = cluster_weight_0_2
  # )

  cluster_weight_0_2 <- declare_cluster_weights(c(0, 0), labels_0_2, n_clust_2,
    weight_0 = cluster_weight_0_2
  )

  # Generate random labels if not declared in advance
  if (is.null(labels_0_1)) {
    labels_0_1 <- sample(1:n_clust_1,
      size = N,
      replace = T #,
      # prob = cluster_weight_0_1
    )
  }

  if (is.null(labels_0_2)) {
    labels_0_2 <- sample(1:n_clust_2,
      size = N,
      replace = T #,
      # prob = cluster_weight_0_2
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
      normalise_1
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






# === OUTPUT ===================================================================

# --- Heatmap ------------------------------------------------------------------
#' @title Annotated Heatmap
#' @description Returns and prints an annotated pheatmap
#'
#' @param input_data A matrix of data to be heat mapped. Needs column and
#' rownames.
#' @param annotation_row A data frame of the annotation variable(s). Names must
#' match with input_data. If NULL returns normal pheatmap.
#' @param sort_by_col The name of the column to sort input_data and
#' annotation_row by when heatmapping. If pheatmap is instructed to sort_rows
#' this has no impact. Default is NULL in which case no sorting occurs.
#' @param train If FALSE returns normal pheatmap.
#' @param ... The usual inputs for pheatmap.
#' @return An annotated pheatmap from the pheatmap package
annotated_heatmap <- function(input_data, annotation_row = NULL,
                              sort_by_col = NULL,
                              train = NULL,
                              ...) {
  if (is.null(annotation_row) & !(isTRUE(train) | is.null(train))) {
    stop("If data")
  }
  dissim <- input_data

  sort_by_col <- attempt::try_catch(
    expr = enquo(sort_by_col),
    .e = NULL,
    .w = NULL
  )

  # print(sort_by_col)

  if (sort_by_col != quo(NULL)) {
    # print("HIII")

    # sort_col <- enquo(sort_by_col)

    # Declare row names ro ensure re-ordered correctly
    row_names <- data.frame(Names = row.names(annotation_row))

    # Select the column of interest
    col_of_interest <- annotation_row %>%
      dplyr::select(!!sort_by_col)

    if (!is.null(annotation_row)) {
      # Combine the datasets to ensure all have common ordering
      combined_data <- bind_cols(dissim, annotation_row, row_names)

      # Arrange based on the user selected columm
      sorted_data <- combined_data %>%
        arrange(!!sort_by_col)

      # Separate out the dataframes
      annotation_row <- sorted_data %>%
        dplyr::select(one_of(names(annotation_row)))

      row_names <- combined_data %>%
        dplyr::select(Names)

      dissim <- sorted_data %>%
        dplyr::select(-one_of(names(annotation_row))) %>%
        dplyr::select(-Names)

      # Redeclare row names (as dplyr strips these)
      row.names(annotation_row) <- row_names$Names
      row.names(dissim) <- row_names$Names
    } else {
      combined_data <- bind_cols(dissim, row_names)

      # Arrange based on the user selected columm
      sorted_data <- combined_data %>%
        arrange(!!sort_by_col)

      row_names <- combined_data %>%
        dplyr::select(Names)

      dissim <- sorted_data %>%
        dplyr::select(-Names)

      # Redeclare row names (as dplyr strips these)
      row.names(dissim) <- row_names$Names
    }

    # rownames(dissim) <- rownames(input_data)
  }

  if (!is.null(annotation_row)) {
    rownames(annotation_row) <- rownames(input_data)
    rownames(dissim) <- rownames(input_data)
  }

  # Colour scheme for heatmap
  col_pal <- RColorBrewer::brewer.pal(9, "Blues")

  if (!is.null(annotation_row) | !(is.null(train) | isTRUE(train))) {
    feature_names <- names(annotation_row)

    # Annotation colours
    new_cols_list <- list()
    my_colours <- list()
    for (feature in feature_names) {
      outlier_present <- FALSE

      # print(feature)
      types_feature_present <- unique(annotation_row[[feature]][!is.na(annotation_row[[feature]])])

      # print(types_feature_present)

      if (feature == "Predicted_class") {
        if ("Outlier" %in% types_feature_present) {
          outlier_present <- TRUE
          types_feature_present <- types_feature_present[types_feature_present != "Outlier"]
        }
      }

      # features_present[[feature]] <- types_feature_present
      new_cols_list[[feature]] <- colorRampPalette(grDevices::rainbow(length(types_feature_present)))



      my_colours[[feature]] <- new_cols_list[[feature]](length(types_feature_present))

      if (outlier_present) {
        my_colours[[feature]] <- c(my_colours[[feature]], "black")
        types_feature_present <- c(types_feature_present, "Outlier")
      }

      names(my_colours[[feature]]) <- types_feature_present
    }

    # Heatmap
    if (is.null(train) | isTRUE(train)) {
      heat_map <- pheatmap(dissim,
        annotation_row = annotation_row,
        annotation_colors = my_colours,
        ...
      )
    }
  }
  else {
    heat_map <- pheatmap(
      dissim,
      ...
    )
  }
  return(heat_map)
}


#' @title Pheatmap cluster by col
#' @description Returns and prints an annotated pheatmap sorted by a given
#' column of the annotation data frame with clustering occuring within the
#' levels of said column
#'
#' @param num_data A matrix of data to be heat mapped. Needs column and
#' rownames.
#' @param annotation_row A data frame of the annotation variable(s). Names must
#' match with input_data. If NULL returns normal pheatmap.
#' @param sort_col The name of the column to sort input_data and
#' annotation_row by when heatmapping. If pheatmap is instructed to sort_rows
#' this has no impact.
#' @param main The default title for the heatmap
#' @param ... The usual inputs for pheatmap.
#' @return An annotated pheatmap from the pheatmap package

pheatmap_cluster_by_col <- function(num_data, annotation_row, sort_col,
                                    main = "sense_check",
                                    use_col_gaps = TRUE,
                                    ...) {

  # save row names as dplyr removes them
  row_names <- row.names(num_data)

  # Enclose the column name using tidy evaluation
  sort_col <- enquo(sort_col)

  # select the variable of interest for sorting
  if (ncol(annotation_row) > 1) {
    col_of_interest <- annotation_row %>%
      dplyr::select(!!sort_col)
  } else {
    col_of_interest <- annotation_row
  }

  # Create new order
  new_order <- order(col_of_interest)

  # Apply new order to annotation row and the data
  num_data <- num_data[new_order, ]
  annotation_row <- annotation_row[new_order, ]
  row_names <- row_names[new_order]

  if (!is.data.frame(annotation_row)) {
    annotation_row <- as.data.frame(annotation_row)
  }

  # print(col_of_interest)

  # Arrange the annotation data frame based on the sort column
  # annotation_row <- annotation_row %>%
  # dplyr::arrange(!!sort_col)

  # Select the sort column to find the location for the gaps

  if (ncol(annotation_row) > 1) {
    col_of_interest <- annotation_row %>%
      dplyr::select(!!sort_col)
  } else {
    col_of_interest <- annotation_row
  }

  # Col of interest has row names, remove these and replace with numbers
  row.names(col_of_interest) <- 1:nrow(col_of_interest)

  # Create the gaps for the heatmap for between certain rows
  gaps <- as.numeric(row.names(unique(col_of_interest))) - 1

  # Create an empty vector to hold the new ordering
  ordering <- c()

  # Create a new gap variable including the end point of the data frame
  loc_gapping <- c(gaps, nrow(num_data))

  # Iterate over the data frame clustering the points between gaps
  for (i in 2:length(loc_gapping)) {
    # Cluster the points in the current gap
    hc <- hclust(dist(num_data[(loc_gapping[i - 1] + 1):loc_gapping[i], ])^2, "cen")

    # Update the order from hclust to match the position in the original data frame
    curr_order <- loc_gapping[i - 1] + hc$order

    # Add the local ordering to the global ordering
    ordering <- c(ordering, curr_order)
  }

  # Update the data and annotation data frames ordering
  num_data <- num_data[ordering, ]
  annotation_row <- annotation_row[ordering, ]
  row_names <- row_names[ordering]

  # If single column data frame this re-ordering converts to a vector
  if (!is.data.frame(annotation_row)) {
    annotation_row <- as.data.frame(annotation_row)
  }

  # Re-impose the row names to enable correct annotation
  row.names(num_data) <- row_names
  row.names(annotation_row) <- row_names

  # Create gaps for the columns of the heatmap
  if (use_col_gaps) {
    num_col_gaps <- ncol(num_data) - 1
    col_gaps <- 1:num_col_gaps
  } else {
    col_gaps <- 0
  }

  # Heatmap
  annotated_heatmap(num_data, annotation_row,
    main = main,
    cluster_row = FALSE,
    cluster_cols = FALSE,
    gaps_row = gaps,
    gaps_col = col_gaps,
    ...
  )
}

# === Wrapper ==================================================================
#' @title MCMC out
#' @description Returns mean, variance and similarity posteriors from Gibbs sampling with
#' option of pheatmap
#'
#' @param MS_object A dataset in the format used by pRolocdata.
#' @param class_labels_0 An optional prior for clusters in MS_object. If NULL
#' defaults to a randomly generated set using the proportion in the labelled
#' data.
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
#' @param train: instruction to include all data (NULL), labelled data (TRUE) or
#' unlabelled data (FALSE). Default is NULL.
#' @param num_iter The number of iterations to sample over.
#' @param burn The number of iterations to record after (i.e. the burn-in).
#' @param thinning The step between iterations for which results are recorded in
#' the mcmc output.
#' @param heat_plot A bool. Instructs saving and printing of heatmap of
#' similarity matrix from Gibbs sampling. Default is TRUE.
#' @param main String. The title for heatmap, default is "heatmap_for_similarity".
#' @param cluster_row A bool. Instructs pheatmap to cluster rows using a tree.
#' @param cluster_cols: A bool. instructs pheatmap to cluster columns using a
#' tree.
#' @param fontsize: The size of font in pheatmap.
#' @param fontsize_row: The fontsize in rows in pheatmap.
#' @param fontsize_col: The fontsize in columns in pheatmap.
#' @param entropy_plot A bool instructing function to save a plot of the entropy
#' across all iterations. Default is TRUE.
#' @param window_length A number. Input to entropy ploy function. Default is
#' min(25, num_iter / 5).
#' @param mean_tolerance Input to entropy ploy function. Default is 0.0005.
#' @param sd_tolerance Input to entropy ploy function. Default is 0.0005.
#' @param outlier A bool instructing the sampler to consider an additional
#' outlier cluster following a t-distribution
#' @param t_df The degrees of freedom for the outlier t-distribution (default
#' is 4)
#' @param prediction_threshold The minimum proportion of recorded iterations
#' for which a point is in its most common cluster for which a prediction is
#' returned. If below this predicted class is NA.
#' @param record_posteriors A bool instructing the mcmc function to record the
#' posterior distributions of the mean and variance for each cluster
#' (default is FALSE)
#' @param normalise Bool instructing normalisation of continuous data (default
#' is false)
#' @return A named list including at least the output from the gibbs sampler,
#' but can include two pheatmaps and a scatter plot of the entropy over
#' iterations.
#' @examples
#' data("hyperLOPIT2015") # MS object from pRolocData
#' mcmc_object <- mcmc_out(hyperLOPIT2015,
#'   num_iter = 10000,
#'   burn = 1000,
#'   thinning = 50,
#'   outlier = TRUE,
#'   heat_plot = TRUE,
#'   main = "Gene clustering by organelle",
#'   prediction_threshold = 0.5
#' )
#'
#' Generate some nonsense categorical data
#' cat_data <- matrix(nrow = 1371, ncol = 5)
#' cat_data[, 1] <- c(rep(1, 300), rep(0, 1071))
#' cat_data[, 2] <- c(rep(0, 300), rep(1, 300), rep(0, 771))
#' cat_data[, 3] <- c(rep(0, 600), rep(1, 300), rep(0, 471))
#' cat_data[, 4] <- c(rep(0, 800), rep(1, 300), rep(0, 271))
#' cat_data[, 5] <- c(rep(0, 1100), rep(1, 271))
#'
#' # Complete noise
#' cat_data <-  matrix(sample(0:1, 1371 * 10, replace=TRUE), nrow = 1371, ncol = 10)
#'
#' stuff <- mcmc_out(HEK293T2011,
#'                   cat_data = cat_data,
#'                   num_iter = 10,
#'                   burn = 1,
#'                   thinning = 1,
#'                   outlier = TRUE,
#'                   heat_plot = T,
#'                   main = "Gene clustering by organelle",
#'                   prediction_threshold = 0.4,
#'                   sense_check_map = F
#' )
mcmc_out <- function(MS_object,
                     labels_0_1 = NULL,
                     args_1 = NULL,
                     cluster_weight_0_1 = NULL,
                     data_2 = NULL,
                     type_1 = "Gaussian",
                     type_2 = "Categorical",
                     cluster_weight_0_2 = 1,
                     args_2 = NULL,
                     labels_0_2 = NULL,
                     # n_clust_1 = NULL,
                     n_clust_2 = 50,
                     fix_vec_1 = NULL,
                     fix_vec_2 = rep(0, nrow(data_2)),
                     a_0 = 1,
                     b_0 = 0.2,
                     outlier_1 = FALSE,
                     t_df_1 = 4.0,
                     normalise_1 = FALSE,
                     outlier_2 = FALSE,
                     t_df_2 = 4.0,
                     normalise_2 = FALSE,
                     record_posteriors = FALSE,
                     train = NULL,
                     num_iter = NULL,
                     burn = floor(num_iter / 10),
                     thinning = 25,
                     heat_plot = TRUE,
                     main = "heatmap_for_similarity",
                     cluster_row = T,
                     cluster_cols = T,
                     fontsize = 10,
                     fontsize_row = 6,
                     fontsize_col = 6,
                     entropy_plot = TRUE,
                     window_length = min(25, num_iter / 5),
                     mean_tolerance = 0.0005,
                     sd_tolerance = 0.0005,
                     sense_check_map = TRUE,
                     prediction_threshold = 0.5) {
  # MS data
  MS_data <- MS_dataset(MS_object, train = train)
  
  mydata <- MS_data$data
  nk <- MS_data$nk
  row_names <- MS_data$row_names
  fix_vec_1 <- MS_data$fix_vec
  
  class_labels <- data.frame(Class = mydata$markers)
  
  classes_present <- unique(fData(markerMSnSet(MS_object))[, "markers"])
  
  rownames(class_labels) <- rownames(mydata)
  
  # Numerical data of interest for clustering
  num_data <- mydata %>%
    dplyr::select(-markers)
  
  # Parameters
  n_clust_1 <- length(classes_present)
  N <- nrow(num_data)
  d <- ncol(num_data)
  
  # Key to transforming from int to class
  class_labels_key <- data.frame(Class = classes_present) # , Class_num = 1:k)
  class_labels_key %<>%
    arrange(Class) %>%
    dplyr::mutate(Class_key = as.numeric(Class))
  
  class_labels %<>%
    mutate(Class_ind = as.numeric(mydata$markers))
  
  labels_0_1 <- cluster_label_prior(labels_0_1, nk, train, MS_object, n_clust_1, N)
  
  
  if (is.null(num_iter)) {
    num_iter <- floor(min((d^2) * 1000 / sqrt(N), 10000))
  }
  
  if (is.null(burn)) {
    burn <- floor(num_iter / 10)
  }
  
  if (burn > num_iter) {
    stop("Burn in exceeds total iterations. None will be recorded.\nStopping.")
  }
  
  # Warning if thinning is too high compared to num_iter - burn
  thinning_warning(thinning, num_iter, burn)
  
  # Prior on mass parameter for cluster  weights
  if (is.null(cluster_weight_0_1)) {
    print(paste0(
      "Creating vector of ", n_clust_1, " repetitions of ", 1,
      " for mass parameter prior for dataset 1."
    ))
    cluster_weight_0_1 <- rep(1, (n_clust_1))
  } else if (length(cluster_weight_0_1) < (n_clust_1)) {
    print(paste0(
      "Creating vector of ", n_clust_1, " repetitions of ", cluster_weight_0_1,
      " for mass parameter prior for dataset 1."
    ))
    cluster_weight_0_1 <- rep(cluster_weight_0_1, n_clust_1)
  }
  
  # Convert to matrix format
  num_data_mat <- as.matrix(num_data)
  data_2_mat <- as.matrix(data_2)
  
  if (is.null(data_2)) {
    
    gibbs <- gibbs_sampling(num_data_mat, n_clust_1, labels_0_1, fix_vec_1,
                            d = d,
                            N = N,
                            num_iter = num_iter,
                            burn = burn,
                            mu_0 = args_1$mu_0,
                            df_0 = args_1$df_0,
                            scale_0 = args_1$scale_0,
                            lambda_0 = args_1$lambda_0,
                            concentration_0 = cluster_weight_0_1,
                            thinning = thinning,
                            outlier = outlier_1,
                            t_df = t_df_1,
                            record_posteriors = record_posteriors,
                            normalise = normalise_1
    )
  } else {

    gibbs <- mdi(num_data_mat, data_2_mat,
                 args_1 = args_1,
                 args_2 = args_2,
                 type_1 = type_1,
                 type_2 = type_2,
                 n_clust_1 = n_clust_1,
                 n_clust_2 = n_clust_2,
                 labels_0_1 = labels_0_1,
                 labels_0_2 = labels_0_2,
                 d = d,
                 N = N,
                 a_0 = 1,
                 b_0 = 0.2,
                 fix_vec_1 = fix_vec_1,
                 fix_vec_2 = fix_vec_2,
                 cluster_weight_0_1 = cluster_weight_0_1,
                 cluster_weight_0_2 = cluster_weight_0_2,
                 num_iter = num_iter,
                 burn = burn,
                 thinning = thinning,
                 outlier_1 = outlier_1,
                 t_df_1 = t_df_1,
                 normalise_1 = normalise_1,
                 outlier_2 = outlier_2,
                 t_df_2 = t_df_2,
                 normalise_2 = normalise_2,
                 record_posteriors = record_posteriors
    )
  }
  
  print("Gibbs sampling complete")
  
  if (is.null(data_2)) {
    class_record <- gibbs$class_record
  } else {
    class_record <- gibbs$class_record_1
  }
  
  # Create a dataframe for the predicted class
  class_allocation_table <- with(
    stack(data.frame(t(class_record))),
    table(ind, values)
  )
  
  eff_iter <- ceiling((num_iter - burn) / thinning)

  
  # Create a column Class_key containing an integer in 1:k representing the most
  # common class allocation, and a Count column with the proportion of times the
  # entry was allocated to said class
  predicted_classes <- data.frame(
    Class_key =
      as.numeric(colnames(class_allocation_table)
                 [apply(
                   class_allocation_table,
                   1,
                   which.max
                 )]),
    Count = apply(class_allocation_table, 1, max) / eff_iter
  )
  
  # Change the prediction to NA for any entry with a proportion below the input
  # threshold
  predicted_classes[predicted_classes$Count < prediction_threshold, ] <- NA
  
  predicted_classes$Class <- class_labels_key$Class[match(
    predicted_classes$Class_key,
    class_labels_key$Class_key
  )]
  
  gibbs$predicted_class <- predicted_classes
  
  # Example input for annotation_row in pheatmap
  annotation_row <- class_labels %>% dplyr::select(Class)
  annotation_row %<>%
    mutate(Predicted_class = predicted_classes$Class)
  
  rownames(num_data) <- row_names
  # print(rownames(mydata[1:10,]))
  
  rownames(annotation_row) <- rownames(num_data)
  
  col_pal <- RColorBrewer::brewer.pal(9, "Blues")
  
  if (sense_check_map) {
    # pauls_heatmap <- annotated_heatmap(num_data, annotation_row,
    #   sort_by_col = Predicted_class,
    #   train = train,
    #   main = "Paul's sense check heatmap",
    #   cluster_row = FALSE,
    #   cluster_cols = FALSE,
    #   color = col_pal,
    #   fontsize = fontsize,
    #   fontsize_row = fontsize_row,
    #   fontsize_col = fontsize_col
    # )

    pauls_heatmap <- pheatmap_cluster_by_col(num_data, annotation_row, Predicted_class,
                                             main = "Paul's sense check heatmap",
                                             color = col_pal,
                                             fontsize = fontsize,
                                             fontsize_row = fontsize_row,
                                             fontsize_col = fontsize_col
    )
    
  }
  
  # return(pauls_heatmap)
  
  all_data <- dplyr::bind_cols(num_data, dplyr::select(gibbs$predicted_class, Class))
  all_data$Fixed <- fix_vec_1
  all_data$Protein <- rownames(num_data)
  
  # rownames(all_data) <- rownames(num_data)
  
  if (heat_plot) {
    
    # dissimilarity matrix
    if (is.null(data_2)) {
      sim <- gibbs$similarity
    } else {
      sim <- gibbs$similarity_1
    }
    
    dissim <- 1 - sim
    
    # Require names to associate data in annotation columns with original data
    colnames(dissim) <- rownames(num_data)
    rownames(dissim) <- rownames(num_data)
    
    col_pal <- RColorBrewer::brewer.pal(9, "Blues")
    
    heat_map <- annotated_heatmap(dissim, annotation_row,
                                  train = train,
                                  main = main,
                                  cluster_row = cluster_row,
                                  cluster_cols = cluster_cols,
                                  color = col_pal,
                                  fontsize = fontsize,
                                  fontsize_row = fontsize_row,
                                  fontsize_col = fontsize_col
    )
  }
  if (entropy_plot) {
    entropy_data <- data.frame(Index = 1:(num_iter + 1), Entropy = gibbs$entropy)
    
    rec_burn <- entropy_window(gibbs$entropy,
                               window_length = window_length,
                               mean_tolerance = mean_tolerance,
                               sd_tolerance = sd_tolerance
    )
    
    # Check if instantly ok
    rec_burn <- ifelse(is.null(rec_burn), 1, rec_burn)
    
    entropy_scatter <- ggplot(data = entropy_data, mapping = aes(x = Index, y = Entropy)) +
      geom_point() +
      geom_vline(mapping = aes(xintercept = rec_burn, colour = "Reccomended"), lty = 2) +
      geom_vline(mapping = aes(xintercept = burn, colour = "Implemented"), lty = 4) +
      ggtitle("Entropy over iterations including recommended and implemented burn") +
      xlab("Iteration") + ylab("Entropy") +
      scale_color_manual(name = "Burn", values = c(
        Reccomended = "red",
        Implemented = "blue"
      ))
  }
  if (heat_plot & entropy_plot) {
    return(list(
      gibbs = gibbs,
      heat_map = heat_map,
      entropy_plot = entropy_scatter,
      rec_burn = rec_burn,
      data = all_data
    ))
  }
  if (heat_plot) {
    return(list(
      gibbs = gibbs,
      heatmap = heat_map,
      data = all_data
    ))
  }
  if (entropy_plot) {
    return(list(
      gibbs = gibbs,
      entropy_plot = entropy_scatter,
      rec_burn = rec_burn,
      data = all_data
    ))
  }
  
  return(list(gibbs = gibbs, data = all_data, data_2 = data_2))
}


















# mcmc_out <- function(MS_object,
#                      class_labels_0 = NULL,
#                      mu_0 = NULL,
#                      df_0 = NULL,
#                      scale_0 = NULL,
#                      lambda_0 = 0.01,
#                      concentration_0 = NULL,
#                      cat_data = NULL,
#                      cluster_weight_priors_categorical = 1,
#                      phi_0 = NULL,
#                      c_clusters_label_0 = NULL,
#                      num_clusters_cat = 100,
#                      train = NULL,
#                      num_iter = NULL,
#                      burn = floor(num_iter / 10),
#                      thinning = 25,
#                      heat_plot = TRUE,
#                      main = "heatmap_for_similarity",
#                      cluster_row = T,
#                      cluster_cols = T,
#                      fontsize = 10,
#                      fontsize_row = 6,
#                      fontsize_col = 6,
#                      entropy_plot = TRUE,
#                      window_length = min(25, num_iter / 5),
#                      mean_tolerance = 0.0005,
#                      sd_tolerance = 0.0005,
#                      sense_check_map = TRUE,
#                      outlier = FALSE,
#                      t_df = 4.0,
#                      prediction_threshold = 0.6,
#                      record_posteriors = FALSE,
#                      normalise = FALSE) {
#   # MS data
#   MS_data <- MS_dataset(MS_object, train = train)
# 
#   mydata <- MS_data$data
#   nk <- MS_data$nk
#   row_names <- MS_data$row_names
#   fix_vec <- MS_data$fix_vec
# 
#   class_labels <- data.frame(Class = mydata$markers)
# 
#   classes_present <- unique(fData(markerMSnSet(MS_object))[, "markers"])
# 
#   rownames(class_labels) <- rownames(mydata)
# 
#   # Numerical data of interest for clustering
#   num_data <- mydata %>%
#     dplyr::select(-markers)
# 
#   # Parameters
#   k <- length(classes_present)
#   N <- nrow(num_data)
#   d <- ncol(num_data)
# 
#   # Key to transforming from int to class
#   class_labels_key <- data.frame(Class = classes_present) # , Class_num = 1:k)
#   class_labels_key %<>%
#     arrange(Class) %>%
#     dplyr::mutate(Class_key = as.numeric(Class))
# 
#   class_labels %<>%
#     mutate(Class_ind = as.numeric(mydata$markers))
# 
#   # if (outlier) {
#   #   outlier_row <- data.frame(Class = c("Outlier"), Class_key = c(k + 1))
#   #   class_labels_key <- suppressWarnings(bind_rows(class_labels_key, outlier_row))
#   # }
#   #
#   # if(outlier){
#   # fixed_labels <- as.numeric(fData(markerMSnSet(MS_object))[, "markers"])
#   # unfixed_labels <- rep(k + outlier, N - length(fixed_labels))
#   # class_labels_0 <- c(fixed_labels, unfixed_labels)
#   # } else {
# 
#   class_labels_0 <- cluster_label_prior(class_labels_0, nk, train, MS_object, k, N)
#   # }
# 
# 
#   if (is.null(num_iter)) {
#     num_iter <- floor(min((d^2) * 1000 / sqrt(N), 10000))
#   }
# 
#   if (is.null(burn)) {
#     burn <- floor(num_iter / 10)
#   }
# 
#   if (burn > num_iter) {
#     stop("Burn in exceeds total iterations. None will be recorded.\nStopping.")
#   }
# 
#   # Warning if thinning is too high compared to num_iter - burn
#   thinning_warning(thinning, num_iter, burn)
# 
#   # Prior on mass parameter for cluster  weights
#   if (is.null(concentration_0)) {
#     # if (is.null(train) | isTRUE(train)) {
#     #   concentration_0 <- cluster_weight_prior(mydata$markers, outlier = outlier)
#     # } else {
#     #   concentration_0 <- rep(0.1, k)
#     # }
#     concentration_0 <- rep(1, (k))
#   } else if (length(concentration_0) < (k)) {
#     print(paste0(
#       "Creating vector of ", k, " repetitions of ", concentration_0,
#       " for concentration prior."
#     ))
#     concentration_0 <- rep(concentration_0, k)
#   }
# 
#   gibbs <- gibbs_sampling(num_data, k, class_labels_0, fix_vec,
#     d = d,
#     N = N,
#     num_iter = num_iter,
#     burn = burn,
#     mu_0 = mu_0,
#     df_0 = df_0,
#     scale_0 = scale_0,
#     lambda_0 = lambda_0,
#     concentration_0 = concentration_0,
#     cat_data = cat_data,
#     cluster_weight_priors_categorical = cluster_weight_priors_categorical,
#     phi_0 = phi_0,
#     c_clusters_label_0 = c_clusters_label_0,
#     num_clusters_cat = num_clusters_cat,
#     thinning = thinning,
#     outlier = outlier,
#     t_df = t_df,
#     record_posteriors = record_posteriors,
#     normalise = normalise
#   )
# 
#   if (!is.null(cat_data)) {
#     gibbs <- mdi(num_data, cat_data,
#       args_1 = args_1,
#       args_2 = args_2,
#       type_1 = "Gaussian",
#       type_2 = "Categorical",
#       n_clust_1 = k,
#       n_clust_2 = num_clusters_cat,
#       labels_0_1 = class_labels_0,
#       labels_0_2 = c_clusters_label_0,
#       d = d,
#       N = N,
#       a_0 = 1,
#       b_0 = 0.2,
#       fix_vec_1 = fix_vec,
#       fix_vec_2 = rep(0, nrow(cat_data)),
#       cluster_weight_0_1 = concentration_0,
#       cluster_weight_0_2 = cluster_weight_priors_categorical,
#       num_iter = num_iter,
#       burn = burn,
#       thinning = thinning,
#       outlier_1 = outlier,
#       t_df_1 = t_df,
#       normalise_1 = normalise,
#       outlier_2 = FALSE,
#       t_df_2 = 4.0,
#       normalise_2 = FALSE,
#       record_posteriors = record_posteriors
#     )
#   }
# 
# 
# 
#   print("Gibbs sampling complete")
# 
#   if (is.null(cat_data)) {
#     class_record <- gibbs$class_record
#   } else {
#     class_record <- gibbs$class_record_1
#   }
# 
#   # Create a dataframe for the predicted class
#   class_allocation_table <- with(
#     stack(data.frame(t(class_record))),
#     table(ind, values)
#   )
# 
#   eff_iter <- ceiling((num_iter - burn) / thinning)
# 
#   # Create a column Class_key containing an integer in 1:k representing the most
#   # common class allocation, and a Count column with the proportion of times the
#   # entry was allocated to said class
#   predicted_classes <- data.frame(
#     Class_key =
#       as.numeric(colnames(class_allocation_table)
#       [apply(
#           class_allocation_table,
#           1,
#           which.max
#         )]),
#     Count = apply(class_allocation_table, 1, max) / eff_iter
#   )
# 
#   # Change the prediction to NA for any entry with a proportion below the input
#   # threshold
#   predicted_classes[predicted_classes$Count < prediction_threshold, ] <- NA
# 
#   predicted_classes$Class <- class_labels_key$Class[match(
#     predicted_classes$Class_key,
#     class_labels_key$Class_key
#   )]
# 
#   gibbs$predicted_class <- predicted_classes
# 
#   # Example input for annotation_row in pheatmap
#   annotation_row <- class_labels %>% dplyr::select(Class)
#   annotation_row %<>%
#     mutate(Predicted_class = predicted_classes$Class)
# 
#   rownames(num_data) <- row_names
#   # print(rownames(mydata[1:10,]))
# 
#   rownames(annotation_row) <- rownames(num_data)
# 
#   col_pal <- RColorBrewer::brewer.pal(9, "Blues")
# 
#   if (sense_check_map) {
#     pauls_heatmap <- annotated_heatmap(num_data, annotation_row,
#       sort_by_col = Predicted_class,
#       train = train,
#       main = "Paul's sense check heatmap",
#       cluster_row = FALSE,
#       cluster_cols = FALSE,
#       color = col_pal,
#       fontsize = fontsize,
#       fontsize_row = fontsize_row,
#       fontsize_col = fontsize_col
#     )
#   }
# 
#   # return(pauls_heatmap)
# 
#   all_data <- dplyr::bind_cols(num_data, dplyr::select(gibbs$predicted_class, Class))
#   all_data$Fixed <- fix_vec
#   all_data$Protein <- rownames(num_data)
# 
#   # rownames(all_data) <- rownames(num_data)
# 
#   if (heat_plot) {
# 
#     # dissimilarity matrix
#     if (is.null(cat_data)) {
#       sim <- gibbs$similarity
#     } else {
#       sim <- gibbs$similarity_1
#     }
# 
#     dissim <- 1 - sim
# 
#     # Require names to associate data in annotation columns with original data
#     colnames(dissim) <- rownames(num_data)
#     rownames(dissim) <- rownames(num_data)
# 
#     col_pal <- RColorBrewer::brewer.pal(9, "Blues")
# 
#     heat_map <- annotated_heatmap(dissim, annotation_row,
#       train = train,
#       main = main,
#       cluster_row = cluster_row,
#       cluster_cols = cluster_cols,
#       color = col_pal,
#       fontsize = fontsize,
#       fontsize_row = fontsize_row,
#       fontsize_col = fontsize_col
#     )
#   }
#   if (entropy_plot) {
#     entropy_data <- data.frame(Index = 1:num_iter, Entropy = gibbs$entropy)
# 
#     rec_burn <- entropy_window(gibbs$entropy,
#       window_length = window_length,
#       mean_tolerance = mean_tolerance,
#       sd_tolerance = sd_tolerance
#     )
# 
#     # Check if instantly ok
#     rec_burn <- ifelse(is.null(rec_burn), 1, rec_burn)
# 
#     entropy_scatter <- ggplot(data = entropy_data, mapping = aes(x = Index, y = Entropy)) +
#       geom_point() +
#       geom_vline(mapping = aes(xintercept = rec_burn, colour = "Reccomended"), lty = 2) +
#       geom_vline(mapping = aes(xintercept = burn, colour = "Implemented"), lty = 4) +
#       ggtitle("Entropy over iterations including recommended and implemented burn") +
#       xlab("Iteration") + ylab("Entropy") +
#       scale_color_manual(name = "Burn", values = c(
#         Reccomended = "red",
#         Implemented = "blue"
#       ))
#   }
#   if (heat_plot & entropy_plot) {
#     return(list(
#       gibbs = gibbs,
#       heat_map = heat_map,
#       entropy_plot = entropy_scatter,
#       rec_burn = rec_burn,
#       data = all_data
#     ))
#   }
#   if (heat_plot) {
#     return(list(
#       gibbs = gibbs,
#       heatmap = heat_map,
#       data = all_data
#     ))
#   }
#   if (entropy_plot) {
#     return(list(
#       gibbs = gibbs,
#       entropy_plot = entropy_scatter,
#       rec_burn = rec_burn,
#       data = all_data
#     ))
#   }
# 
#   return(list(gibbs = gibbs, data = all_data))
# }


# === Cross-Validation =========================================================
