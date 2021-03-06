#!/usr/bin/env Rscript

# Functions to prepare the data and various parameters for the sampling

# === PRE-PROCESSING ===========================================================

# --- Parameters  --------------------------------------------------------------

#' @title Check Thinning Factor
#' @description Generates an error or warning if the relative values of
#' thinning, num_iter and burn (the arguments controlling the number of recorded
#' MCMC samples) are strange.
#' @return None.
checkThinningFactor <- function(thinning, num_iter, burn) {
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

# Old name: empirical_bayes_gaussian
#' @title Generate Gaussian Empirical Bayes Prior
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
#'  @export
generateGaussianEmpiricalBayesPrior <- function(data,
                                                mu_0,
                                                df_0,
                                                scale_0,
                                                N,
                                                k,
                                                d,
                                                lambda_0 = 0.01) {
  parameters <- list()
  if (is.null(mu_0)) {
    mu_0 <- colMeans(data)
  }

  if (is.null(df_0)) {
    df_0 <- d + 2
  }

  if (is.null(scale_0)) {
    scale_0 <- diag(colSums((data - mean(data))^2) / N) / (k^(1 / d))
  }
  parameters$mu_0 <- mu_0
  parameters$df_0 <- df_0
  parameters$scale_0 <- scale_0
  parameters$lambda_0 <- lambda_0

  return(parameters)
}

# Old name: phi_prior
#' @title Phi prior
#' @description Generates a prior for the phi vector for each variable for the
#' Dirichlet distribution
#' @param matrix_data A matrix of data.
#' @return A list of vectors of the proportion of each level across all of
#' matrix_data.
#' @export
generateCategoricalEmpiricalBayesPrior <- function(matrix_data) {

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

# Old name: gaussian_arguments
#' @title Create Gaussian Arguments
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
#' args <- CreateGaussianArguments(data, 5)
createGaussianArguments <- function(data, n_clust,
                                    mu_0 = NULL,
                                    scale_0 = NULL,
                                    lambda_0 = NULL,
                                    df_0 = NULL) {
  N <- nrow(data)
  d <- ncol(data)

  if (is.null(lambda_0)) {
    lambda_0 <- 0.01
  }

  args <- generateGaussianEmpiricalBayesPrior(data,
    mu_0,
    df_0,
    scale_0,
    N,
    n_clust,
    d,
    lambda_0 = lambda_0
  )

  args
}

# Old name: categorical_arguments
#' @title Create Categorical Arguments
#' @description Creates a named list of the priors required for a mixture of
#' Dirichlets. The priors are generated using an empirical method (hence the
#' need for inputting the data and the number of clusters).
#' @param data A matrix of continuous values.
#' @param n_clust A unsigned integer; the number of clusters to be used in the
#' clustering method.
#' @return A named list ready to be used as input into the MDI function.
createCategoricalArguments <- function(data, n_clust) {
  args <- list(phi = generateCategoricalEmpiricalBayesPrior(data))
  args
}

#' @title Cluster weight prior
#' @description Produces a vector of prior weights for clusters
generateClusterWeightPrior <- function(train_data, outlier = FALSE) {
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

# declare_cluster_weights
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
createClusterWeights <- function(fix_vec, clust_labels, n_clust,
                                 weight_0 = NULL) {
  if (any(fix_vec == 1) & is.null(weight_0)) {
    relevant_labels <- clust_labels[fix_vec == 1]
    weight_0 <- unname(table(relevant_labels) / sum(table(relevant_labels)))
  } else {
    if (is.null(weight_0)) {
      weight_0 <- 1
    }
    if (length(weight_0) < n_clust) {
      # print(paste0(
      #   "Creating vector of ",
      #   n_clust,
      #   " repetitions of ",
      #   weight_0,
      #   " for mass parameter prior."
      # ))
      weight_0 <- rep(weight_0, n_clust)
    }
  }
  weight_0
}

# Old name: cluster_label_prior
#' @title Generate Cluster Label Prior
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
#' @importFrom MSnbase fData
#' @importFrom pRoloc markerMSnSet
generateClusterLabelPrior <- function(class_labels_0,
                                      train,
                                      MS_object,
                                      k,
                                      N) {
  # Generate class labels
  if (is.null(class_labels_0)) {
    if (is.null(train)) {
      fixed_labels <- as.numeric(fData(markerMSnSet(MS_object))[, "markers"])
      class_labels_0 <- c(fixed_labels, sample(1:k, N - length(fixed_labels),
        replace = T
      ))
    } else if (isTRUE(train)) {
      class_labels_0 <- as.numeric(MSnbase::fData(pRoloc::markerMSnSet(MS_object))[, "markers"])
    } else {
      class_labels_0 <- sample(1:k, N, replace = T)
    }
  }
  return(class_labels_0)
}
