#!/usr/bin/env Rscript

# Function for judging convergence, specifically the entropy plot

# === Convergence ==============================================================

# Old name: entropy_window
#' @title Plot Entropy Diagnostic
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
PlotEntropyDiagnostic <- function(entropy_vec,
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
    win_2 <- entropy_vec[(i + window_length):min((i + 2 * window_length - 1), n)]

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
