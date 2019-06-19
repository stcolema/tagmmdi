#!/usr/bin/env Rscript

# Wrapper function for the package - this is the key function and carries out
# the appropriate pre-sampling, sampling and post-sampling steps

# === WRAPPER FUNCTION==========================================================
#' @title MCMC out
#' @description Returns mean, variance and similarity posteriors from Gibbs
#' sampling with option of pheatmap
#' @param MS_object A dataset in the format used by pRolocdata.
#' @param labels_0_1 An optional prior for clusters in MS_object. If NULL
#' defaults to a randomly generated set using the proportion in the labelled
#' data.
#' @param args_1 A named list containing the relevant objects for the data type
#' of MS_object. If NULL this is generated based on the input type_1 and using
#' the function gaussian_arguments or categorical arguments.
#' @param cluster_weight_0_1 The prior for dirichlet distribution of cluster
#' weights.
#' @param data_2 Second dataset for MDI. Either a matrix of continuous data,
#' categorical data or NULL (in which case MDI is not used).
#' Default is NULL in which case a generic mixture of Gaussians is used.
#' @param type_1 String denoting the type of data in MS_object. One of
#' "Gaussian" or "Categorical" (shorthand "G" and "C" also suffice). Default is
#' "Gaussian".
#' @param type_2 String denoting the type of data in MS_object. One of
#' "Gaussian" or "Categorical" (shorthand "G" and "C" also suffice). Default is
#' "Categorical".
#' @param cluster_weight_0_2 Vector of the prior on cluster
#' weights in the categorical data. If a single number is given a uniform vector
#' is generated.
#' @param args_2 A named list containing the relevant objects for the data type
#' of data_2. If NULL this is generated based on the input type_1 and using
#' the function gaussian_arguments or categorical arguments.
#' @param labels_0_2 Vector of labels for the prior clustering of the
#' categorical data.
#' @param n_clust_2 Integer of the number of clusters to have as a
#' maximum in the categorical dataset. Default is 100.
#' @param fix_vec_1 A vector of 0's and 1's indicating which points are fixed in
#' their initial allocation in dataset 1. Default of NULL means that the
#' information in MS_object is used to define this (any labelled points in the
#' dataset are held in this allocation).
#' @param fix_vec_2 A vector of 0's and 1's indicating which points are fixed in
#' their initial allocation for dataset 2. Default of all points unfixed.
#' @param a_0 Prior shape parameter for gamma distribution of strategic latent
#' variable, v, in MDI. Default is 1.
#' @param b_0 Prior rate parameter for gamma distribution of strategic latent
#' variable, v, in MDI. Default is 0.2.
#' @param outlier_1 A bool instructing the sampler to consider an outlier
#' cluster following a t-distribution (i.e. implementing TAGM). Only
#' applicable to gaussian datasets.
#' @param t_df_1 The degrees of freedom for the outlier t-distribution (default
#' is 4).
#' @param normalise_1 Bool instructing normalisation of continuous data in
#' dataset 1 (default is false).
#' @param outlier_2 A bool instructing the sampler to consider an outlier
#' cluster following a t-distribution (i.e. implementing TAGM). Only
#' applicable to gaussian datasets.
#' @param t_df_2 The degrees of freedom for the outlier t-distribution (default
#' is 4).
#' @param normalise_2 Bool instructing normalisation of continuous data in
#' dataset 2 (default is false).
#' @param record_posteriors A bool instructing the mcmc function to record the
#' posterior distributions of the mean and variance for each cluster
#' (default is TRUE)
#' @param save_results Bool instructing program to save results to file. Default
#' is FALSE.
#' @param overwrite Bool instructing program to overwrite pre-existing results
#' if they exist in the current directroy. Default is FALSE with warnings anyway.
#' @param train instruction to include all data (NULL), labelled data (TRUE) or
#' unlabelled data (FALSE). Default is NULL.
#' @param num_iter The number of iterations to sample over.
#' @param burn The number of iterations to record after (i.e. the burn-in).
#' @param thinning The step between iterations for which results are recorded in
#' the mcmc output.
#' @param heat_plot A bool. Instructs saving and printing of heatmap of
#' similarity matrix from Gibbs sampling. Default is TRUE.
#' @param main String. The title for heatmap, default is "heatmap_for_similarity".
#' @param cluster_row A bool. Instructs pheatmap to cluster rows using a tree.
#' @param cluster_cols A bool. instructs pheatmap to cluster columns using a
#' tree.
#' @param fontsize The size of font in pheatmap.
#' @param fontsize_row The fontsize in rows in pheatmap.
#' @param fontsize_col The fontsize in columns in pheatmap.
#' @param entropy_plot A bool instructing function to save a plot of the entropy
#' across all iterations. Default is TRUE.
#' @param window_length A number. Input to entropy ploy function. Default is
#' min(25, num_iter / 5).
#' @param mean_tolerance Input to entropy ploy function. Default is 0.0005.
#' @param sd_tolerance Input to entropy ploy function. Default is 0.0005.
#' @param sense_check_map A bool instructing function to save a heatmap of the
#' componenet level clusterings.
#' @param sense_check_main String. Title for sense_check_map. Default is
#' "component_level_clustering".
#' @param prediction_threshold The minimum proportion of recorded iterations
#' for which a point is in its most common cluster for which a prediction is
#' returned. If below this predicted class is NA.
#' @return A named list including at least the output from the gibbs sampler,
#' but can include two pheatmaps and a scatter plot of the entropy over
#' iterations.
#' @examples
#' Single dataset clustering
#' data("hyperLOPIT2015") # MS object from pRolocData
#' mcmc_object <- mcmc_out(hyperLOPIT2015,
#'   num_iter = 10000,
#'   burn = 1000,
#'   thinning = 50,
#'   outlier_1 = TRUE,
#'   heat_plot = TRUE,
#'   main = "Gene clustering by organelle",
#'   prediction_threshold = 0.5
#' )
#'
#' Use some categorical data
#' cat_data <- as.matrix(exprs(tan2009r1goCC))
#'
#' Implement MDI
#' stuff <- mcmc_out(tan2009r1,
#'                   data_2 = cat_data,
#'                   num_iter = 10,
#'                   burn = 1,
#'                   thinning = 1,
#'                   outlier_1 = TRUE,
#'                   heat_plot = T,
#'                   main = "Gene clustering by organelle",
#'                   prediction_threshold = 0.4,
#'                   sense_check_map = F
#' )
#' @importFrom dplyr select arrange mutate bind_cols
#' @importFrom ggplot2 ggplot aes geom_point geom_vline ggtitle xlab ylab scale_color_manual
#' @importFrom MSnbase fData
#' @importFrom pRoloc markerMSnSet
#' @importFrom RColorBrewer brewer.pal
#' @importFrom attempt try_catch
mcmc_out <- function(MS_object,
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
                     t_df_1 = 4.0,
                     normalise_1 = FALSE,
                     outlier_2 = FALSE,
                     t_df_2 = 4.0,
                     normalise_2 = FALSE,
                     record_posteriors = TRUE,
                     save_results = FALSE,
                     load_results = FALSE,
                     overwrite = FALSE,
                     train = NULL,
                     num_iter = NULL,
                     burn = floor(num_iter / 10),
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
                     window_length = min(25, num_iter / 5),
                     mean_tolerance = 0.0005,
                     sd_tolerance = 0.0005,
                     sense_check_map = TRUE,
                     sense_check_map_2 = FALSE,
                     sense_check_main = "component_level_clustering",
                     prediction_threshold = 0.5,
                     verbose = FALSE) {

  # Consider allowing input of fcol in which case don't deselect markers, rather use:

  # Expression data
  num_data <- MS_object %>%
    Biobase::exprs()

  # Marker data
  marker_data <- getMarkers(MS_object, fcol = "markers", names = T, verbose = F)

  # Order by known and unknown
  row_names <- names(marker_data)[c(which(marker_data != "unknown"), which(marker_data == "unknown"))]

  # Reorder data
  num_data <- num_data[match(row_names, row.names(num_data)), ]
  marker_data <- marker_data[match(row_names, names(marker_data))]

  # Generate a vector to indicate known and unknown labels if not given
  if (is.null(fix_vec_1)) {
    fix_vec_1 <- (marker_data != "unknown") * 1
  }


  if (!is.null(data_2)) {
    # Unpack pRoloc dataset

    # Expression data
    data_2 <- exprs(data_2)

    # Common order with original data
    data_2 <- data_2[match(row_names, row.names(data_2)), ]

    # I think this can never happen
    row_names_2 <- row.names(data_2)

    if (sum(row_names != row_names_2) > 0) {
      print(sum(row_names != row_names_2))
      stop("Row names in datasets are not the same. Check compatible data.")
    }
  }

  # If a fix_vec is not given for dataset 2, allow all points to move
  # (i.e. unsupervised clustering)
  if (!is.null(data_2) & is.null(fix_vec_2)) {
    fix_vec_2 <- rep(0, nrow(data_2))
  }

  # Put marker data in a data.frame
  class_labels <- marker_data %>%
    as.data.frame() %>%
    magrittr::set_colnames("Class")

  # Find the organelles present in the current data
  classes_present <- MS_object %>% getMarkerClasses()

  # As this is called a few times for match()
  sorted_classes <- sort(classes_present)

  if (sum(row.names(class_labels) != row.names(num_data))) {
    stop("Row names in disagreement")
  }

  # Parameters for clustering
  n_clust_1 <- length(classes_present)
  N <- nrow(num_data)
  d <- ncol(num_data)

  # Folders to save to
  output_folders(n_clust_1, n_clust_2,
    save_results = save_results,
    load_results = load_results,
    overwrite = overwrite
  )

  # Try and load previous runs of the same implementation
  num_load <- 0

  if (load_results) {
    num_load <- attempt::try_catch(
      expr = length(list.files("./output/dataset_1/allocation")),
      .e = NULL,
      .w = NULL
    )
  }

  if (is.null(num_load)) {
    num_load <- 0
  }


  # Key to transforming from int to class
  class_labels_key <- data.frame(Class = sorted_classes) %>%
    magrittr::inset2(., "Class_key", value = as.numeric(.$Class))

  # Find the number to represent the classes based on the above class_labels_key
  class_numerical <- class_labels_key$Class_key[match(class_labels$Class, class_labels_key$Class)]

  class_labels <- class_labels %>%
    magrittr::inset2("Class_ind", value = class_numerical)

  # Generate initial labels to begin with in clustering
  labels_0_1 <- cluster_label_prior(
    labels_0_1,
    train,
    MS_object,
    n_clust_1,
    N
  )

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
    cluster_weight_0_1 <- rep(1, (n_clust_1))
  } else if (length(cluster_weight_0_1) < (n_clust_1)) {
    cluster_weight_0_1 <- rep(cluster_weight_0_1, n_clust_1)
  }



  # Convert to matrix format
  num_data_mat <- as.matrix(num_data)


  if (!is.null(data_2)) {
    data_2_mat <- as.matrix(data_2)
    if (!all.equal(row.names(num_data_mat), row.names(data_2_mat))) {
      print(sum(row.names(num_data_mat) != row.names(data_2_mat)))
      stop("Row names are not the same")
    }
  }

  if(verbose){
    cat("Beginning clustering.\n")
  }
  
  # If no second dataset, run a single mixture model
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
  }

  # If second dataset is present, implement MDI
  else {
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
      record_posteriors = record_posteriors,
      save_results = save_results,
      load_results = load_results,
      num_load = num_load
    )
  }
  
  if(verbose){
    cat("Clustering finished.\n")
  }
  
  # print("Gibbs sampling complete")

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

  # The number of iterations for which results are recorded (+1 from start at 0)
  eff_iter <- ceiling((num_iter + 1 - burn) / thinning)

  # Create a column Class_key containing an integer in 1:k representing the most
  # common class allocation, and a Count column with the proportion of times the
  # entry was allocated to said class
  predicted_classes <- data.frame(
    Class_key =
      as.numeric(colnames(class_allocation_table)
      [apply(
          # gibbs$allocation_mat_1,
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

  gibbs$predicted_class <- predicted_classes %>%
    magrittr::set_rownames(row.names(class_labels))


  # Check if known labels are the same in predicted version - if not somethign is broken
  if (sum(as.character(class_labels$Class[which(as.logical(fix_vec_1))]) != as.character(predicted_classes$Class[which(as.logical(fix_vec_1))]))) {
    print(sum(as.character(class_labels$Class[which(as.logical(fix_vec_1))]) != as.character(predicted_classes$Class[which(as.logical(fix_vec_1))])))
    stop("Classes not matching for known")
  }

  # Example input for annotation_row in pheatmap
  annotation_row <- data.frame(
    Class = class_labels$Class,
    Predicted_class = predicted_classes$Class
  ) %>%
    magrittr::set_rownames(row.names(class_labels))

  # Use NAs for blank space on annotation row (rather than an additional level)
  annotation_row$Class[annotation_row$Class == "unknown"] <- NA

  # Colour palette for heatmap (white for 0, blue for 1)
  col_pal <- grDevices::colorRampPalette(c("white", "#146EB4"))(100)

  # Define the breaks in the heatmap (i.e. that 0 is white and 1 is blue)
  palette_length <- length(col_pal)
  my_breaks <- c(
    seq(1 / palette_length, 1, length.out = palette_length)
  )

  # col_pal <- RColorBrewer::brewer.pal(9, "Blues")

  if (sense_check_map) {
    if(verbose){
      cat("Constructing heat plot of expression data annotated by clustering.\n")
    }
    
    component_heat_map <- pheatmap_cluster_by_col(num_data,
      annotation_row,
      Predicted_class,
      main = sense_check_main,
      color = col_pal,
      fontsize = fontsize,
      fontsize_row = fontsize_row,
      fontsize_col = fontsize_col,
      breaks = my_breaks
    )
    if (!is.null(data_2) & sense_check_map_2) {
      component_heat_map <- pheatmap_cluster_by_col(data_2,
        annotation_row,
        Predicted_class,
        main = sense_check_main,
        color = col_pal,
        fontsize = fontsize,
        fontsize_row = fontsize_row,
        fontsize_col = fontsize_col,
        breaks = my_breaks
      )

      # annotated_heatmap(data_2, annotation_row,
      #                   train = train,
      #                   main = main,
      #                   cluster_row = cluster_row,
      #                   cluster_cols = cluster_cols,
      #                   color = col_pal,
      #                   fontsize = fontsize,
      #                   fontsize_row = fontsize_row,
      #                   fontsize_col = fontsize_col,
      #                   breaks = my_breaks)
      #
      #   pheatmap_cluster_by_col(data_2,
      #                                               annotation_row,
      #                                               Predicted_class,
      #                                               main = sense_check_main,
      #                                               color = col_pal,
      #                                               fontsize = fontsize,
      #                                               fontsize_row = fontsize_row,
      #                                               fontsize_col = fontsize_col,
      #                                               breaks = my_breaks
      # )
    }
  }

  # Create a data.frame for the output
  all_data <- cbind(
    num_data,
    gibbs$predicted_class$Class,
    fix_vec_1,
    rownames(num_data),
    class_labels$Class
  ) %>% as.data.frame()
  
  # print(head(all_data))
  # print(head(all_data))
  
  # Set the column names
  colnames(all_data)[(d + 1):(d + 4)] <- c(
    "Prediction",
    "Fixed",
    "Protein",
    "Organelle"
  )

  if (heat_plot) {
    if(verbose){
      cat("Constructing heat plot of posterior similarity matrix.\n")
      # print("")
    }

    # dissimilarity matrix
    if (is.null(data_2)) {
      sim <- gibbs$similarity
    } else {
      sim <- gibbs$similarity_1
    }

    # dissim <- 1 - sim

    # Require names to associate data in annotation columns with original data
    colnames(sim) <- rownames(num_data)
    rownames(sim) <- rownames(num_data)

    # col_pal <- RColorBrewer::brewer.pal(9, "Blues")

    heat_map <- annotated_heatmap(sim, annotation_row,
      train = train,
      main = main,
      cluster_row = cluster_row,
      cluster_cols = cluster_cols,
      color = col_pal,
      fontsize = fontsize,
      fontsize_row = fontsize_row,
      fontsize_col = fontsize_col,
      breaks = my_breaks
    )
    if (!is.null(data_2) & heat_plot_2) {
      sim_2 <- gibbs$similarity_2

      dissim_2 <- 1 - sim_2

      # Require names to associate data in annotation columns with original data
      colnames(dissim_2) <- rownames(num_data)
      rownames(dissim_2) <- rownames(num_data)

      # col_pal <- RColorBrewer::brewer.pal(9, "Blues")

      heat_map_2 <- annotated_heatmap(dissim_2, annotation_row,
        train = train,
        main = main,
        cluster_row = cluster_row,
        cluster_cols = cluster_cols,
        color = col_pal,
        fontsize = fontsize,
        fontsize_row = fontsize_row,
        fontsize_col = fontsize_col,
        breaks = my_breaks
      )
    }
  }

  # print("Entropy")

  if (entropy_plot) {
    if(verbose){
      cat("Constructing entropy plot.\n")
    }
    entropy_data <- data.frame(
      Index = 1:(num_iter + 1),
      Entropy = gibbs$entropy
    )

    rec_burn <- entropy_window(gibbs$entropy,
      window_length = window_length,
      mean_tolerance = mean_tolerance,
      sd_tolerance = sd_tolerance
    )

    # Check if instantly ok
    rec_burn <- ifelse(is.null(rec_burn), 1, rec_burn)

    entropy_scatter <- ggplot2::ggplot(
      data = entropy_data,
      mapping = ggplot2::aes(x = Index, y = Entropy)
    ) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(
        mapping = ggplot2::aes(
          xintercept = rec_burn,
          colour = "Reccomended"
        ),
        lty = 2
      ) +
      ggplot2::geom_vline(
        mapping = ggplot2::aes(
          xintercept = burn,
          colour = "Implemented"
        ),
        lty = 4
      ) +
      ggplot2::ggtitle("Entropy over iterations including recommended and implemented burn") +
      ggplot2::xlab("Iteration") + ggplot2::ylab("Entropy") +
      ggplot2::scale_color_manual(name = "Burn", values = c(
        Reccomended = "red",
        Implemented = "blue"
      ))
  }

  # print("Sense")

  if (sense_check_map & heat_plot & entropy_plot) {
    return(list(
      gibbs = gibbs,
      cluster_map = component_heat_map,
      heat_map = heat_map,
      entropy_plot = entropy_scatter,
      rec_burn = rec_burn,
      data = all_data
    ))
  }
  if (sense_check_map & heat_plot) {
    return(list(
      gibbs = gibbs,
      cluster_map = component_heat_map,
      heat_map = heat_map,
      data = all_data
    ))
  }
  if (sense_check_map & entropy_plot) {
    return(list(
      gibbs = gibbs,
      cluster_map = component_heat_map,
      entropy_plot = entropy_scatter,
      rec_burn = rec_burn,
      data = all_data
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
  if (sense_check_map) {
    return(list(
      gibbs = gibbs,
      cluster_map = component_heat_map,
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

  return(list(gibbs = gibbs, data = all_data))
}
