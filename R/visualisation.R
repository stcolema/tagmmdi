#!/usr/bin/env Rscript

# Functions related to visualisation of results - heatmaps of clusterings and
# PCA

# === VISUALISATION ============================================================

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
#' @importFrom attempt try_catch
#' @importFrom dplyr select arrange one_of bind_cols
#' @importFrom grDevices rainbow
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang enquo quo_text
annotated_heatmap <- function(input_data, annotation_row = NULL,
                              sort_by_col = NULL,
                              train = NULL,
                              ...) {
  if (is.null(annotation_row) & !(isTRUE(train) | is.null(train))) {
    stop("If data")
  }

  dissim <- input_data

  sort_by_col <- attempt::try_catch(
    expr = rlang::enquo(sort_by_col),
    .e = NULL,
    .w = NULL
  )

  # If the sort_by_col contains something other than NULL

  # has rlang been updated? this line no longer works
  # if(! is.null(!!sort_by_col)){
  if (rlang::quo_text(sort_by_col) != "NULL") {

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
        dplyr::arrange(!!sort_by_col)

      # Separate out the dataframes
      annotation_row <- sorted_data %>%
        dplyr::select(dplyr::one_of(names(annotation_row)))

      row_names <- combined_data %>%
        dplyr::select(Names)

      dissim <- sorted_data %>%
        dplyr::select(-dplyr::one_of(names(annotation_row))) %>%
        dplyr::select(-Names)

      # Redeclare row names (as dplyr strips these)
      row.names(annotation_row) <- row_names$Names
      row.names(dissim) <- row_names$Names
    } else {
      combined_data <- dplyr::bind_cols(dissim, row_names)

      # Arrange based on the user selected columm
      sorted_data <- combined_data %>%
        dplyr::arrange(!!sort_by_col)

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
    heat_map <- pheatmap::pheatmap(
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
#' @importFrom dplyr select
#' @importFrom rlang enquo
#' @importFrom stats hclust
#' @importFrom dplyr select
#' @importFrom rlang enquo
#' @importFrom stats hclust
pheatmap_cluster_by_col <- function(num_data, annotation_row, sort_col,
                                    main = "sense_check",
                                    use_col_gaps = TRUE,
                                    ...) {

  # save row names as dplyr removes them
  row_names <- row.names(num_data)

  # Enclose the column name using tidy evaluation
  sort_col <- rlang::enquo(sort_col)

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
    hc <- stats::hclust(dist(num_data[(loc_gapping[i - 1] + 1):loc_gapping[i], ])^2, "cen")

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

# --- PCA Plot -----------------------------------------------------------------
#' @title PCA MS object
#' @description Creates a PCA plot of a given MS object from pRolocdata grouping
#' based on a vector of integers.
#' @param MS_object A dataset from pRolocdata.
#' @param test_pred A vector of allocations.
#' @param ellipses A bool indicating if clusters should be contained within
#' appropriately coloured ellipses on the plot.
#' @param alpha.ind A number or vector of numbers between 0 and 1 inddicating
#' the transparency of points (1 indicates opaque, 0 transparent).
#' @return PCA plot
#' @importFrom FactoMineR PCA
#' @importFrom factoextra fviz_pca_ind
pca_ms_obj <- function(MS_object, test_pred, ellipses = FALSE, alpha.ind = 1) {
  ms_data <- MS_dataset(MS_object)
  rel_data <- ms_data$data
  last_col <- ncol(rel_data)
  data_pca <- rel_data[, -last_col]
  data_pca$class <- test_pred
  data_pca$class <- as.factor(data_pca$class)
  test_pca <- FactoMineR::PCA(data_pca[, -last_col], graph = FALSE)
  pca_plot <- factoextra::fviz_pca_ind(test_pca,
    geom.ind = "point", # show points only (nbut not "text")
    col.ind = data_pca$class, # color by groups
    addEllipses = ellipses, # Concentration ellipses
    legend.title = "Groups",
    alpha.ind = alpha.ind,
    rotate = T
  )
  pca_plot
}
