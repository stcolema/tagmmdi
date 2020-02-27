#!/usr/bin/env Rscript

# script containing functions specifically related to cross-validation

# === Cross-Validation =========================================================

# Old name: mdi_cross_validate
#' @title Cross validate MDI
#' @description Returns the allocation confusion matrix and quadratic loss score
#' of MDI applied to the expression data from MS object and MS cat object.
#'
#' @param MS_object A dataset in the format used by pRolocdata.
#' @param MS_cat_object A second dataset in the format used by pRolocdata.
#' @param times An integer; the number of folds to run.
#' @param test_size A fraction between 0 and 1; the proportion of the known
#' points to use as a test set in each fold.
#' @param num_iter An integer; the number of iterations to use in each MDI
#' sampling.
#' @param burn The number of iterations to record after (i.e. the burn-in).
#' @param thinning The step between iterations for which results are recorded in
#' the mcmc output.
#' @param ... Any other inputs applicable to mcmc_out.
#' @return A named list of the confusion matrix for allocaitons across each fold
#' and the associated quadratic loss.
#' @examples
#' cv_metrics <- CrossValidateMDI(tan2009r1, tan2009r1goCC)
#' @importFrom BiocGenerics combine
#' @importFrom caret confusionMatrix
#' @importFrom MSnbase fData MSnSet exprs pData
#' @importFrom pRoloc markerMSnSet getMarkerClasses
#' @importFrom sampling strata
#' @export
crossValidateMDI <- function(MS_object,
                             MS_cat_object = NULL,
                             times = 10,
                             test_size = 0.2,
                             num_iter = 1000,
                             burn = floor(num_iter / 10),
                             thinning = 25,
                             n_clust_cat = 50,
                             ...) {
  
  marker.data <- pRoloc::markerMSnSet(MS_object)

  if (!is.null(MS_cat_object)) {
    marker.data.cat <- pRoloc::markerMSnSet(MS_cat_object)
  }

  X <- pRoloc:::subsetAsDataFrame(marker.data, "markers", train = TRUE)

  K <- length(pRoloc::getMarkerClasses(MS_object))

  .testPartitions <- .cmMatrices <- vector("list", length = times)

  f1score <- matrix(0, times)
  .f1Matrices <- matrix(0, times)
  cmlist <- vector("list", length = times)
  quadloss <- vector("list", length = times)


  # Create a data frame of the classes present and their associated number
  classes_pres <- pRoloc::getMarkerClasses(MS_object)
  class_key <- data.frame(Class = classes_pres, Key = 1:length(classes_pres))

  for (i in seq_along(cmlist)) {
    # get sizes
    .size <- ceiling(table(MSnbase::fData(marker.data)$markers) * test_size)

    # strata needs size to be ordered as they appear in data
    .size <- .size[unique(MSnbase::fData(marker.data)$markers)]

    # get strata indices
    test.idx <- sampling::strata(X, "markers",
      size = .size,
      method = "srswor"
    )$ID_unit


    ## 'unseen' test set
    .test <- MSnbase::MSnSet(
      MSnbase::exprs(marker.data)[test.idx, ],
      MSnbase::fData(marker.data[test.idx, ]),
      MSnbase::pData(marker.data)
    )

    ## 'seen' training set
    .train <- MSnbase::MSnSet(
      MSnbase::exprs(marker.data)[-test.idx, ],
      MSnbase::fData(marker.data[-test.idx, ]),
      MSnbase::pData(marker.data)
    )
    
    # save true marker labels
    test.markers <- MSnbase::fData(.test)$markers

    # create new combined MSnset
    mydata <- BiocGenerics::combine(.train, .test)

    # print("Combined")
    
    # Set levels of markers cateogries
    levels(MSnbase::fData(mydata)$markers) <- c(
      levels(
        MSnbase::fData(mydata)$markers
      ),
      "unknown"
    )

    # hide marker labels
    MSnbase::fData(mydata)[rownames(.test), "markers"] <- "unknown"
    
  
    if (!is.null(MS_cat_object)) {
      # create new combined MSnset
      cat_data <- BiocGenerics::combine(marker.data.cat[-test.idx, ],
        marker.data.cat[test.idx, ]
      )
      
      # Set levels of markers cateogries
      levels(MSnbase::fData(cat_data)$markers) <- c(
        levels(
          MSnbase::fData(cat_data)$markers
        ),
        "unknown"
      )
      
      # hide marker labels
      MSnbase::fData(cat_data)[rownames(.test), "markers"] <- "unknown"
      
    } else {
      cat_data <- NULL
    }

    # Fix training points, allow test points to move component
    fix_vec_1 <- c(
      rep(1, nrow(MSnbase::exprs(.train))),
      rep(0, nrow(MSnbase::exprs(.test)))
    )
   
    # I think this is a better way of doing the fix-vec
    fix_vec_test <- (MSnbase::fData(mydata)[, "markers"] != "unknown") * 1
    
    if(sum(fix_vec_test != fix_vec_1) > 0){
      stop("Fix vec odd")
    }
    
    # MDI
    params <- clusteringWrapper(mydata,
      data_2 = cat_data,
      fix_vec_1 = fix_vec_1,
      n_clust_2 = n_clust_cat,
      num_iter = num_iter,
      burn = burn,
      thinning = thinning,
      heat_plot = F,
      sense_check_map = F,
      entropy_plot = F,
      outlier_1 = T,
      prediction_threshold = 0.0001,
      record_posteriors = FALSE,
      save_results = FALSE,
      load_results = FALSE,
      overwrite = FALSE,
      ...
    )

    # Allocation for test data
    
    # Find the relevant indices from the prediction and output of the sampler
    indices_for_prediction <- match(
      rownames(MSnbase::fData(.test)),
      params$data$Protein
    )

    # MAP prediction
    map_pred <- params$data$Prediction[indices_for_prediction]

    # MCMC prediction
    mcmc_pred <- apply(params$gibbs$allocation_mat_1, 1, max)
    prediction_vec <- params$gibbs$allocation_mat_1 == mcmc_pred

    # print("define allocation matrix")
    
    # Allocation matrix
    test_alloc <- matrix(
      nrow = nrow(prediction_vec),
      as.numeric(prediction_vec)
    )[indices_for_prediction, ]

    # Find predictions - order not as new dataset as wrapper function orders
    # based on fixing points
    mcmc_pred <- apply(
      params$gibbs$allocation_mat_1[indices_for_prediction, ],
      1,
      which.max
    )

    mcmc_pred_mat <- params$gibbs$allocation_mat_1[indices_for_prediction, ]

    # Make predictions on test data
    mcmc_predictions <- class_key$Class[match(
      mcmc_pred,
      class_key$Key
    )]

    map_predictions <- class_key$Class[match(
      map_pred,
      class_key$Key
    )]

    # True allocation for test data
    reference <- factor(test.markers, levels = pRoloc::getMarkerClasses(mydata))

    # Confusion matrix for current fold
    cmlist[[i]] <- conf <- caret::confusionMatrix(
      data = mcmc_predictions,
      reference = reference
    )$table

    # Create allocation matrices for truth, filled initially with 0's
    allocmatrix <- matrix(0,
      nrow = length(test.idx),
      ncol = length(pRoloc::getMarkerClasses(mydata))
    )

    test_alloc_2 <- allocmatrix

    # The numbers associated with the classes (i.e. numerical representation of
    # the classes)
    class_numerics <- seq(1, length(unique(test.markers)))

    # create allocation matrix
    for (j in seq_along(test.idx)) {
      # The class the current individual belongs to
      alloc <- as.numeric(test.markers, class_numerics)[j]

      # Enter this in the allocation matrix
      allocmatrix[j, alloc] <- 1
      test_alloc_2[j, map_predictions[j]] <- 1
    }

    # Compute quadratic loss
    quadloss[[i]] <- sum((allocmatrix - mcmc_pred_mat)^2)
  }
  list(cmlist = cmlist, quadloss = quadloss)
}
