#!/usr/bin/env Rscript

# script containing functions specifically related to cross-validation

# === Cross-Validation =========================================================

#' @title MDI cross validate
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
#' cv_metrics <- mdi_cross_validate(tan2009r1, tan2009r1goCC)
#' @importFrom BiocGenerics combine
#' @importFrom caret confusionMatrix
#' @importFrom MSnbase fData MSnSet exprs pData
#' @importFrom pRoloc markerMSnSet getMarkerClasses
#' @importFrom sampling strata
mdi_cross_validate <- function(MS_object,
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

  # print("Here at least")
  
  f1score <- matrix(0, times)
  .f1Matrices <- matrix(0, times)
  cmlist <- vector("list", length = times)
  quadloss <- vector("list", length = times)

  # print("Vectors declared")
  
  # Create a data frame of the classes present and their associated number
  classes_pres <- pRoloc::getMarkerClasses(MS_object)
  class_key <- data.frame(Class = classes_pres, Key = 1:length(classes_pres))

  # print("Looping")
  
  for (i in seq_along(cmlist)) {
    # get sizes
    .size <- ceiling(table(MSnbase::fData(marker.data)$markers) * test_size)

    # strata needs size to be ordered as they appear in data
    .size <- .size[unique(MSnbase::fData(marker.data)$markers)]

    # print("Test indices")
    
    # get strata indices
    test.idx <- sampling::strata(X, "markers",
      size = .size,
      method = "srswor"
    )$ID_unit

    # print("Test data")
    
    ## 'unseen' test set
    .test1 <- MSnbase::MSnSet(
      MSnbase::exprs(marker.data)[test.idx, ],
      MSnbase::fData(marker.data[test.idx, ]),
      MSnbase::pData(marker.data)
    )

    # print("Training data")
    
    ## 'seen' training set
    .train1 <- MSnbase::MSnSet(
      MSnbase::exprs(marker.data)[-test.idx, ],
      MSnbase::fData(marker.data[-test.idx, ]),
      MSnbase::pData(marker.data)
    )
    
    # print("Test markers")

    # save true marker labels
    test.markers <- MSnbase::fData(.test1)$markers

    # print("Combine")
    
    # print(str(test.idx))
    # 
    # print(str(exprs(.test1)))
    # print(str(exprs(.train1)))
    
    # print(str(exprs(BiocGenerics::combine(.train1, .test1))))
    
    # print("OK")
    
    # print(str(rbind(exprs(.train1), exprs(.test1))))
    # 
    # print(str(MSnbase::fData(.train1)))
    # print(str(MSnbase::fData(.test1)))
    
    # create new combined MSnset
    mydata <- BiocGenerics::combine(.train1, .test1)

    # print("Combined")
    
    # Set levels of markers cateogries
    levels(MSnbase::fData(mydata)$markers) <- c(
      levels(
        MSnbase::fData(mydata)$markers
      ),
      "unknown"
    )
    
    # print("Let them be red")
    
    # hide marker labels
    MSnbase::fData(mydata)[rownames(.test1), "markers"] <- "unknown"
    
    
    # print("I'm alive")
    
    if (!is.null(MS_cat_object)) {
      # cat_data <- MSnbase::exprs(marker.data.cat)
      
      ## 'unseen' test set
      # .test2 <- MSnbase::MSnSet(
      #   MSnbase::exprs(marker.data.cat)[test.idx, ],
      #   MSnbase::fData(marker.data.cat[test.idx, ]),
      #   MSnbase::pData(marker.data.cat)
      # )
      # 
      # ## 'seen' training set
      # .train2 <- MSnbase::MSnSet(
      #   MSnbase::exprs(marker.data.cat)[-test.idx, ],
      #   MSnbase::fData(marker.data.cat[-test.idx, ]),
      #   MSnbase::pData(marker.data.cat)
      # )
      
      # create new combined MSnset
      # cat_data <- BiocGenerics::combine(.test2, .train2)
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
      comp_cat <- cat_data
      # MSnbase::fData(comp_cat)[rownames(.test2), "markers"] <- "unknown"
      MSnbase::fData(cat_data)[rownames(.test1), "markers"] <- "unknown"
      
      # if(sum(MSnbase::fData(comp_cat)[, "markers"] != MSnbase::fData(cat_data)[, "markers"])){
      #   print("Wrong inex used for cat test data")
      # }
      
    } else {
      cat_data <- NULL
    }

    # print("Are you?")
    
    # Fix training points, allow test points to move component
    fix_vec_1 <- c(
      rep(1, nrow(MSnbase::exprs(.train1))),
      rep(0, nrow(MSnbase::exprs(.test1)))
    )
    
    # print(MSnbase::fData(mydata)[, "markers"])
    
    fix_vec_test <- (MSnbase::fData(mydata)[, "markers"] != "unknown") * 1
    
    # print("Wagg")
    
    if(sum(fix_vec_test != fix_vec_1) > 0){
      stop("Fix vec odd")
    }

    # print("Cat data")
    # print(head(cat_data))
    # 
    # print("data")
    # print(head(mydata))
    
    # MDI
    params <- mcmc_out(mydata,
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

    # print("Out")
    
    # print(head(params$data))
    
    # Find the relevant indices from the prediction and output of the sampler
    indices_for_prediction <- match(
      rownames(MSnbase::fData(.test1)),
      params$data$Protein
    )

    # print("Further out")
    
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

    # print("SHOHE")
    
    # Find predictions - order not as new dataset as wrapper function orders
    # based on fixing points
    mcmc_pred <- apply(
      params$gibbs$allocation_mat_1[indices_for_prediction, ],
      1,
      which.max
    )

    # print("More $$$$")
    
    mcmc_pred_mat <- params$gibbs$allocation_mat_1[indices_for_prediction, ]

    # print("Past $$$$")
    
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
