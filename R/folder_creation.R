#!/usr/bin/env Rscript

# Script containing functions to create output folders

# === Folder creation ==========================================================

#' @title Create folder
#' @description Creates a folder in the appropriate location (default is the
#' current directory) with an option to overwrite.
#' @param folder_name A string. The name for the folder.
#' @param root_dir A string. The name of the root directory to create the folder
#' in. Default is the current working directory (i.e. ".").
#' @param overwrite A bool. Instructs overwriting if the desired directory is
#' found to already exist. Default is FALSE (in which case an error is thrown
#' if the directory does already exist).
create_folder <- function(folder_name,
                          root_dir = ".",
                          overwrite = FALSE,
                          load_results = FALSE) {
  # Create appropriate file path
  file_path <- file.path(root_dir, folder_name)

  if (dir.exists(file_path) & overwrite & !load_results) {
    unlink(file_path, T)
  }

  # If the directory does not already exist or we have instructed overwriting
  # create the directory. Else throw an error.
  if (!dir.exists(file_path) | (overwrite & !load_results)) {
    dir.create(file_path)
  } else if (dir.exists(file_path) & !load_results) {
    stop("Folder already exists - please check")
  }
}


#' @title Create dataset folders
#' @description Creates the directories to hold the dataset specific outputs
#' (i.e. parameter posteriors and allocation probabilities)
#'
#' @param n_clust Integer. The number of clusters in the dataset.
#' @param type String. The type of muixture model used, Gaussian or Categorical.
#' @param dataset_num Integer. The number identifying this dataset.
#' @param overwrite Bool. Instructs overwriting previously existing directories.
create_dataset_folders <- function(n_clust, type, dataset_num,
                                   overwrite = FALSE) {

  # Create the dataset specific folder
  dataset_folder <- paste0("dataset_", dataset_num)

  # create_folder(dataset_folder, overwrite = overwrite)

  # Create vector covering each of the components
  cluster_numbers <- 1:n_clust

  # Create folder and subfolder names for parameters based on the type of data
  if (type == "Gaussian") {
    mean_folder <- paste0(dataset_folder, "/", "mean")
    # create_folder(mean_folder, overwrite = overwrite)

    var_folder <- paste0(dataset_folder, "/", "variance")
    # create_folder(dataset_folder, overwrite = overwrite)

    mean_subfolders <- paste0(
      dataset_folder,
      "/",
      "mean/mean",
      "_",
      cluster_numbers
    )

    var_subfolders <- paste0(
      dataset_folder,
      "/",
      "variance/var",
      "_",
      cluster_numbers
    )

    param_subfolders <- c(
      mean_folder,
      mean_subfolders,
      var_folder,
      var_subfolders
    )
  } else { # If Categorical
    class_probs_folder <- paste0(dataset_folder, "/", "class_probs")
    # create_folder(class_probs_folder, overwrite = overwrite)

    class_probs_subfolder <- paste0(
      dataset_folder,
      "/",
      "class_probs/comp",
      "_",
      cluster_numbers
    )
    param_subfolders <- c(class_probs_folder, class_probs_subfolder)
  }

  # Folder for dataset allocation
  allocation_folder <- paste0(
    dataset_folder,
    "/",
    "allocation"
  )

  # Collect all the relevant pathnames in one vector
  dataset_folders <- c(dataset_folder, allocation_folder, param_subfolders)
}

#' @title Output folders
#' @description Creates all the relevant directories to save outputs from the
#' sampling
#' @param n_clust_1 An integer. The number of clusters in the first dataset.
#' @param n_clust_2 An integer. The number of clusters in the second dataset.
#' @param save_results Bool instructing program to save results to file. Default
#' is FALSE.
#' @param overwrite Bool instructing program to overwrite pre-existing results
#' if they exist in the current directroy. Default is FALSE with warnings anyway.
#' @param type_1 String indicating the type of mixture used in dataset 1. This
#' determines the folder for the appropriate parameters (i.e. if Gaussian
#' instructs creation of folders for a mean posterior and a variance posterior,
#' and a class probability posterior for "Categorical"). Default is "Gaussian".
#' @param type_2 String indicating the type of mixture used in dataset 2.
#' Default is "Categorical".
output_folders <- function(n_clust_1, n_clust_2,
                           save_results = TRUE,
                           load_results = FALSE,
                           overwrite = FALSE,
                           type_1 = "Gaussian",
                           type_2 = "Categorical") {
  if (save_results) {

    # Dataset parameter folder names
    dataset_1_folders <- create_dataset_folders(n_clust_1, type_1, 1,
      overwrite = overwrite
    )
    dataset_2_folders <- create_dataset_folders(n_clust_2, type_2, 2,
      overwrite = overwrite
    )

    dataset_folders <- c(dataset_1_folders, dataset_2_folders)

    # Create the root folder
    output_folder <- "output"
    suppressWarnings(create_folder(output_folder,
      overwrite = overwrite,
      load_results = load_results
    ))

    # Create the pathname to all the subdirectories of interest
    subdirectories <- c(
      c(
        "gaussian_allocation",
        "categorical_allocation",
        "outlier_allocation"
      ),
      dataset_folders
    )

    # Cycle theough these creating the folders
    for (subdir in subdirectories) {
      folder_name <- paste0(output_folder, "/", subdir)
      create_folder(folder_name,
        overwrite = overwrite,
        load_results = load_results
      )
    }
  }
}
