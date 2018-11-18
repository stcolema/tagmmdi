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
                          overwrite = FALSE){
  # Create appropriate file path
  file_path <- file.path(root_dir, folder_name)
  
  # If the directory does not already exist or we have instructed overwriting
  # create the directory. Else throw an error.
  if (!dir.exists(file_path) | overwrite) {
    dir.create(file_path)
  } else {
    stop("Folder already exists - please check")
  }
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
                           overwrite = FALSE,
                           type_1 = "Gaussian",
                           type_2 = "Categorical"){
  if(save_results){
    
    # Create vector covering each of the components
    cluster_numbers_1 <- 1:n_clust_1
    cluster_numbers_2 <- 1:n_clust_2
    
    # Create folder and subfolder names for parameters
    if(type_1 == "Gaussian"){
      mean_subfolders <- paste0("mean/mean",
                                "_", 
                                cluster_numbers_1)
      
      var_subfolders <- paste0("variance/var",
                               "_", 
                               cluster_numbers_1)
      
      type_1_folders <- c(mean_subfolders, var_subfolders)
    }
    
    if(type_2 == "Categorical"){
      type_2_folders <- paste0("class_probs/comp",
                               "_", 
                               cluster_numbers_2)
    }
    
    
    param_folders <- c(type_1_folders, type_2_folders)
    
    # Create the root folder
    output_folder <- "output"
    suppressWarnings(create_folder(output_folder, overwrite))
    
    # Create the pathname to all the subdirectories of interest
    subdirectories <- c(c("gaussian_allocation",
                          "categorical_allocation", 
                          "outlier_allocation",
                          "class_probs",
                          "allocation_1",
                          "allocation_2",
                          "mean",
                          "variance"),
                        param_folders)
    
    # Cycle theough these creating the folders
    for (subdir in subdirectories) {
      folder_name <- paste0(output_folder, "/", subdir)
      create_folder(folder_name, overwrite)
    }
  }
}
