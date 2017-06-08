#' An example data set
#'
#' This simulated data list is for demonstration.
#' @docType data
#' @name sim_dat
#' @return A list containing
#' \item{Ym}{A N by p outcome data from N clusters/batches/experiments; and p is the number of samples within each cluster. 
#' The first sample within each cluster is a reference sample with a different error variance than other samples. Missing values are coded as NAs. Note the model allows unbalanced data.}
#' \item{X}{A covariate array of dimension of N by k by p, where k is the number of covariates.}
#' @examples data(sim_dat)
NULL
