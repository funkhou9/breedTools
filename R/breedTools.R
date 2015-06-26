#' breedTools: Predict ancestry of unknown animals
#' 
#' The breedTools package contains functions used to determine the breed/ancestry composition
#' for test animals from their genotyping data alone. It does so using genotyping data from 
#' a set of reference populations.
#' \cr\cr
#' There are a few ways in which breedTools can infer ancestry. These are:
#' \itemize{
#'    \item Predict overall breed/ancestry proportions genome wide.
#'    \item Predict overall breed/ancestry proportions genome wide
#'          using additional constraints when one of the parents of the test
#'          animal is genotyped.
#'    \item Estimate breed posterior probabilities based only on looking at
#'          limited regions of the genome.
#' }
#' @docType package
#' @name breedTools
NULL
