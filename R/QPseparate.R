# Function to apply QPsolve across a population of animals whose names are contained
# in listBreed
# 
# @param listBreed vector of animal IDs for which to test using QPsolve
# @param Y numeric matrix of genotypes from all animals in population
# @param X numeric matrix of allele frequencies from reference animals
# @return data.frame of breed composition estimates
QPseparate <- function(listBreed, Y, X) {
  
  # For each column of Y that matches a name in listBreed, perform QPsolve
  result <- t(apply(Y[, colnames(Y) %in% listBreed], 2, QPsolve, X))
  
  return(result)
}