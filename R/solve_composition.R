#' solve_composition
#' 
#' Computes genome-wide breed/ancestry composition using quadratic programming on a
#'  batch of animals.
#' 
#' @param Y numeric matrix of genotypes (rows) from all animals (columns) in population
#'  coded as dosage of allele B {0, 1, 2}
#' @param X numeric matrix of allele frequencies from reference animals. Frequencies are
#'  relative to allele B.
#' @param names character vector of names for each breed/ancestral population in the order
#'  that they appear in X.
#' @param ped data.frame giving pedigree information. Must be formatted "ID", "Sire", "Dam"
#' @param groups list of IDs catagoriezed by breed/population. If specified, output will be a list
#'  of results categorized by breed/population.
#' @return A data.frame or list of data.frames (if groups is !NULL) with breed/ancestry compostion
#'  results
#' @import quadprog
#' @export
solve_composition <- function(Y, X, names, ped = NULL, groups = NULL) {
  
  # At the moment can only have ped OR groups. They are not currently compatable
  #   although they should be
  
  # Functions require Y to be animals x snps. Transpose
  Y <- t(Y)
  
  # If ped is supplied, use QPsolve_par to compute genomic composition using
  #   only animals who have genotyped parents (by incorporating Sire genotype).
  if (!is.null(ped)) {
    mat_results <- lapply(colnames(Y), QPsolve_par, Y, X, ped)
    mat_results_tab <- do.call(rbind, mat_results)
    
    return (mat_results_tab)
    
  # Else if groups supplied - perform regular genomic computation
  #   and list results by groups
  } else if (!is.null(groups)) {
    
    # When using regular genomic computation - dosage must first be halved
    Y <- Y/2
    
    grouped_results <- lapply(groups, QPseparate, Y, X)
    return (grouped_results)
  
  # If neither using the ped or grouping option - just perform normal, unsegregated
  #   calculation
  } else {
    
    # When using regular genomic computation - dosage must first be halved
    Y <- Y/2
    
    results <- t(apply(Y, 2, QPsolve, X))
    return (results)
  }
}




