#' Solve generic breed probabilities
#' 
#' Given a breed probability matrix (provided by \code{build_KBP()} or \code{build_MBP()}), and a set of test animals'
#' genotypes, use this function to derive breed probabilities for each test animal genotype based on
#' a set of chosen SNPs
#' 
#' @param snps_ character vector containing names of SNPs to be extracted in geno, for use
#' in breed probability estimation.
#' @param geno_ a matrix or data.frame with genotypes of test animals. SNPs in columns and 
#' individuals in rows. Genotypes must be coded as dosage of allele 'b' {0, 1, 2}. 
#' @param b_mat_ a data.frame that contains genotypes in rows and associated breed probabilities in
#' columns. This is the output of build_KBP(), build_MBP(), etc.
#' @param white_ a logical indicating whether "white" breeds (Yorkshire and Landrace, assuming present)
#' should be grouped together into a white breed probability.
#' @return a data.frame containing test animals in rows and breed probabilities in columns
solve_bp_generic <- function(snps_, geno_, b_mat_, white_) {
  
  if (!any(snps_ %in% colnames(geno_)))
    stop(paste0("Cannot detect KIT SNPs in provided genotypes. \n
                Check to make sure geno has SNP colnames, and that all SNPs are present: \n",
                snps_))
  
  # Subset and recode missing genotypes
  geno_subset <- geno_[, snps_]
  geno_subset[is.na(geno_subset)] <- "?"
  
  # Concatenate 7 SNP genotypes into a single string
  geno_subset <- apply(geno_subset, 1, paste, collapse = '')
  
  # For each test animal genotype, search through b_mat for appropriate KBPs
  results <- t(sapply(geno_subset, breed_lookup, b_mat_, white_freq = white_))
  
  return(results)
}