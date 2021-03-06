#' Solve KIT-Based Breed Probabilities (KBP)
#' 
#' Given a breed probability matrix (provided by \code{build_KBP()}), and a set of test animals
#' genotypes, use this function to derive breed probabilities for each test animal genotype based on
#' 7 SNPs around the KIT gene.
#' 
#' @param geno a matrix or data.frame with genotypes of test animals. SNPs in columns and 
#' individuals in rows. Genotypes must be coded as dosage of allele 'b' {0, 1, 2}. NOTE, only 7 SNPs
#' surrounding KIT will actually be used.
#' @param b_mat a data.frame that contains genotypes in rows and associated breed probabilities in
#' columns. This is the output of build_KBP().
#' @param white a logical indicating whether "white" breeds (Yorkshire and Landrace, assuming present)
#' should be grouped together into a white breed probability.
#' @return a data.frame containing test animals in rows and KBPs in columns
#' @export
solve_KBP <- function(geno, b_mat, white = TRUE) {
  
  solve_bp_generic(snps_ = kit_snps,
                   geno_ = geno,
                   b_mat_ = b_mat,
                   white_ = white)
  
}



