#' Solve MC1R-Based Breed Probabilities (MBP)
#' 
#' Given a breed probability matrix (provided by \code{build_MBP()}), and a set of test animals
#' genotypes, use this function to derive breed probabilities for each test animal genotype based on
#' 8 SNPs within a 2MB window around the MC1R gene.
#' 
#' @param geno a matrix or data.frame with genotypes of test animals. SNPs in columns and 
#' individuals in rows. Genotypes must be coded as dosage of allele 'b' {0, 1, 2}. NOTE, only 8 SNPs
#' surrounding MC1R will actually be used.
#' @param b_mat a data.frame that contains genotypes in rows and associated breed probabilities in
#' columns. This is the output of build_MBP().
#' @param white a logical indicating whether "white" breeds (Yorkshire and Landrace, assuming present)
#' should be grouped together into a white breed probability. For MC1R analysis, this is set to FALSE
#' by default.
#' @return a data.frame containing test animals in rows and KBPs in columns
#' @export
solve_MBP <- function(geno, b_mat, white = FALSE) {
  
  solve_bp_generic(snps_ = mc1r_snps,
                   geno_ = geno,
                   b_mat_ = b_mat,
                   white_ = white)
  
}