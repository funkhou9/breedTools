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
  
  # Figure out how to share this data with build_KBP, using an environment?
  kit_snps <- c("ALGA0047798",
                "ALGA0047807",
                "ALGA0047809",
                "ALGA0102731",
                "ALGA0115258",
                "ALGA0123881",
                "MARC0034580")
  
  if (!any(kit_snps %in% colnames(geno)))
    stop("Cannot detect KIT SNPs in provided genotypes. \n
         Check to make sure geno has SNP colnames, and that all SNPs are present: \n
         ALGA0047798, ALGA0047807, ALGA0047809, ALGA0102731, ALGA0115258, ALGA0123881, MARC0034580")
  
  # Subset and recode missing genotypes
  kit_geno <- geno[, kit_snps]
  kit_geno <- kit_geno[is.na(kit_geno)] <- "?"
  
  # Concatenate 7 SNP genotypes into a single string
  kit_geno <- apply(kit_geno, 1, paste, collapse = '')
  
  # For each test animal genotype, search through b_mat for appropriate KBPs
  results <- t(sapply(kit_geno, breed_lookup, b_mat, white_freq = white))
  
  return(results)
}



