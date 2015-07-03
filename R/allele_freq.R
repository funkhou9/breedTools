#' Computes allele frequencies for specified populations given SNP array data
#'
#' @param geno matrix of genotypes coded as the dosage of allele B {0, 1, 2}
#'  with individuals in rows (named) and SNPs in columns (named)
#' @param populations list of named populations. Each population has a vector of IDs
#'  that belong to the population. Allele frequencies will be derived from all animals
#'   
#' @return data.frame consisting of allele_frequencies for populations (columns) for
#'  each SNP (rows)
#' @export
allele_freq <- function(geno, populations) {
  
  # Initialize returned df
  X <- matrix(NA, nrow = ncol(geno), ncol = length(populations))
  
  # Subset geno into different populations
  for (i in 1:length(populations)) {

    # Get name of ith item in the list (population name)
    pop_name <- names(populations[i])
    
    # Subset geno to only include genotypes of IDs in pop
    pop_geno <- geno[rownames(geno) %in% populations[[i]], ] 
    
    # Calculate allele frequencies
    al_freq <- colMeans(pop_geno, na.rm = TRUE) / 2
    
    # Add to X
    X[, i] <- al_freq
  }
  
  # Label X with populations and SNPs
  colnames(X) <- names(populations)
  rownames(X) <- colnames(geno)
  
  return(X)
}