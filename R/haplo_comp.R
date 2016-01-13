#' Function to calculate haplotype frequencies and use
#'	them to estimate genotype probabilities
#' 
#' @param x A data.frame that contains haplotypes for individuals
#' (rows) of a population, for a number of SNPs (columns). Haplotypes coded
#' 0 or 1 as type character.
#' @param y A data.frame that contains haplotypes for individuals
#' (rows) of a population, for a number of SNPs (columns). Haplotypes coded
#' 0 or 1 as type character.
#' @param haplotypes boolean. If true, return haplotype frequencies of population x.
#'  If false, return breed probability matrix (posterior probabilities for each breed,
#'  rows, to have particular genotype, columns)
#' @return Either haplotype frequencies or breed probability matrix
#' @export
haplo_comp <-  function(x, y = NULL, haplotypes = FALSE) {
  
  # Check inputs. Must be a data.frame or matrix
  if (class(x) != "matrix" & class(x) != "data.frame") {
    stop("x must be a data.frame or matrix")
  } 
  
  if (!is.null(y)) {
    if (class(y) != "matrix" & class(y) != "data.frame") {
      stop("y must be a data.frame or matrix")
    } 
    
  } else {
    # If y is not supplied, set y and x to be the same
    y <- x
  }
  # Consolidate haplotypes into a single string
  hapx <- apply(x, 1, paste, collapse = '')
  hapy <- apply(y, 1, paste, collapse = '')
  
  # Constuct haplotype frequencies
  hapx_count <- table(hapx)
  hapy_count <- table(hapy)
  hapx_freq <- c(hapx_count / sum(hapx_count))
  hapy_freq <- c(hapy_count / sum(hapy_count))
  
  # If haplotypes is TRUE, simply return haplotype frequencies
  #	Only returns haplotype of x if y is provided
  if (haplotypes) return(hapx_freq)
  
  # Compute all possible genotype probabilities using outer product of
  #	hap_freq * hap_freq'
  G <- outer(hapx_freq, hapy_freq, FUN = "*")
  
  # Similarly, compute all possible genotypes by summing the each combination of
  #	haplotypes
  mat <- outer(as.numeric(rownames(G)), as.numeric(colnames(G)), FUN = "+")
  
  # Leading 0s will be lost unless genotypes are padded with 0s
  mat <- sprintf('%07d', mat)
  
  # Which genotypes are duplicated?
  dup_geno <- names(table(mat))[table(mat) > 1]
  
  # Which genotypes aren't duplicated? Should just be diagonal
  single_geno <- as.character(mat[!mat %in% dup_geno])
  
  # For each geno - find the position in mat and sum the 
  #	corresponding positions in G 
  sum_Gi <- function(geno) {
    idx <- which(mat %in% geno, arr.ind = TRUE)
    genotype <- sum(G[idx])
    return(genotype)
  }
  
  # Apply sum_Gi across dup_geno - a vector of genotypes which
  #	have duplicates
  dup_geno_sum <- sapply(dup_geno, sum_Gi)
  
  # Append genotypes which only appear once
  single_geno_sum <- sapply(single_geno, sum_Gi)
  result <- c(dup_geno_sum, single_geno_sum)
  
  return(result)
}