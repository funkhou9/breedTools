#' Lookup haplotypes from a breed probability matrix
#' 
#' Breed probability matrices can be constructed with the help of haplo_comp(), they contain 
#' genotypes for a region in rows and probabilities that that genotype belongs to a breed / breed hybrid
#' in the columns.
#' 
#' @param geno string containing genotypes for a region represented as (ex. 2011210?1211?)
#' @param bmat a breed probability matrix as described from above
#' @param white_freq logical. If true, combines breed probabilities of Yorkshire and Landrace breeds
#'  when those are present.
#' @return probabilities of being particular breeds based on stretch of genotype information provided
#'  in geno argument.
breed_lookup <- function(geno, bmat, white_freq = FALSE) {
  
  # Create regex expression for geno
  geno_rx <- glob2rx(geno); names(geno_rx) <- names(geno)
  
  # Identify row(s) of bmat that correspond to geno
  rows <- grep(geno_rx, rownames(bmat))
  
  # Hits
  hits <- bmat[rows, ]
  
  # Average if more than one hit was found
  if (length(rows) > 1) {
    hits_avg <- colMeans(hits)	
  } else {
    hits_avg <- hits
  }
  
  # If breed "white" breed composition is desired
  #	instead of yorkshire or landrace
  if (white_freq) {
    white <- sum(hits_avg["Landrace"], hits_avg["Yorkshire"], hits_avg["Landrace/Yorkshire"])
    duroc_white <- sum(hits_avg["Duroc/Landrace"], hits_avg["Duroc/Yorkshire"])
    hampshire_white <- sum(hits_avg["Hampshire/Landrace"], hits_avg["Hampshire/Yorkshire"])
    
    hits_avg <- c(hits_avg[1:2], "White" = white, hits_avg[5], "Duroc/White" = duroc_white,
                  "Hampshire/White" = hampshire_white)
  }
  
  return(hits_avg)
}
