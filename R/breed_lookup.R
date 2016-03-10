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
    white <- sum(hits_avg["landrace"], hits_avg["yorkshire"], hits_avg["landrace/yorkshire"])
    duroc_white <- sum(hits_avg["duroc/landrace"], hits_avg["duroc/yorkshire"])
    hampshire_white <- sum(hits_avg["hampshire/landrace"], hits_avg["hampshire/yorkshire"])
    
    hits_avg <- c(hits_avg[1:2], "white" = white, hits_avg[5], "duroc/white" = duroc_white,
                  "hampshire/white" = hampshire_white)
  }
  
  return(hits_avg)
}
