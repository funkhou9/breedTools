#' Detect if each ID has progeny in pedigree and how many.
#' 
#' Genotyping data can also be provided to know whether
#'  an available sire has genotyping data available
#'  
#' @param id character naming ID of interest
#' @param ped data.frame consisting of pedigree information
#'  with columns {ID, Sire, Dam} in that order
#' @param geno matrix of SNP array data
#' @return integer indicating how many progeny the id of interest
#'  has, either in the provided pedigree or in the pedgree with
#'  genotyping data. If none found, returns NA.
#' @export
has_progeny <- function(id, ped, geno = NULL) {
  
  # Check for proper pedigree format
  if (!all(c("Sire", "Dam") %in% names(ped))) {
    stop("Can't recognize pedigree format. Requires 'Sire' and 'Dam' columns")
  }
  
  # Check if id is a sire
  if (id %in% ped[, "Sire"]) {
    proj <- ped[ped[, "Sire"] == id, 1]
  
  # Else check if id is a dam
  } else if (id %in% ped[, "Dam"]) {
    proj <- ped[ped[, "Dam"] == id, 1]
  
  # Case when no projeny present
  } else {
    proj <- NULL
    return(NA)
  }
  
  # If genotyping data is provided - return the number
  #   of progeny genotyped
  if (!is.null(geno)) {
    return (sum(proj %in% rownames(geno)))  
  
  # Otherwise just return the number of progeny
  } else {
    return (length(proj))
  }
}