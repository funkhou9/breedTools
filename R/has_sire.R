#' Detect if each ID has a Father/Sire in pedigree
#' 
#' Genotyping data can also be provided to know whether
#'  an available sire has genotyping data available
#'  
#' @param id character naming ID of interest
#' @param ped data.frame consisting of pedigree information
#'  with columns {ID, Sire, Dam} in that order
#' @param geno matrix of SNP array data
#' @return boolean. TRUE if animal has a Father/Sire in the
#'  pedigree. If geno is provided, only returns TRUE if that
#'  Father/Sire has pedigree information
#' @export
has_sire <- function(id, ped, geno = NULL) {
  # Check for proper pedigree format
  if (!all(c("Sire", "Dam") %in% names(ped))) {
    stop("Can't recognize pedigree format. Requires 'Sire' and 'Dam' columns")
  }
  
  # Does test animal have a sire in ped and does the sire have
  if (id %in% ped[, 1]) {
    sire_id <- ped[ped[, 1] == id, "Sire"]
    
    # If geno is provided, is the sire genotyped?
    if(!is.null(geno)) {
      if (sire_id %in% rownames(geno)) {  
        return(TRUE)
      } else return(FALSE)
    
    # If no geno provided - Sire still in pedigree!
    } else {
      return(TRUE)
    }
    
  # If no sire in pedigree   
  } else return(FALSE)
}


