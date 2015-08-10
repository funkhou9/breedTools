#' Detect if each ID has a full-sib or half-sib in the pedigree, as specified
#' 
#' Genotyping data can also be provided to know whether
#'  the sib has genotyping data available
#'  
#' @param id character naming ID of interest
#' @param ped data.frame consisting of pedigree information
#'  with columns {ID, Sire, Dam} in that order
#' @return boolean. TRUE if animal has a sib of specified sib_type
#' @export
has_sib <- function(id, ped) {
  # Check for proper pedigree format
  if (!all(c("Sire", "Dam") %in% names(ped))) {
    stop("Can't recognize pedigree format. Requires 'Sire' and 'Dam' columns")
  }
  
  # Obtain parents of id
  if (id %in% ped[, 1]) {
    sire_id <- ped[ped[, 1] == id, "Sire"]
    dam_id <- ped[ped[, 1] == id, "Dam"]
  } else
    return (FALSE)
    
  
  # Get sib lengths
  sib_from_sire <- ped[ped[, 2] == sire_id, 1]
  sib_from_dam <- ped[ped[, 3] == dam_id, 1]
  
  # Check if sib lengths are greater than 1 (implying sibs are present)
  if (length(sib_from_sire) == 1 & length(sib_from_dam) == 1)
    return (FALSE)
  else 
    return (TRUE)
}