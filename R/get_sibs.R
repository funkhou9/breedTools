#' Obtain siblings of id, if any
#' 
#'  
#' @param id character naming ID of interest
#' @param ped data.frame consisting of pedigree information
#'  with columns {ID, Sire, Dam} in that order
#' @return character vector of sibling IDs
#' @export
get_sibs <- function(id, ped) {
  
  # Check for proper pedigree format
  if (!all(c("Sire", "Dam") %in% names(ped)))
    stop("Can't recognize pedigree format. Requires 'Sire' and 'Dam' columns")
  
  # Obtain parents of id
  if (id %in% ped[, 1]) {
    sire_id <- ped[ped[, 1] == id, "Sire"]
    dam_id <- ped[ped[, 1] == id, "Dam"]
  } else
      return (NULL)
  
  # Get sib IDs
  sibs_from_sire <- ped[ped[, 2] == sire_id, 1]
  sibs_from_dam <- ped[ped[, 3] == dam_id, 1]
  
  # Obtain sibs! (including id)
  unique(c(sibs_from_sire, sibs_from_dam))
}