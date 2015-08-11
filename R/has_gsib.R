#' Detect if ID has a 'great sibling' in the pedigree, as specified
#' 
#' 'great sibling' would mean siblings related by their grandparents
#'  
#' @param id character naming ID of interest
#' @param ped data.frame consisting of pedigree information
#'  with columns {ID, Sire, Dam} in that order
#' @return boolean. TRUE if animal has a great sibling
#' @export
has_gsib <- function(id, ped) {
  # Check for proper pedigree format
  if (!all(c("Sire", "Dam") %in% names(ped))) {
    stop("Can't recognize pedigree format. Requires 'Sire' and 'Dam' columns")
  }
  
  # Obtain parents of id
  if (id %in% ped[, 1]) {
    sire_id <- ped[ped[, 1] == id, "Sire"]
    dam_id <- ped[ped[, 1] == id, "Dam"]
    
    g_sire_id <- ped[ped[, 1] == sire_id, "Sire"]
    g_dam_id <- ped[ped[, 1] == sire_id, "Dam"]
  } else
    return (FALSE)
  
  
  # Get sib lengths
  g_sire_proj <- ped[ped[, 2] == g_sire_id, 1]
  g_dam_proj <- ped[ped[, 3] == g_dam_id, 1]
  
  sib_sires <- ped[ped[, 2] %in% g_sire_proj, 1]
  sib_dams <- ped[ped[, 2] %in% g_dam_proj, 1]
  
  
  
  # Check if sib lengths are greater than 1 (implying sibs are present)
  if (length(sib_sires) == 1 & length(sib_dams) == 1)
    return (FALSE)
  else 
    return (TRUE)
}