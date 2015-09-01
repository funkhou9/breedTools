#' Obtain siblings of id, if any
#' 
#'  
#' @param id character naming ID of interest
#' @param ped data.frame consisting of pedigree information
#'  with columns {ID, Sire, Dam} in that order
#' @param gsib logical. If true, finds siblings related by parents
#'  and grand parents
#' @return character vector of sibling IDs
#' @export
get_sibs <- function(id, ped, gsib = TRUE) {
  
  # Check for proper pedigree format
  if (!all(c("Sire", "Dam") %in% names(ped)))
    stop("Can't recognize pedigree format. Requires 'Sire' and 'Dam' columns")
  
  # Obtain parents of id
  if (id %in% ped[, 1]) {
    sire_id <- ped[ped[, 1] == id, "Sire"]
    dam_id <- ped[ped[, 1] == id, "Dam"]
    
      
    if (gsib) {
      g_sire_1 <- ped[ped[, 1] == sire_id, "Sire"]
      g_sire_2 <- ped[ped[, 1] == dam_id, "Sire"]
      g_dam_1 <- ped[ped[, 1] == sire_id, "Dam"]
      g_dam_2 <- ped[ped[, 1] == dam_id, "Dam"]
    }
    
  } else
      return (NULL)
  
  # Get sib IDs
  sibs_from_sire <- ped[ped[, 2] == sire_id, 1]
  sibs_from_dam <- ped[ped[, 3] == dam_id, 1]
  
  # Get siblings from grandparents
  if (gsib) {
    g_sire_proj <- unique(c(ped[ped[, 2] == g_sire_1, 1],
                            ped[ped[, 2] == g_sire_2, 1]))
    
    g_dam_proj <- unique(c(ped[ped[, 3] == g_dam_1, 1],
                           ped[ped[, 3] == g_dam_2, 1]))
                    
    
    # Combine "grand siblings"
    g_proj <- unique(c(g_sire_proj, g_dam_proj))
    
    
    gsibs_from_sire <- ped[ped[, 2] %in% g_proj, 1]
    gsibs_from_dam <- ped[ped[, 3] %in% g_proj, 1]
  }
  

  # Obtain sibs! (including id)
  if (gsib)
    siblings <- unique(c(sibs_from_sire,
                         sibs_from_dam,
                         gsibs_from_sire,
                         gsibs_from_dam))
  else
    siblings <- unique(c(sibs_from_sire,
                         sibs_from_dam))
  
  siblings
}