#' Constructs pedigree information for select individuals
#' 
#' @param ids character vector containing individuals of interest
#' @param ped data.frame of full pedigree to subset from
#' @param gmax numeric max number of generations to go back
#' @return data.frame of subsetted pedigree
#' @export
ped_build <- function(ids, ped, gmax=1000) {
  
  selectIDs <- ids
  g <- 0 
  
  repeat {
    selectPedigree <- ped[ped[, 1] %in% selectIDs, 1:3] 
    newSelectIDs <- unique(unlist(selectPedigree)) 
    print(nrow(selectPedigree)) 
    if(length(newSelectIDs) == length(selectIDs)) break
    if(g == gmax) break
    selectIDs <- newSelectIDs
    g <- g + 1
  }
  return(selectPedigree)
}