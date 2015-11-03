#' Given a list of ids to choose from and a pedigree, returns a subset of the list that are as
#'  distantly related as possible.
#' 
#' @param ids vector of ids to choose from
#' @param ped data.frame containing pedigree information. Should contain columns "ID" "Sire" and "Dam"
#'  in that order.
#' @param coeff numeric giving maximum coefficient of relatedness among chosen animals.
#' @import kinship2
#' @export
pick_unrelated <- function(ids, ped, coeff = 0.0625) {
  
  # Create A matrix (pedigree relationship matrix)
  A <- kinship(ped[, 1], ped[, 2], ped[, 3])
  A_red <- A[rownames(A) %in% ids, colnames(A) %in% ids]
  
  # Sort matrix according to which animals have the least relationships among IDs
  idx <- c()
  for (i in 1:nrow(A_red)) {
    idx <- c(idx, sum(A_red[i, -i] > 0))
  }
  
  # Run through rows of sorted matrix, picking unrelated individuals
  nochoose <- c()
  pick <- c()
  for (i in order(idx)) {
    kins <- A_red[i, ]
    
    # current animal
    id <- rownames(A_red)[i]
    
    # is the current animal related to an animal chosen so far? or has it already been picked?
    if (id %in% nochoose)
      next
    else
      pick <- c(pick, id)
    
    # any related animals added to nochoose
    rel <- names(kins[kins > coeff])
    nochoose <- c(nochoose, rel)
  }
  
  return (list("A_mat" = A_red, "IDs" = pick))
}

