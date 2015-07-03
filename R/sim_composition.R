#' sim_composition
#' 
#' Use synthetic offspring to test the accuracy of breed/ancestry calculations among test population
#' 
#' @param Y numeric matrix of genotypes from all animals in population
#' @param X numeric matrix of allele frequencies from reference animals
#' @param rep integer indicating how many repititions of the simulation to perform
#' @param par1 list of IDs catagorized by breed/population to be used as "Parent 1" in a simulated
#'  progeny test
#' @param par2 list of IDs catagorized by breed/population to be used as "Parent 2" in a simulated
#'  progeny test. If null, two parents from sim1 will be chosen
#' @import quadprog
#' @export
sim_composition <- function(Y, X, rep = 10000, par1, par2 = NULL) {
  
  # Functions require Y to be animals x snps. Transpose
  Y <- t(Y)
  
  # If no par2 provided, just re-use par1
  if (is.null(par2)) {
    par2 <- par1
  }
  
  # First, ensure par1 and par2 contain breeds of different names
  names_par1 <- names(par1)
  names_par2 <- names(par2)
  
  if (length(unique(c(names_par1, names_par2))) < 2)
    stop("There must be at least two breeds to choose from")
  
  # Use replicate to call QP_SimCross rep times
  sim <- replicate(rep, QP_SimCross(Y, X, par1, par2))
  
  # Modify to obtain final tabulated results 
  sim_tab <- do.call(rbind, sim)
  
  # Format with appropriate colnames
  colnames(sim_tab) <- c("idx", "Breed1.qp",
                         "Breed2.qp", "XBreed1.qp", "XBreed2.qp",
                         "ActBreed1", "rr")
  
  return(sim_tab)
}