# Performs QPsolve on a simulated cross between a random animal in blist and a random
#  animal in plist.
#  
# @param Y numeric matrix of genotypes from all animals in population
# @param X numeric matrix of allele frequencies from reference animals
# @param blist list containing breed(s) and animal IDs for each breed
# @param plist list containing breed(s) and animal IDs for each breed
# @return list whose name is the animals IDs involved in the simulated cross
#  and value is a vector of QPsolve results.
QP_SimCross <- function(Ymat, Xmat, blist, plist) {

  # Grab random animals - ensuring they are from different breeds
  #   Obviously, each breed must have a different name.
  repeat {
    # Randomly selecting one breed from each blist and plist
    idx1 <- sample(length(blist), 1)
    idx2 <- sample(length(plist), 1)
    
    # Determine if breeds are the same - we need them to be different
    breeds <- c(blist[idx1], plist[idx2])
    
    # As long as breeds are different - use those breeds
    if (names(breeds)[1] != names(breeds)[2]) break
  }
  
  # Sample random animals from each breed - obtain ID and genotype
  name1 <- sample(breeds[1][[1]], 1)
  name2 <- sample(breeds[2][[1]], 1)
  
  gen1 <- Ymat[, name1]
  gen2 <- Ymat[, name2]
  
  # Creating randomly drawn cutpoints for use in shuffling
  ncuts <- rpois(1, 10)
  nsnp <- length(gen1)
  posc <- sort(sample(nsnp, ncuts))
  
  # The first item of posc must be 1
  if (min(posc) != 1) {posc <- c(1, posc)}
  
  # j and k are counters that will trade between 0 and 1
  # 	Resulting in alternating segments of gen1 and gen2
  #	present in genn.
  j <- 1
  genn <- gen1 * 0
  
  # Shuffling gen1 and gen2 into genn
  for (i in 2:length(posc)) {
    k <- 1 - j
    genn[posc[i - 1] : posc[i]] <- j * gen1[posc[i-1] : posc[i]] + k * gen2[posc[i - 1] : posc[i]]
    j <- k
  }
  
  # Solve composition for both original and synthetic genotypes
  qpnew <- QPsolve(genn, Xmat)
  qp1 <- QPsolve(gen1, Xmat)
  qp2 <- QPsolve(gen2, Xmat)
  
  # Compute the actual composition (of breed1) of the synthetic genotype 
  compos <- sum((posc[-1] - posc[-length(posc)])[seq(1, (length(posc) - 1), 2)]) / length(gen1)
  
  # Assemble results - index1, compos of par1, compos of par2, compos (breed1) of synthetic,
  #   compos (breed2) of synthetic, actual compos (breed1) of synthetic, R2
  result <- c(idx1[1], qp1[names(breeds)[1]], qp2[names(breeds)[2]], qpnew[names(breeds)[1]],
              qpnew[names(breeds)[2]], compos, qpnew["R2"])
  
  # Reset names, shape into list, and name the list
  names(result) <- NULL
  result <- list(result)
  names(result) <- paste(names(breeds)[1], ".", name1, "/", names(breeds)[2], ".", name2, sep='')
  
  return(result)
}

