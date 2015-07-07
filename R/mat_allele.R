# Computes maternally inherited allele, given
# genotype of progeny and sire for a single locus.

# @param geno character {0, 1, 2} dosage of allele B
#		for target animal for each SNP
# @param sire_geno character {0, 1, 2} dosage of allele B
#		for sire of target animal for each SNP
# @param missing character used to represent missing values
# @return character {0, 0.5, 1, "?", NA} representing maternally
#		inherited allele. "?" indicates maternal inconsistancy.
#		0.5 indicates a maternally inherited allele of either 0 or 1.
#   NA if the progeny genotype is missing.
mat_allele <- function(geno, sire_geno, missing = '-') {
  
  # Two possibilities for maternal allele if sire_geno
  #	dosage is 0
  if (sire_geno == 0) {
    if (geno == 0)
      mat_allele <- 0 
    if (geno == 1)
      mat_allele <- 1
    if (geno == 2)
      mat_allele <- '?'
    if (geno == missing)
      mat_allele <- NA
  }
  
  # Three possibilities for maternal allele if sire_geno
  #	dosage is 1. Ambiguous mat_allele coded as 0.5	
  if (sire_geno == 1) {
    if (geno == 0)
      mat_allele <- 0
    if (geno == 1)
      mat_allele <- 0.5
    if (geno == 2)
      mat_allele <- 1
    if (geno == missing)
      mat_allele <- NA
  }
  
  # Two possibilities for maternal allele if sire_geno
  #	dosage is 2
  if (sire_geno == 2) {
    if (geno == 0)
      mat_allele <- '?'
    if (geno == 1)
      mat_allele <- 0
    if (geno == 2)
      mat_allele <- 1
    if (geno == missing)
      mat_allele <- NA
  }
  
  # If sire_geno is unknown, we return 'expected' mat_allele
  # If neither target or sire geno is known, returns NA
  #	(nothing is known!)
  if (sire_geno == missing) {
    if (geno == 0)
      mat_allele <- 0
    if (geno == 1)
      mat_allele <- 0.5
    if (geno == 2)
      mat_allele <- 1
    if (geno == missing)
      mat_allele <- NA
  }
  
  return (mat_allele)
}