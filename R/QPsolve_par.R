# Performs breed composition prediction based on knowledge 
# of sire genotypes. Done by predicting only maternally
# inherited haplotype instead of target animal genotype. 
# Requires pedigree file.
# 
# @param id character of a single target animal
# @param Y numeric matrix of genotypes from all animals in population
# @param X numeric matrix of allele frequencies from reference animals
# @param ped data.frame of pedigree info for animal Y
# @param p numeric indicating number of breeds represented in X
# @param names character names of breeds
# @param mia logical. If true, rather than returning breed compsosition estimates, return
#   inferred maternally inherited alleles for each locus.
# @param sire logical. If true, provides sire genotypes rather than breed composition results.
# @param dam logical. If true, provides dam genotypes rather than breed composition results.
# @return data.frame of breed composition estimates
# @import quadprog
# @export
QPsolve_par <- function(id, Y, X, ped,
                        p = 4,
                        names = c("Duroc", "Hampshire", "Landrace", "Yorkshire"),
                        mia = FALSE,
                        sire = FALSE,
                        dam = FALSE) {
  
  # Check for proper pedigree format
  if (!all(c("Sire", "Dam") %in% names(ped))) {
    stop("Can't recognize pedigree format. Requires 'Sire' and 'Dam' columns")
  }
  

  # Does test animal have a sire in ped and does the sire have
  # 	a genotype in Y? (it must)
  if (id %in% ped[, 1]) {
    
    # If dam genotypes are requested, obtain and return those
    if (dam) {
      dam_id <- ped[ped[, 1] == id, "Dam"]
      
      if (dam_id %in% colnames(Y)) {
        dam_geno <- Y[, colnames(Y) == dam_id]
        
        dam_geno <- t(as.data.frame(dam_geno))
        rownames(dam_geno) <- id
        
        return (dam_geno)
      }
    }
    
    
    
    # Assuming dam geno isn't requested, proceed with using sire genotype
    sire_id <- ped[ped[, 1] == id, "Sire"]
    
    if (sire_id %in% colnames(Y)) {
    
      # Get genotype of target animal and convert missing genotypes
      #	to the appropriate missing character.
      geno <- Y[, colnames(Y) == id]
      geno[is.na(geno)] <- '-'
      
      # Get sire genotype and convert missing characters
      sire_geno <- Y[, colnames(Y) == sire_id]
      
      # If sire genotype is requested, return that istead of proceeding
      if (sire) {
        sire_geno <- t(as.data.frame(sire_geno))
        rownames(sire_geno) <- id
        
        return (sire_geno)
      }
      
      sire_geno[is.na(sire_geno)] <- '-'
      
      # Perform calculation of maternally inherited haplotype
      #	for each id genotype with mapply()
      mat_hap <- mapply(mat_allele, geno, sire_geno)
      
      # If only maternally inherited alleles are desired, convert
      #   modified Y to a df[1, ] with name equal to supplied id
      if (mia) {
        mat_hap <- t(as.data.frame(mat_hap))
        rownames(mat_hap) <- id
        
        return (mat_hap)
      }
      
      # Remove non-numerics from maternal haplotype
      # ("?") and convert to numeric. 
      Ymod <- mat_hap[mat_hap != '?']
      Ymod_num <- as.numeric(Ymod)
      
      # Ensure names are kept after conversion to numeric
      names(Ymod_num) <- names(Ymod)
      
      # Perform QPsolve with the new modified Ymod and X
      result <- QPsolve(Ymod_num, X)
      
      # Convert result to a single row df and give name of id
      result <- t(as.data.frame(result))
      rownames(result) <- id
      
      return(result)
    }
  }
}
