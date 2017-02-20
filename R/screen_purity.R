#' Calculate both GWBC and KBP for a test animal
#' 
#' GWBC (Genome-wide breed composition) and KBP (KIT-based breed probabilities) are calculations
#' done together to assess the breed composition of a test animal and the probability that it is
#' purebred in the KIT region of the genome, resposible for influencing color segregation patters.
#' This function will assemble all GWBC and KBP results for any number of genotyped animals.
#' 
#' @param geno matrix containing genotyping data (dosage of allele B) for animals (in rows) for
#' any number of SNPs (in columns)
#' @param gwbc_ref data.frame containing allele frequencies for any number of SNPs (in rows) for
#' any number of reference populations (in columns). See an example with `data("GWBC_ref_A")`.
#' @param kbp_ref data.frame containing genotype probabilities for various multi-locus genotypes
#' surrounding the KIT region (the output of build_KBP).
#' @param mbp_ref data.frame containing genotype probabilities for various multi-locus genotypes
#' surrounding the MC1R region (the output of build_MBP). (OPTIONAL)
#' @return data.frame with various GWBC and KBP values for each animal in geno
#' @export
screen_purity <- function(geno,
                          gwbc_ref,
                          kbp_ref,
                          mbp_ref = NULL) {
  
  # Solve GWBC
  gwbc_res <- solve_composition(geno, gwbc_ref)
  
  # Solve KBP
  kbp_res <- solve_KBP(geno, kbp_ref)
  
  # Join GWBC and KBP results keeping only values related to 'white' KBPs.
  results <- data.frame(gwbc_res[rownames(kbp_res), 1:4],
                        kbp_res[, c(3, 5, 6)],
                        kit_analysis_inconclusive = rep("no", nrow(geno)),
                        stringsAsFactors = FALSE)
  colnames(results) <- c("GWBC_Duroc",
                         "GWBC_Hampshire",
                         "GWBC_Landrace",
                         "GWBC_Yorkshire",
                         "KBP_White",
                         "KBP_Duroc_White",
                         "KBP_Hampshire_White",
                         "KIT_inconclusive")
  
  # Find which animals, if any, have an incomplete set of SNPs in the KIT region
  kit_geno <- geno[, kit_snps]
  inconclusives <-
    kit_geno[rowSums(is.na(kit_geno)) > 0, ] %>%
    rownames()
  
  # Fill in which animals were "KIT inconclusive"
  results[inconclusives, 8] <- "yes"
  
  # If provided a reference panel for MBP, solve MBP and attach MBP results.
  # Repeat the process of identifying 'MC1R_inconclusive' results.
  if (!is.null(mbp_ref)) {
    mbp_res <- solve_MBP(geno, mbp_ref)
    results_mbp <- cbind(results,
                         mbp_res[rownames(kbp_res), c(2, 5, 8, 9)],
                         "MC1R_inconclusive" = rep("no", nrow(geno)),
                         stringsAsFactors = FALSE)
    colnames(results_mbp) <- c(colnames(results),
                               "MBP_Hampshire",
                               "MBP_Hampshire_Duroc",
                               "MBP_Hampshire_Landrace",
                               "MBP_Hampshire_Yorkshire",
                               "MC1R_inconclusive")
    
    mc1r_geno <- geno[, mc1r_snps]
    inconclusives <-
      mc1r_geno[rowSums(is.na(mc1r_geno)) > 0, ] %>%
      rownames()

    results_mbp[inconclusives, "MC1R_inconclusive"] <- "yes"
    return(results_mbp)
  }
  
  return (results)
}