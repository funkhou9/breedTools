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
#' surrounding the MC1R region (the output of build_MBP). (OPTIONAL) If supplied, it will
#' assume Hampshire animals are being evaluated in \code{geno}.
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
  
  # If mbp_ref is NULL, Yorkshire or Landrace animals are being evaluated.
  if (is.null(mbp_ref)) {
    # Join GWBC and KBP results keeping only values related to 'white' KBPs.
    results <- data.frame(gwbc_res[rownames(kbp_res), 1:4],
                          kbp_res[, c(3, 5, 6)],
                          rep("no", nrow(geno)),
                          stringsAsFactors = FALSE)
    colnames(results) <- c("GWBC_Duroc",
                           "GWBC_Hampshire",
                           "GWBC_Landrace",
                           "GWBC_Yorkshire",
                           "KBP_White",
                           "KBP_Duroc_White",
                           "KBP_Hampshire_White",
                           "KIT_inconclusive")  
  } else {
    # Here, Hampshire are presumably being evaluated. First estimate MC1R-based
    # probabilities
    mbp_res <- solve_MBP(geno, mbp_ref)
    
    # Join GWBC and KBP/MBP results. For KBP/MBP, provide purebred Hampshire
    # probability and the probability that is the highest (or second highest,
    # if purebred Hampshire is the highest)
    kbp_max_idx <- apply(kbp_res, 1, function(x) which(x == max(x)))
    kbp_max <- colnames(kbp_res)[kbp_max_idx]
    kbp_max_prob <- kbp_res[cbind(1:nrow(kbp_res), kbp_max_idx)]
    
    mbp_max_idx <- apply(mbp_res, 1, function(x) which(x == max(x)))
    mbp_max <- colnames(mbp_res)[mbp_max_idx]
    mbp_max_prob <- mbp_res[cbind(1:nrow(mbp_res), mbp_max_idx)]
    
    results <- data.frame(gwbc_res[rownames(kbp_res), 1:4],
                          kbp_res[, "Hampshire"],
                          kbp_max,
                          kbp_max_prob,
                          rep("no", nrow(geno)),
                          mbp_res[, "Hampshire"],
                          mbp_max,
                          mbp_max_prob,
                          rep("no", nrow(geno)),
                          stringsAsFactors = FALSE)
    colnames(results) <- c("GWBC_Duroc",
                           "GWBC_Hampshire",
                           "GWBC_Landrace",
                           "GWBC_Yorkshire",
                           "KBP_Hampshire",
                           "Max_KBP",
                           "KBP_Max",
                           "KIT_inconclusive",
                           "MBP_Hampshire",
                           "Max_MBP",
                           "MBP_Max",
                           "MC1R_inconclusive")  
  }
  
  # Find which animals, if any, have an incomplete set of SNPs in the KIT region
  kit_geno <- geno[, kit_snps]
  inconclusives <-
    rownames(kit_geno)[rowSums(is.na(kit_geno)) > 0]
  
  # Fill in which animals were "KIT inconclusive"
  results[inconclusives, "KIT_inconclusive"] <- "yes"
  
  if (!is.null(mbp_ref)) {
    mc1r_geno <- geno[, mc1r_snps]
    inconclusives <-
      rownames(mc1r_geno)[rowSums(is.na(mc1r_geno)) > 0]

    results[inconclusives, "MC1R_inconclusive"] <- "yes"
    return(results)
  }
  
  return (results)
}