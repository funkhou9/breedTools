#' Build reference panel for KBP
#' 
#' KBP (Kit-based breed probabilities) utilize a reference panel of purebred animals' genotypes.
#' The genotypes are phased, and the resulting haplotypes are used to construct haplotype frequences,
#' genotype frequencies given breed \eqn{P(g | b)}, and finally breed probabilities given a genotype
#' \eqn{P(b | g)}. This function is used to construct a reference panel of breed probabilities.
#' 
#' @param geno a matrix or data.frame with genotypes of reference animals. SNPs in columns and 
#' individuals in rows. Genotypes must be coded as dosage of allele 'b' {0, 1, 2}. 
#' @param map a data.frame containing SNP map information for each SNP present in geno
#' @param ped a data.frame pedigree providing family information for each individual in geno. The first
#' column of the pedigree is for ID, second is for sire/father ID, third is for dam/mother ID, fourth
#' is sex of ID {"M" or "F"}. If ids in geno are not present in pedigree, they will be added to the
#' end with missing parent information and a "M" sex. Correct sex information should only be required
#' when imputing/phasing the sex chromosomes.
#' @param path a character represting the path to the FImpute binary. If omitted, assumes FImpute binary
#' resides along PATH.
#' @param groups a list of character vectors with the names of IDs meant to be processed as groups. 
#' Names of list should be the names of the groups, with each element containing IDs.
#' @param haplotypes logical indicating whether build_KBP should return haplotype frequencies (TRUE)
#' or breed probabilities (FALSE).
#' @param keep_fimpute logical if TRUE, output files generated from snpTools::fimpute_run are kept. 
#' By default, build_KBP will use fimpute_run to build haplotypes, but as soon as haplotypes are 
#' extracted, fimpute output will be deleted so not to clutter hard drive space.
#' @import snpTools
#' @import magrittr
#' @export
#' @return a data.frame containing breed probabilities for each genotype made possible from the
#' reference animal haplotypes
build_KBP <- function(geno,
                      map,
                      ped = NULL,
                      path = NULL,
                      groups = NULL,
                      haplotypes = FALSE,
                      keep_fimpute = FALSE) {
  
  if (is.null(groups))
    stop("build_KBP() is designed to be used with a groups arugment, corresponding to a list of reference breeds")
  
  # SNPs to be extracted from haplotypes
  kit_snps <- c("ALGA0047798",
                "ALGA0047807",
                "ALGA0047809",
                "ALGA0102731",
                "ALGA0115258",
                "ALGA0123881",
                "MARC0034580")
  
  # Run FImpute to obtain reference haplotypes, but only along chromosome 8
  invisible({
    snpTools::fimpute_run(geno = geno,
                          map = map,
                          ped = ped,
                          path = path,
                          groups = groups,
                          exclude_chr = "1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21")
  })
  message("Obtaining reference haplotypes from FImpute")
  
  # Obtain haplotypes and haplotype frequencies from each reference breed
  hap_list <- list()
  hap_freq_list <- list()
  for (i in 1:length(groups)) {
    
    # Get SNP names
    snps <- read.table(paste0(names(groups)[i], "_fimpute_run/hap_library.txt"),
                       nrow = 1) %>%
              unlist() %>%
                as.character()
    
    # Read haplotype file as a single string
    hap_list[[i]] <- read.table(paste0(names(groups)[i], "_fimpute_run/hap_library.txt"),
                                  skip = 1, 
                                  colClasses = "character")
    
    # Convert to a vector with snpTools::string_to_vec()
    hap_list[[i]] <- t(apply(hap_list[[i]], 1, snpTools::string_to_vec))
    
    # Subset KIT SNPs
    colnames(hap_list[[i]]) <- snps
    hap_list[[i]] <- hap_list[[i]][, kit_snps]
    
    # Remove any haplotypes with missing data (marked by the int 5)
    hap_list[[i]] <- hap_list[[i]][apply(hap_list[[i]] != "5", 1, all), ]
    
    # Convert coding to dosage of allele 2
    hap_list[[i]] <- t(apply(hap_list[[i]], 1, function(x) as.numeric(x) - 1))
    hap_list[[i]] <- t(apply(hap_list[[i]], 1, as.character))
    
    # Compute haplotype frequencies and tabulate results
    hap_freq_list[[i]] <- haplo_comp(hap_list[[i]], haplotypes = TRUE)
    hap_freq_list[[i]] <- data.frame("haplotypes" = names(hap_freq_list[[i]]),
                                     "frequency" = hap_freq_list[[i]],
                                     stringsAsFactors = FALSE)
  }
  
  if (!keep_fimpute)
    system('rm -rf *fimpute_run')
  
  if (haplotypes) {
    # Merge haplotype frequencies for optional return of reference haplotype information
    hap_data <- hap_freq_list[[1]]
    for (i in 2:length(hap_freq_list)) {
      hap_data <- merge(hap_data, hap_freq_list[[i]], by = "haplotypes", all = TRUE)
    }
    colnames(hap_data) <- c("Haplotype", names(groups))
    return(hap_data)  
  }
  
  
  # Compute "within breed" genotype probabilities
  within_breed <- list()
  for (i in 1:length(hap_list)) {
    hap_freq_i <- haplo_comp(hap_list[[i]])
    within_breed[[i]] <- data.frame("G" = names(hap_freq_i),
                                    "pG" = hap_freq_i,
                                    stringsAsFactors = FALSE)
  }
  
  # Compute "between breed" genotype probabilities (with counter to count the number of combinations)
  between_breed <- list()
  n_between <- 0
  between_names <- c()
  for (i in 1:(length(hap_list) - 1)) {
    for (j in (i + 1):length(hap_list)) {
      hap_cross_freq_ij <- haplo_comp(hap_list[[i]], hap_list[[j]])
      between_breed <- 
        c(between_breed, 
          list(data.frame("G" = names(hap_cross_freq_ij),
                          "pG" = hap_cross_freq_ij,
                          stringsAsFactors = FALSE)))
      n_between <- n_between + 1
      between_names <- c(between_names,
                         paste0(names(groups)[i],
                                '/',
                                names(groups[j])))
    }
  }
  
  # Combine probabilities and join by "G" (genotype)
  g_mat_list <- c(within_breed, between_breed)
  g_data <- g_mat_list[[1]]
  for (i in 2:length(g_mat_list)) {
    g_data <- merge(g_data, g_mat_list[[i]], by = "G", all = TRUE)
  }
  
  # Change NAs in g_data to 0
  g_data[is.na(g_data)] <- 0
  
  # Adjust g_data so that genotypes (first column) are in rownames
  # g_data should be a data.frame with (1 + length(groups) + n_between) columns, corresponding to the 
  # genotypes + number of "within breed" and "between breed" combinations
  rownames(g_data) <- g_data[, 1]
  g_data <- g_data[, 2 : (1 + length(groups) + n_between)]
  
  # Set column names from groups input
  colnames(g_data) <- c(names(groups), between_names)
  
  # Calculate b_mat - breed probabilities
  b_data <- t(apply(g_data, 1, function(x) x / sum(x)))
  
  return(b_data)
}



