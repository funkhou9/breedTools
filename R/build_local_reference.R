#' Build reference panel of breed probabilities 
#' 
#' Generalized function to build a reference panel, used to estimate breed probabilities for a set 
#' of multi-locus genotypes from a specified region.
#' 
#' @param chr character providing the chromosome requiring phasing
#' @param markers character vector containing names of probes used to construct multi-locus haplotypes
#' @param geno_ a matrix or data.frame with genotypes of reference animals. SNPs in columns and 
#' individuals in rows. Genotypes must be coded as dosage of allele 'b' {0, 1, 2}. 
#' @param map_ a data.frame containing SNP map information for each SNP present in geno
#' @param ped_ a data.frame pedigree providing family information for each individual in geno. The first
#' column of the pedigree is for ID, second is for sire/father ID, third is for dam/mother ID, fourth
#' is sex of ID {"M" or "F"}. If ids in geno are not present in pedigree, they will be added to the
#' end with missing parent information and a "M" sex. Correct sex information should only be required
#' when imputing/phasing the sex chromosomes.
#' @param path_ a character represting the path to the FImpute binary. If omitted, assumes FImpute binary
#' resides along PATH.
#' @param groups_ a list of character vectors with the names of IDs meant to be processed as groups. 
#' Names of list should be the names of the groups, with each element containing IDs.
#' @param keep_fimpute_ logical if TRUE, output files generated from snpTools::fimpute_run are kept. 
#' By default, build_KBP will use fimpute_run to build haplotypes, but as soon as haplotypes are 
#' extracted, fimpute output will be deleted so not to clutter hard drive space.
#' @param parent_ logical if TRUE, breed probabilities or haplotypes will be returned only from
#' parents and not from progeny,
#' @param reference_ logical use TRUE if you want `groups` to be used to specify the reference animals
#' used by FImpute. Not compatable with parent TRUE.
#' @return a data.frame containing breed probabilities for each genotype made possible from the
#' reference animal haplotypes
build_local_reference <- function(chr,
                                  markers,
                                  geno_,
                                  map_,
                                  ped_ = NULL,
                                  path_ = NULL,
                                  groups_ = NULL,
                                  keep_fimpute_ = FALSE,
                                  parent_ = FALSE,
                                  reference_ = FALSE) {
  
  if (is.null(groups_))
    stop("build_local_reference() is designed to be used with a groups arugment, corresponding to a list of reference breeds")
  
  # Obtain a single string containing all chromosomes to ignore (all chromosomes except the one to phase)
  chr_to_ignore <- paste(paste(1:21)[as.numeric(paste0("-", chr))], sep = "", collapse = " ")
  
  # Run FImpute to obtain reference haplotypes, but only along specified chromosome
  invisible({
    snpTools::fimpute_run(geno = geno_,
                          map = map_,
                          ped = ped_,
                          path = path_,
                          groups = groups_,
                          exclude_chr = chr_to_ignore,
                          parent = parent_,
                          reference = reference_)
  })
  message("Obtaining reference haplotypes from FImpute")
  
  # Obtain haplotypes and haplotype frequencies from each reference breed
  hap_list <- list()
  hap_freq_list <- list()
  for (i in 1:length(groups_)) {
    
    # Get SNP names
    snps <- read.table(paste0(names(groups_)[i], "_fimpute_run/hap_library.txt"),
                       nrow = 1) %>%
      unlist() %>%
      as.character()
    
    # Read haplotype file as a single string
    hap_list[[i]] <- read.table(paste0(names(groups_)[i], "_fimpute_run/hap_library.txt"),
                                skip = 1, 
                                colClasses = "character")
    
    # Convert to a vector with snpTools::string_to_vec()
    hap_list[[i]] <- t(apply(hap_list[[i]], 1, snpTools::string_to_vec))
    
    # Subset chosen SNPs
    colnames(hap_list[[i]]) <- snps
    hap_list[[i]] <- hap_list[[i]][, markers]
    
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
  
  if (!keep_fimpute_)
    system('rm -rf *fimpute_run')
  
  ## if (haplotypes) {
  # Merge haplotype frequencies for optional return of reference haplotype information
  hap_data <- hap_freq_list[[1]]
  for (i in 2:length(hap_freq_list)) {
    hap_data <- merge(hap_data, hap_freq_list[[i]], by = "haplotypes", all = TRUE)
  }
  colnames(hap_data) <- c("Haplotype", names(groups_))
  ## return(hap_data)
  ## }
  
  
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
                         paste0(names(groups_)[i],
                                '/',
                                names(groups_[j])))
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
  # g_data should be a data.frame with (1 + length(groups_) + n_between) columns, corresponding to the 
  # genotypes + number of "within breed" and "between breed" combinations
  rownames(g_data) <- g_data[, 1]
  g_data <- g_data[, 2 : (1 + length(groups_) + n_between)]
  
  # Set column names from groups_ input
  colnames(g_data) <- c(names(groups_), between_names)
  
  # Calculate b_mat - breed probabilities
  b_data <- t(apply(g_data, 1, function(x) x / sum(x)))
  
  # Return a list containing haplotype frequencies and resulting breed probabilities
  return(list("BP" = b_data, "haplotype_freq" = hap_data))
}
