#' Build reference panel for KBP
#' 
#' KBP (Kit-based breed probabilities) utilize a reference panel of purebred animals' genotypes.
#' The genotypes are phased, and the resulting haplotypes are used to construct haplotype frequences,
#' genotype frequencies given breed $P(g|b)$, and finally breed probabilities given a genotype
#' $P(b | g)$. This function is used to construct a reference panel of breed probabilities.
#' 
#' @param geno a matrix or data.frame with genotypes of reference animals. SNPs in columns and 
#' individuals in rows. Genotypes must be coded as dosage of allele 'b' {0, 1, 2}. 
#' @param map a data.frame containing SNP map information for each SNP present in geno
#' @param ped a data.frame pedigree providing family information for each individual in geno. The first
#' column of the pedigree is for ID, second is for sire/father ID, and third is for dam/mother ID.
#' @param path a character represting the path to the FImpute binary. If omitted, assumes FImpute binary
#' resides along PATH.
#' @param groups a list of character vectors with the names of IDs meant to be processed as groups. 
#' Names of list should be the names of the groups, with each element containing IDs.
#' @import snpTools
#' @import magrittr
#' @export
#' @return a data.frame containing breed probabilities for each genotype made possible from the
#' reference animal haplotypes
build_KBP <- function(geno,
                      map,
                      ped = NULL,
                      path = NULL,
                      groups = NULL) {
  
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
  snpTools::fimpute_run(geno = geno,
                        map = map,
                        ped = ped,
                        path = path,
                        groups = groups,
                        exclude_chr = "1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21")
  message("Obtaining reference haplotypes from FImpute")
  
  breed_list <- list()
  # Obtain haplotypes from each group (each group should correspond to each reference breed)
  for (i in 1:length(groups)) {
    
    # Get SNP names
    snps <- read.table(paste0(names(groups)[i], "_geno_fimpute.txt"),
                       nrow = 1) %>%
              unlist() %>%
                as.character()
    
    # Read haplotype file as a single string
    breed_list[[i]] <- read.table(paste0(names(groups)[i], "_geno_fimpute.txt"),
                                  skip = 1, 
                                  colClasses = "character")
    
    # Convert to a vector with snpTools::string_to_vec()
    breed_list[[i]] <- t(apply(breed_list[[i]], 1, snpTools::string_to_vec))
    
    # Subset KIT SNPs
    colnames(breed_list[[i]]) <- snps
    breed_list[[i]] <- breed_list[[i]][, kit_snps]
    
    # Remove any haplotypes with missing data (marked by the int 5)
    breed_list[[i]] <- breed_list[[i]][apply(breed_list[[i]] != "5", 1, all), ]
    
    # Convert coding to dosage of allele 2
    breed_list[[i]] <- t(apply(breed_list[[i]], 1, function(x) as.numeric(x) - 1))
    breed_list[[i]] <- t(apply(breed_list[[i]], 1, as.character))
    
    # Compute haplotype frequencies and tabulate results
    breed_list[[i]] <- haplo_comp(breed_list[[i]], haplotypes = TRUE)
    breed_list[[i]] <- data.frame("haplotypes" = names(breed_list[[i]]),
                                  "frequency" = breed_list[[i]],
                                  stringsAsFactors = FALSE)
  }
  
  # Merge haplotypes
  hap_data <- breed_list[[1]]
  for (i in 2:length(breed_list)) {
    hap_data <- merge(hap_data, breed_list[[i]], by = "haplotypes", all = TRUE)
  }
  
  colnames(hap_data) <- c("Haplotype", names(groups))
  return(hap_data)
}