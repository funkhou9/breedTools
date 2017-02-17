#' Build reference panel for MBP (MC1R-based breed probability) calculation
#' 
#' MBP (MC1R-based breed probabilities) utilize a reference panel of purebred animals' genotypes.
#' The genotypes are phased, and the resulting haplotypes are used to construct haplotype frequences,
#' genotype frequencies given breed \eqn{P(g | b)}, and finally breed probabilities given a genotype
#' \eqn{P(b | g)}. This method is used to construct a reference panel of breed probabilities for which
#' one can use to test the probability that an animal is of a certain breed, in the region of MC1R.
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
#' @param keep_fimpute logical if TRUE, output files generated from snpTools::fimpute_run are kept. 
#' By default, build_KBP will use fimpute_run to build haplotypes, but as soon as haplotypes are 
#' extracted, fimpute output will be deleted so not to clutter hard drive space.
#' @param parent logical if TRUE, breed probabilities or haplotypes will be returned only from
#' parents and not from progeny,
#' @param reference logical use TRUE if you want `groups` to be used to specify the reference animals
#' used by FImpute. Not compatable with parent TRUE.
#' @import snpTools
#' @import magrittr
#' @export
#' @return a data.frame containing breed probabilities for each genotype made possible from the
#' reference animal haplotypes
build_MBP <- function(geno,
                      map,
                      ped = NULL,
                      path = NULL,
                      groups = NULL,
                      keep_fimpute = FALSE,
                      parent = FALSE,
                      reference = FALSE) {
  
  build_local_reference(chr = "6",
                        markers = mc1r_snps,
                        geno_ = geno,
                        map_ = map,
                        ped_ = ped,
                        path_ = path,
                        groups_ = groups,
                        keep_fimpute_ = keep_fimpute,
                        parent_ = parent,
                        reference_ = reference)
}



