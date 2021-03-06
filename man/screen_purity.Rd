% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/screen_purity.R
\name{screen_purity}
\alias{screen_purity}
\title{Calculate both GWBC and KBP for a test animal}
\usage{
screen_purity(geno, gwbc_ref, kbp_ref, mbp_ref = NULL)
}
\arguments{
\item{geno}{matrix containing genotyping data (dosage of allele B) for animals (in rows) for
any number of SNPs (in columns)}

\item{gwbc_ref}{data.frame containing allele frequencies for any number of SNPs (in rows) for
any number of reference populations (in columns). See an example with `data("GWBC_ref_A")`.}

\item{kbp_ref}{data.frame containing genotype probabilities for various multi-locus genotypes
surrounding the KIT region (the output of build_KBP).}

\item{mbp_ref}{data.frame containing genotype probabilities for various multi-locus genotypes
surrounding the MC1R region (the output of build_MBP). (OPTIONAL) If supplied, it will
assume Hampshire animals are being evaluated in \code{geno}.}
}
\value{
data.frame with various GWBC and KBP values for each animal in geno
}
\description{
GWBC (Genome-wide breed composition) and KBP (KIT-based breed probabilities) are calculations
done together to assess the breed composition of a test animal and the probability that it is
purebred in the KIT region of the genome, resposible for influencing color segregation patters.
This function will assemble all GWBC and KBP results for any number of genotyped animals.
}
