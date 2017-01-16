# Methods to build reference panels for GWBC and KBP, which are loaded upon breedTools attachment
# using data(). A summary of "raw" data objects in data-raw/ that were used to generate reference
# panels:

# LowD_chip_maps.RData - a SNP chip map for the GGP-LD.
# trio_gpData_fix.RData - genotypes from trio animals.
# trio_IDs.RData - a list of animals in trio_gpData_fix that are parents.
# trio_ped_fimpute.txt - a pedigree file for trio animals.
# additional_ref_geno.RData - genotypes from additional animals beyond the Trio set

# First reference panel for GWBC and KBP (GWBC_ref_A and KBP_ref_A, respectively) ------------------
# GWBC_ref_A used parents of original Trio study
# KBP_ref_A used all animals from the original Trio study
# Both GWBC_ref_A and KBP_ref_A were used to generate the TAS manuscript by Funkhouser et. al

# GWBC_ref_A ---------------------------------------------------------------------------------------
# Load 8K SNP chip map
load("data-raw/LowD_chip_maps.RData")

# Load trio genotypes (contained in gp.Trio) and IDs of parents
load("data-raw/trio_gpData_fix.RData")
load("data-raw/trio_IDs.RData")

# Subset parent genotypes, filter and calculate allele frequencies with snpTools
trio_geno_par <- gp.Trio$geno[rownames(gp.Trio$geno) %in% parIDs, ]
trio_geno_par_lowD <- trio_geno_par[, colnames(trio_geno_par) %in% rownames(org_LowD_chip)]
trio_par_ids <- list("Duroc" = as.character(DurocIDs)[as.character(DurocIDs) %in% parIDs],
                     "Hampshire" = as.character(HampshireIDs)[as.character(HampshireIDs) %in% parIDs],
                     "Landrace" = as.character(LandraceIDs)[as.character(LandraceIDs) %in% parIDs],
                     "Yorkshire" = as.character(YorkshireIDs)[as.character(YorkshireIDs) %in% parIDs])
  
GWBC_ref_A <- snpTools::filter_geno(trio_geno_par_lowD) %>%
                breedTools::allele_freq(trio_par_ids)

save(GWBC_ref_A, file = "data/GWBC_ref_A.RData")

# KBP_ref_A ----------------------------------------------------------------------------------------
# Load pedigree information for trios
trio_ped <- read.table("data-raw/trio_ped_fimpute.txt", header = TRUE)

# Assemble trio genotypes, map, and list of ids for build_KBP
trio_geno <- gp.Trio$geno
map_60K <- gp.Trio$map
trio_ids <- list("Duroc" = as.character(DurocIDs),
                 "Hampshire" = as.character(HampshireIDs),
                 "Landrace" = as.character(LandraceIDs),
                 "Yorkshire" = as.character(YorkshireIDs))

KBP_ref_A <- breedTools::build_KBP(geno = trio_geno, 
                                   map = map_60K, 
                                   ped = trio_ped,
                                   path = "~/Programs/bin/",
                                   groups = trio_ids,
                                   parent = TRUE)

save(KBP_ref_A, file = "data/KBP_ref_A.RData")

# Second reference panel for GWBC and KBP (GWBC_ref_B and KBP_ref_B, respectively) -----------------
# GWBC_ref_B used parents of original Trio study, plus all marc animals and a subset of Yorkshire 
# sires (see SF_PG_I/breed_compos/6-update_ref_and_add_sire_const.R for sire selection)
# KBP_ref_B used all animals from the original Trio study plus those added above

# GWBC_ref_B ---------------------------------------------------------------------------------------
# Load additional animals for reference panel
load("data-raw/additional_ref_geno.RData")

# Merge additional animals with trio parents and trim to 8K density
trioPar_marc_sire_geno <- 
  snpTools::merge_geno(durocMarcGenoDose,
                       landraceMarcGenoDose,
                       hampshireMarcGenoDose,
                       yorkshireMarcGenoDose,
                       sires_ref_geno,
                       trio_geno_par)

trioPar_marc_sire_geno_lowD <- 
  trioPar_marc_sire_geno[, colnames(trioPar_marc_sire_geno) %in% rownames(org_LowD_chip)]

# Check if animals have duplicate genotypes (there are duplicate IDs). We don't want to double count
# any animals for the allele frequency calculation
dups <- rownames(trioPar_marc_sire_geno_lowD)[duplicated(rownames(trioPar_marc_sire_geno_lowD))]
sapply(dups, function(x) {
  dup_geno <- trioPar_marc_sire_geno_lowD[rownames(trioPar_marc_sire_geno_lowD) %in% x, ]
  diffs <- dup_geno[1, ] - dup_geno[2, ]
  sum(abs(diffs), na.rm = TRUE)
})
# 471827005      8400      8670     40040   2484769   2485411 308563005 308567002 252222006 390769001 
# 1         0         1         3         0         0         0         0         1         0 
# 445958005 285546006 285725004 287528006 437397015 449710003 452583002 461433004 467888002 468964007 
# 2         2         0         0         0        11         1         1         0         1 
# 468985008 469345001 470569011 470915001      8670 
# 6         0        10         0         1 


# How do these compare to two animals chosen at random?
rand_names <- sample(rownames(trioPar_marc_sire_geno_lowD), 2)
rand_geno <- trioPar_marc_sire_geno_lowD[rownames(trioPar_marc_sire_geno_lowD) %in% rand_names, ]
sum(abs(rand_geno[1, ] - rand_geno[2, ]), na.rm = TRUE)
# [1] 5518 (will be slightly different each time)

# Remove duplicate IDs from genotypes (the same will be needed to generate KBP_ref_B)
trioPar_marc_sire_geno_lowD <- 
  trioPar_marc_sire_geno_lowD[!duplicated(rownames(trioPar_marc_sire_geno_lowD)), ]

# Assemble names of animals in reference panel of each breed
trioPar_marc_sire_names <- 
  list("Duroc" = c(as.character(DurocIDs)[as.character(DurocIDs) %in% parIDs], 
                   rownames(durocMarcGenoDose)),
       "Hampshire" = c(as.character(HampshireIDs)[as.character(HampshireIDs) %in% parIDs],
                       rownames(hampshireMarcGenoDose)),
       "Landrace" = c(as.character(LandraceIDs)[as.character(LandraceIDs) %in% parIDs],
                      rownames(landraceMarcGenoDose)),
       "Yorkshire" = c(as.character(YorkshireIDs)[as.character(YorkshireIDs) %in% parIDs],
                       rownames(yorkshireMarcGenoDose),
                       rownames(sires_ref_geno)))

# Ensure no duplicate names in `trioPar_marc_sire_names`. There should be a total of 1179 names
# after filtering.
trioPar_marc_sire_names <- 
  lapply(trioPar_marc_sire_names, function(x) {
    x[!duplicated(x)]
})

# Filter and calculate allele frequencies
GWBC_ref_B <- snpTools::filter_geno(trioPar_marc_sire_geno_lowD) %>%
                breedTools::allele_freq(trioPar_marc_sire_names)

save(GWBC_ref_B, file = "data/GWBC_ref_B.RData")

# KBP_ref_B ----------------------------------------------------------------------------------------
# Merge additional animals with trios. Animals in this dataset must not be duplicates or have
# duplicate IDs
trio_marc_sire_geno <- 
  snpTools::merge_geno(trio_geno_par,
                       yorkshireMarcGenoDose,
                       landraceMarcGenoDose,
                       hampshireMarcGenoDose,
                       durocMarcGenoDose,
                       sires_ref_geno)

trio_marc_sire_geno <- 
  trio_marc_sire_geno[!duplicated(rownames(trio_marc_sire_geno)), ]

trio_marc_sire_names <- 
  list("Duroc" = c(as.character(DurocIDs)[as.character(DurocIDs) %in% parIDs], 
                   rownames(durocMarcGenoDose)),
       "Hampshire" = c(as.character(HampshireIDs)[as.character(HampshireIDs) %in% parIDs], 
                       rownames(hampshireMarcGenoDose)),
       "Landrace" = c(as.character(LandraceIDs)[as.character(LandraceIDs) %in% parIDs], 
                      rownames(landraceMarcGenoDose)),
       "Yorkshire" = c(as.character(YorkshireIDs)[as.character(YorkshireIDs) %in% parIDs], 
                       rownames(yorkshireMarcGenoDose),
                       rownames(sires_ref_geno)))

KBP_ref_B <- breedTools::build_KBP(geno = trio_marc_sire_geno, 
                                   map = map_60K, 
                                   ped = trio_ped,
                                   path = "~/Programs/bin/",
                                   groups = trio_marc_sire_names,
                                   parent = FALSE,
                                   reference = TRUE)

save(KBP_ref_B, file = "data/KBP_ref_B.RData")
