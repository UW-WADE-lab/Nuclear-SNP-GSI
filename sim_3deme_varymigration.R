### Preliminary Mixture Simulations
### AVC January 2024

### set up environment ---------------------------------------------------------

library(strataG)
library(tidyverse)

### simulate 3 pops, sample 1 individual, vary migration -----------------------

# set up 3 populations, sample 1 individual per pop
num.demes <- 3

demes <- fscSettingsDemes(fscDeme(15000, 1), fscDeme(15000, 1), fscDeme(10000, 1))

# questions - set population sizes based on recent abundance estimates?

# set up SNPs of sequence length 1 on 1000 chromosomes (no linkage), mutation rate 1e-6 
genetics <- fscSettingsGenetics(fscBlock_snp(1, 1e-6), num.chrom = 1000)

# questions - should sequence.length be longer per recommendation from Excoffier?
#           - set number of chromosomes to the number of SNPs in the NWFSC panel?
#           - get mutation rate estimate from reference panel from each population? or est for salmon?

# set up mutation rate variability
m.vec <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
mig.rate <- m.vec / (num.demes - 1)

# question - what migration rates would be realistic to test for Puget Sound salmon?

zero.mig <- matrix(data = 1, nrow = num.demes, ncol = num.demes)
diag(zero.mig) <- 0
  
migration <- fscSettingsMigration(matrix(rep(mig.rate[1], num.demes ^ 2), nrow = num.demes) * zero.mig,
                       matrix(rep(mig.rate[2], num.demes ^ 2), nrow = num.demes) * zero.mig,
                       matrix(rep(mig.rate[3], num.demes ^ 2), nrow = num.demes) * zero.mig,
                       matrix(rep(mig.rate[4], num.demes ^ 2), nrow = num.demes) * zero.mig,
                       matrix(rep(mig.rate[5], num.demes ^ 2), nrow = num.demes) * zero.mig,
                       matrix(rep(mig.rate[6], num.demes ^ 2), nrow = num.demes) * zero.mig)
  
# write settings file and run simulation
p <- fscWrite(
    demes = demes,
    migration = migration,
    genetics = genetics, 
    label = "sim3deme.varymig",
    use.wd = TRUE
  )
  
p <- fscRun(p, all.sites = F, inf.sites = T, exec = "fsc28")

# question - is there a way to set this up so that it loops through all of the matrices in the migration 
# settings object? 
  
# read the genetic data
  
arp.file <- fscReadArp(p)

maf <- arp.file %>% 
  pivot_longer(-c(id, deme), names_to = "locus", values_to= "allele") %>% 
  separate(locus, sep = "\\.", into = c("locus","position")) %>% 
  group_by(locus) %>% 
  count(allele) %>% 
  mutate(prop.allele = n/sum(n))


