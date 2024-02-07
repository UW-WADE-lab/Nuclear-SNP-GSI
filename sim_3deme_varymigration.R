### Preliminary Mixture Simulations
### AVC January 2024

### set up environment ---------------------------------------------------------

library(strataG)
library(tidyverse)

### simulate 3 pops, sample 1 individual, vary migration -----------------------

# set up 3 populations, sample 1 individual per pop
num.demes <- 3

demes <- fscSettingsDemes(fscDeme(15000, 100), fscDeme(15000, 100), fscDeme(10000, 100))

# parameters - to calculate MAF for each population, sample 10% of abundance

# questions - set population sizes based on recent abundance estimates?

# set up SNPs of sequence length 1 on 1000 chromosomes (no linkage), mutation rate 1e-6 
genetics <- fscSettingsGenetics(fscBlock_snp(1, 1e-6), num.chrom = 300)

# parameters - setting 300 SNPs for comparability with NWFSC salmon SNP panel

# questions - get mutation rate estimate from reference panel from each population? or est for salmon?

# set up mutation rate variability
m.vec <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
mig.rate <- m.vec / (num.demes - 1)

# question - what migration rates would be realistic to test for Puget Sound salmon?

zero.mig <- matrix(data = 1, nrow = num.demes, ncol = num.demes)
diag(zero.mig) <- 0
  
p.list <- lapply(1:length(mig.rate), function(i) {

# write settings file and run simulation
p <- fscWrite(
    demes = demes,
    migration = fscSettingsMigration(matrix(rep(mig.rate[i], num.demes ^ 2), nrow = num.demes) * zero.mig),
    genetics = genetics, 
    label = "sim3deme.varymig",
    use.wd = TRUE
  )
  
p <- fscRun(p, all.sites = F, inf.sites = FALSE, exec = "fsc28")

})

#### Use simulated data to estimate maf and generate random mixture samples ----

# read the genetic data
  
arp.list <- lapply(1:length(p.list), function(i) {
  
  arp <- fscReadArp(p.list[[i]])
  
})

# generate allele frequencies at each locus for all values of migration

maf <- bind_rows(arp.list, .id = "runID") %>% 
  pivot_longer(-c(runID, id, deme), names_to = "locus", values_to= "allele") %>% 
  separate(locus, sep = "\\.", into = c("locus","position")) %>% 
  group_by(runID, deme, locus) %>% 
  count(allele) %>% 
  mutate(prop.allele = n/sum(n))

# generate mixture samples
# 10 per migration value

mixture.sample <- data.frame()
mixture.samples_all <- data.frame()

for (k in 1:length(arp.list)){
  for (i in 1:10) {
    mixture.sample.i <- arp.list[[k]] %>% 
      sample_n(2, replace = FALSE) %>% 
      pivot_longer(-c(id, deme), names_to = "locus", values_to= "allele") %>% 
      separate(locus, sep = "\\.", into = c("locus","position")) %>% 
      count(locus,allele) %>% 
      # mutate(prop.allele = n/sum(n)) %>% 
      mutate(sampleID = i)
    
    mixture.sample <- bind_rows(mixture.sample, mixture.sample.i)
  }
  
  runID <- k
  
  mixture.sample <- mixture.sample %>% 
    mutate(runID = k)
  
  mixture.samples_all <- bind_rows(mixture.samples_all, mixture.sample)
  
}

#what input format do the two algorithms take? do i need pivot wider?
mixture.samples_wide <- mixture.samples_all %>% 
  unite(col = "locus", locus:allele, sep = "_") %>% 
  pivot_wider(id_cols = c(sampleID, runID), names_from = "locus", values_from = "n", values_fill = 0)
