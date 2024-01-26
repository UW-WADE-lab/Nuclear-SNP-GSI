### First test run fastsimcoal
### AVC January 2024

### set up environment ---------------------------------------------------------

library(strataG)
library(tidyverse)

### strataG fastsimcoal example ------------------------------------------------

#define deme(s)

deme0 <- fscDeme(deme.size = 1000, sample.size = 4, sample.time = 0, inbreeding = 0, growth = 0)

#supply demes to fastsimcoal settings

demes <- fscSettingsDemes(deme0, ploidy = 2)

#specify microsatellites to simulate

msats <- fscBlock_microsat(num.loci = 1, mut.rate = 1e-3)

#create input for genetic data

genetics <- fscSettingsGenetics(msats, num.chrom = 5) #creates 5 unlinked loci

#write fastsimcoal parameters file

ex1.params <- fscWrite(demes = demes, genetics = genetics, label = "ex1", use.wd = TRUE)

#run the simulation

ex1.params <- fscRun(ex1.params, num.sim = 1, exec = "fsc28")

#read Arlequin output file

arp.file <- fscReadArp(ex1.params)
str(arp.file)
head(arp.file)

### SNPs exercise --------------------------------------------------------------

rm(list = ls())


#simulate SNPs

# deme with 1000 individuals and a sample size of 10 individuals
demes <- fscSettingsDemes(fscDeme(deme.size = 1000, sample.size = 10))

# SNP loci 10bp long, mutation rate 1e-6, 1000 total chromosomes
genetics <- fscSettingsGenetics(fscBlock_snp(sequence.length = 10, mut.rate = 1e-6), num.chrom = 1000)

# write fastsimcoal settings file
p <- fscWrite(demes = demes, genetics = genetics, label = "ex2.snps.1k", use.wd = TRUE)

# run simulation
p <- fscRun(p, all.sites = F, exec = "fsc28")

# read arp file

snp.df <- fscReadArp(p)
head(snp.df)
snp.df[1:6,1:6]

# count SNPs per chromosome

snp.occ.freq <- snp.df %>% 
  pivot_longer(-c(id,deme), names_to = "snp", values_to = "allele") %>% 
  separate(snp, into = c("chr.locus", "allele.pos", "L"), sep = "_") %>% 
  group_by(id) %>% 
  count(chr.locus, name = "num.snp.chr") %>% 
  ungroup() %>% 
  distinct(chr.locus, .keep_all = TRUE) %>% 
  count(num)

### Now let's look at migration -----------------------------------------------

# input for migration is a matrix showing migration beteween pairs of populations
# 2 - population model:

m <- 0.00001 # this translates to 1 migrant per generation
mig.mat <- matrix(c(0, m, m, 0), nrow = 2)
mig.mat

# SNP simulation where both populations have 1000 individuals and mutation rate 1e-6:

demes <- fscSettingsDemes(fscDeme(1000, 10), fscDeme(1000, 10))
genetics <- fscSettingsGenetics(fscBlock_snp(1, 1e-6), num.chrom = 1000)
p <- fscWrite(
  demes = demes,
  migration = fscSettingsMigration(mig.mat),
  genetics = genetics,
  label = "ex3.mig.ex",
  use.wd=TRUE
)
p <- fscRun(p, all.sites = F, exec = "fsc28")

# read arp file
snp.df <- fscReadArp(p, one.col = F)

# convert to gtypes object to calculate Fst
snp.g <- df2gtypes(snp.df, ploidy = 2)

# fst
overallTest(snp.g, nrep = 0, quietly = T)$result["wcFst", "estimate"]

# What happens as migration increases

# vector of migration rates
m.vec <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)

# simulate across all migration rates
fst <- sapply(m.vec, function(m) {
  mig.mat <- matrix(c(0, m, m, 0), nrow = 2)
  p <- fscWrite(
    demes = demes,
    migration = fscSettingsMigration(mig.mat),
    genetics = genetics, 
    label = "mig.test",
    use.wd = TRUE
  )
  p <- fscRun(p, all.sites = F, inf.sites = T, exec = "fsc28")
  snp.df <- fscReadArp(p, one.col = F)
  snp.g <- df2gtypes(snp.df, ploidy = 2)
  overallTest(snp.g, nrep = 0, quietly = T)$result["wcFst", "estimate"]
})

# bind results together

cbind(
  m = m.vec, 
  Nm = 1000 * m.vec, 
  expected.Fst = 1 / ((4 * 1000 * m.vec) + 1), 
  observed.Fst = fst
)

# the strataG fastsimcoal vignette has code and examples for multiple types of migration:
# island model, stepping stone model, and spatially-explicit model

### Historical Events ----------------------------------------------------------

# Multiple historical events can be defined
# Refer to strataG vignette and fastsimcoal manual for more information

events <- fscSettingsEvents(
  fscEvent(
    event.time = 2000, 
    source = 1, 
    sink = 2, 
    prop.migrants = 0.05, 
    new.size = 1,
    new.growth = 0,
    migr.mat = 0
  ),
  fscEvent(
    event.time = 2980, 
    source = 1, 
    sink = 1, 
    prop.migrants = 0, 
    new.size = 0.04
  ),
  fscEvent(3000, 1, 0, 1, 1),
  fscEvent(15000, 0, 2, 1, 3)
)

### Predefined values ----------------------------------------------------------

# some starter code for looping through multiple values

param.df <- data.frame(
  N = 10 ^ runif(5, 2, 4),
  MIG = 10 ^ runif(5, -8, -5)
)

p <- fscWrite(
  demes = fscSettingsDemes(fscDeme("N", 5), fscDeme("N", 5)),
  migration = fscSettingsMigration(matrix(c(0, "MIG", "MIG", 0), nrow = 2)),
  genetics = fscSettingsGenetics(fscBlock_snp(100, 1e-6), num.chrom = 1000),
  def = fscSettingsDef(param.df),
  label = "param.sim"
)
p <- fscRun(p)

# from what I can tell, this does not loop through all combinations, instead
# it steps through each position in all vectors at once.
