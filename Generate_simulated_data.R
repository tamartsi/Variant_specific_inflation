########################################################################################################
# Generating simulated data to demonstrate stratified variance models.
#
# Written by Tamar Sofer
# 2019-11-19
######################################################################################################## 
rm(list = ls())
setwd("~/Documents/GitHub/Variant_specific_inflation/")

require(GWASTools)
require(dummies)
## read genotype data

gds <- GdsGenotypeReader("Genotypes.gds")
scanID <- getScanID(gds)
snpID <- getSnpID(gds)
geno <- getGenotype(gds)
rownames(geno) <- snpID
colnames(geno) <- scanID
close(gds)


## simulate phenotype data	
group1_inds <- 1:100
group2_inds <- 101:350
group3_inds <- 351:500

phen <- data.frame( scanID = scanID,
					group = c(rep("g1", length(group1_inds)), 
							 rep("g2", length(group2_inds)),
							 rep("g3", length(group3_inds))),
					age = rnorm(length(scanID), mean = 40, sd = 5), 
					stringsAsFactors = FALSE)
resid_err_vec <- rnorm(length(scanID), mean = 0, sd = c(rep(1, length(group1_inds)), 
														 rep(2, length(group2_inds)),
														 rep(3, length(group3_inds))))					
phen$trait <- phen$age*1.4 + as.numeric(dummy(phen$group) %*% c(-1, 2, 1.5))  + resid_err_vec

					
save(phen, file = "Phenotypes.RData")