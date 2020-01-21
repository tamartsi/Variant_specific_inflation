########################################################################################################
# functions for
# Computing approximate inflation factors based on allele frequencies and trait standard deviations in groups
# Generating a figure QQ-plots by categories of variants, according to expected inflation factors
# 		under the homoengeoue variances model.
#
# Written by Tamar Sofer and Kenneth M Rice
# 2019-11-19
######################################################################################################## 

require(tidyr)
require(data.table)
require(ggplot2)

# A function that computes the naive and the true test standard error based on 
#                    sample sizes, frequencies, and phenotype standard deviations in multiple groups.
#         Uses a "sandwitch" formula under heterogeneous variance model
#         Uses standardformulat under homogeneous residual variance models.
#     To avoid constructions of large matrices, proportions are used in  lieu of repeated rows.
# eafs: vector of effect allele frequencies in the k studies
# ns: sample sizes in each of the k studies
# sigmas: Residual standard deviations in each of the k studies
# beta0: intercept, default is zeo.
# study_main_effects: length k-1 vector of study main effects; the first study is chosen to be the reference
#     (This can usually just be a vector of zeros)

compute_variant_inflation_factor <- function(eafs, 
                                             ns, 
                                             sigmas, 
                                             beta0=0,
                                             study_main_effects = rep(0, length(sigmas)-1)){
   k <- length(eafs)

   if( (length(study_main_effects) !=(k-1) ) | (length(ns) !=k) | (length(sigmas) !=k) ) stop("length mismatch")

  ## compute sample proportions
  sample_proportion <- ns/sum(ns)

  xmat <- data.frame(
                    intercept <- rep(1, 3*k),
                    geno      <- rep(0:2, k),
	                  Z         <-  model.matrix( ~factor(rep(1:k, each=3)))[,-1]	
                    )

  beta <- c(beta0,0,study_main_effects)
  pr.z <- sample_proportion[rep(1:k, each=3)]

  xmat <- as.matrix(xmat)
  pr.x <- dbinom(xmat[,2], 2, rep(eafs, each=3))

  props <- pr.x*pr.z

  Bmat <- t(xmat) %*% diag(props) %*% xmat
  Amat <- t(xmat) %*% diag(props) %*% diag( sigmas[rep(1:k, each=3)] ) %*% xmat

  # sandwitch formula for computing the SE of effect size allowing for heterogeneous vairances:
  SE_true <- sqrt( (solve(Bmat) %*% Amat %*% solve(Bmat))[2,2] )
   
  # standard SE formula:
  SE_naive <- sqrt( solve(Bmat)[2,2] * sum(props*(sigmas[rep(1:k, each=3)]) ) )

  c(SE_true=SE_true, SE_naive=SE_naive, Inflation_factor = SE_true/SE_naive)
}


# pval_df is a data.frame with each column corresponding to a different analysis. rows correspond to variants.
#   it can have only a single column.
# expected_inf is a vector with the expected inflation per variant. variants correspond to those in pval_df
# inflated_val is the minimum values for which a variant expected_inf is categorized under "expected inflation". 
# deflated_val is the maximum values for which a variant expected_inf is categorized under "expected inflation". 
# values between deflated_val and inflated_val will be plotted in the "about right" region of the figure. 
qq_plot_by_region <- function(pval_df, expected_inf,inflated_val = 1.03, deflated_val= 0.97, 
								thin_high_pvals = 1e5, figure_file_name = "qq_plots_by_region.pdf"){
	stopifnot(nrow(pval_df) == length(expected_inf))
	
	n_analyses <- ncol(pval_df)
	analyses_names <- colnames(pval_df)
	
	pval_df$expected_inf <- expected_inf
	pval_df$inflation_category <- "About right"
	pval_df$inflation_category[pval_df$expected_inf <= deflated_val] <- "Deflated"
	pval_df$inflation_category[pval_df$expected_inf >= inflated_val] <- "Inflated"
	
	# structure data for figure:
	if (n_analyses == 1){
		dat <- data.frame(expected_inf = pval_df$expected_inf,
						  inflation_category = pval_df$inflation_category,
						  analysis = analyses_names,
						  pvalue = pval_df[,analyses_names])
	} else{
		dat <- gather(pval_df, analysis, pvalue, analyses_names)
		}	
	
	dat <- data.table(dat)
	
	### prepare info for "all variants" panel:
	dat_a <- .prepare_values_for_qqplot(dat, thin_high_pvals)
	dat_a$inflation_category <- "All variants"
		
	# prepare info for "deflated" panel:
	dat_d <- .prepare_values_for_qqplot(dat[inflation_category  == "Deflated"], thin_high_pvals)
	
	# prepare info for "inflated" panel:
	dat_i <- .prepare_values_for_qqplot(dat[inflation_category  == "Inflated"], thin_high_pvals)
	
	# prepare info for "About right" panel:
	dat_r <- .prepare_values_for_qqplot(dat[inflation_category  == "About right"], thin_high_pvals)
	
	dat_plot <- rbind(dat_a, dat_d, dat_i, dat_r)
	
	p <- ggplot(dat_plot, 
				aes(x=qnorm, y = qval)) + 
				geom_point(aes(colour = analysis), alpha=0.6) + 
				geom_abline(intercept = 0, slope = 1) + 
				labs(x = "Expected p-value", y = "Observed p-value")  
				
	ggsave(filename = figure_file_name, plot = p + facet_wrap(~ inflation_category, ncol = 2))
			
}


.prepare_values_for_qqplot <- function(dat, thin_high_pvals){
	dat[,qval :=  -log(pvalue, 10), by = "analysis"]
	dat[,qnorm := -log(1:length(pvalue)/length(pvalue), 10)[rank(pvalue)], by = "analysis"]
	
	if (!is.null(thin_high_pvals)){
		inds.1 <- which(dat$qval < 1)
		inds.2 <- which(dat$qval > 1 & dat$qval < 2)
		inds.3 <- c(which(dat$qval == 1), which(dat$qval >= 2))
		### sample 1000 points from each 
		dat <- dat[c(sample(inds.1, min(thin_high_pvals, length(inds.1))), 
							sample(inds.2, min(thin_high_pvals, length(inds.2))), 
							inds.3),]
	}

	return(dat)
	
}




