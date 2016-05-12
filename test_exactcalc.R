# Run this just to test that the c++ implementation is giving the correct answer
# Only use with p<=5!
# Run from the command line with Rscript test_exactcalc.R BOUNDS_FILE COR_FILE


library(combinat)
library(mvtnorm)



args <- commandArgs(trailingOnly=TRUE)

bounds_name <- args[1]
cor_name <- args[2]
bounds <- unlist(read.table(bounds_name))
cor_vec <- unlist(read.table(cor_name))


# Make the matrix from a vector
p <- length(bounds)
sig_mat <- matrix(data=1, nrow=p, ncol=p)
for (i in 2:p) {
	for (j in 1:(i-1)) {
		sig_mat[i,j] <- cor_vec[j + (i-2)*(i-1)/2] 
		sig_mat[j,i] <- sig_mat[i,j]
	}
}

# Make the delta_p matrix to extend integration dimension from p to 2p-1
Delta_p <- diag(rep(1,p))
for(i in 1:(p-1)) {
	new_row <- c(rep(0, i-1), -1, 1, rep(0, p-1-i))
	Delta_p <- rbind(Delta_p, new_row)
}



############################################################
############################################################
# Calculates the noncrossing probability for a specific set of bounds, given
# the correlation structure.
# Pass in \Delta_p matrix because it's unwieldly to generate each time
# Call it separately for each test (i.e. three times for HC, GHC, BJ)
# Outer loop is over p! permutations of a, inner loop is over 2^p permutations of s
# Returns the p-value and total pmvnorm error

calc_pvalue <- function(bounds, sig_mat, Delta_p) {
	
	p <- nrow(sig_mat)
	
	# Holds the sum of each loop over class_A (order permutation) as well as error sum
	loop_sums <- matrix(data=NA, nrow=(factorial(p)), ncol=2)
	
	class_A <- permn(1:p)			# a list of vectors
	class_S <- expand.grid(rep(list(0:1), p))			# a matrix
	class_S[class_S==0] <- -1							# each row is one permutation of +-
	lower_bounds <- c(0, rep(-Inf, (p-1)), rep(0, (p-1)))
	upper_bounds <- c(bounds, rep(Inf, (p-1)))
	
	# Outer loop over all permutations of a, i.e. shuffling the order of observations
	for(i in 1:(factorial(p))) {
		
		# P_a is how we permute the variance matrix to account for the ordering of the Y
		a <- unlist(class_A[i])
		I_p <- diag(x=1, nrow=p)
		P_a <- I_p[a,]
		
		# Permute the variance matrix to get \Sigma^(a)
		sig_mat_a <- P_a %*% sig_mat %*% t(P_a)
		
		# Inner loop over all permutations of s, to get the multivariate folded normal pdf
		temp_sum_prob <- 0
		temp_sum_err <- 0
		for(j in 1:(2^p)) {
			
			# Multiply our variance matrix by s to account for +- transformation to Y
			s <- class_S[j,]
			Lambda_s <- diag(s, nrow=p)
			sig_mat_a_s <- Lambda_s %*% sig_mat_a %*% t(Lambda_s)
			
			# Multiply the new variance matrix by \Delta_p to get the variance for T
			# T is the final multivariate vector we are interested in, has length 2p-1
			T_var <- Delta_p %*% sig_mat_a_s %*% t(Delta_p)
			
			# Calculate the multiple integral, a multivariate normal cdf
			integral_result <- pmvnorm(lower=lower_bounds, upper=upper_bounds, sigma=T_var)
			temp_sum_prob <- temp_sum_prob + integral_result[1]
			temp_sum_err <- temp_sum_err + attributes(integral_result)$err			
		}
		
		# Record our results in the loop_sums matrix
		loop_sums[i,1] <- temp_sum_prob
		loop_sums[i,2] <- temp_sum_err
	}

	# The p-value is 1-minus the probability that we stay within the bounds
	pvalue <- 1 - sum(loop_sums[,1])
	err_sum <- sum(loop_sums[,2])
		
	return(list(pvalue=pvalue, err_sum=err_sum))
}

results <- calc_pvalue(bounds=bounds, sig_mat=sig_mat, Delta_p=Delta_p)
print(results)
