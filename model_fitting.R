# ***README***
# The following script is used to fit MeRIP-Seq data to the Bayesian negative binomial 
# hidden Markov model proposed in the submitted manuscript titled "A Bayesian 
# Hierarchical Model for Analyzing Methylated RNA Immunoprecipitation Sequencing Data."

# Before running the following code, please first load MeRIP-Seq data. The necessary 
# inputs should be 
# (1) a n-by-W count matrix Y, where n is the number of samples and W is the number 
# of bins
# (2) a n-dimensional logical vector ip, indiating the allocation of IP samples.

# Note that the notations in the code and data follow the notations in the manuscript.
# ***END***

# Load functions
source('functions.R');
Rcpp::sourceCpp('zinb_hmm_x2.cpp');



# ========================================================================================
# ========================================================================================
# Load data
# ========================================================================================
# ========================================================================================
Y <- read.table("real_data_example.txt", header = TRUE);
ip <- grepl("ip", names(Y));
Y <- as.matrix(t(Y));
n <- dim(Y)[1];
W <- dim(Y)[2];


# ========================================================================================
# ========================================================================================
# Load algorithm setting
# ========================================================================================
# ========================================================================================
iter <- 10*W;
# choose any one of the three size factor estimators: "total", "median", "quantile"
s_mode <- "total";   


# ========================================================================================
# ========================================================================================
# Load hyperparameters, please adjust the other hyperparameters in "zinb_hmm_x2.cpp"
# ========================================================================================
# ========================================================================================
upp_eta_1 <- 1;   # upp_eta_1 = lwr_eta_2


# ========================================================================================
# ========================================================================================
# Load initial configuration
# ========================================================================================
# ========================================================================================
z_start <- rbinom(W, 1, 0.5);



# ========================================================================================
# ========================================================================================
# Implement MCMC algorithm
# ========================================================================================
# ========================================================================================
s <- size_factors_estimator(Y, s_mode);
g <- colSums(Y);
start_time <- proc.time();
M <- zinb_hmm(Y, ip, s, g, upp_eta_1, z_start, iter);
end_time <- proc.time();
time <- end_time - start_time;
# The MCMC outputs are stored in M
# $d_0:        the mean of MCMC samples of d_0
# $delta:         the mean of MCMC samples of delta
# $phi:           the mean of MCMC samples of phi
# $pi:            the mean of MCMC samples of pi
# $H:             the marginal posterior probability of inclusion (PPI) of H
# $accept_d_0: the acceptance rate for updating d_0
# $accept_delta:  the acceptance rate for updating delta
# $accept_phi:    the acceptance rate for updating phi
# $A:             the mean of MCMC samples of A
# $mu:            the mean of MCMC samples of mu
# $sigma2:        the mean of MCMC samples of sigma2
# $z:             the marginal posterior probability of inclusion (PPI) of z
# $z_sum:         the total number of selected methylated bins



# ========================================================================================
# ========================================================================================
# Summarize result
# ========================================================================================
# ========================================================================================
print(paste0("Runtime = ", round(time[3], 1), "s"));
# print(paste0("Acceptance rate (d_0): ", round(M$accept_d_0, 3)));
# print(paste0("Acceptance rate (delta): ", round(M$accept_delta, 3)));
# print(paste0("Acceptance rate (phi): ", round(M$accept_phi, 3)));

plot(M$z_sum, type = "l", xlab = "Iterations", ylab = "Number of selected methylated bins");
plot(M$z[2, ], type = "h", xlab = "Bins", ylab = "Posterior probabilities of inclusion");


