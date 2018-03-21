This repository is used for accessing the performance of the Bayesian negative binomial hidden Markov model, which was proposed in the submitted manuscript titled "A Bayesian hierarchical model for analyzing methylated RNA immunoprecipitation sequencing data." Before running the code, please install R packages Rcpp, DESeq2, and zoo.

When running "model_fitting.R" to fit the proposed model, please first load MeRIP-Seq data data. You can load your own data or one of the two examples in the repository. The necessary inputs should be a n-by-W count matrix denoted by Y, where n is the number of samples and W is the number of bins, and a n-dimensional logical vector ip, indiating the allocation of IP samples.

*Note 1, the notations in the code follow the notations in the manuscript.

*Note 2, the results obtained by running the code in this repository may not be exactly the same as the results reported in the manuscript, because we reported the results by pooling multiple MCMC chains in the manuscript.
