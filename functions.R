require(zoo);
require(DESeq2);



# =============================================================================================
# The function of generating a state sequence from Markov model
# ---------------------------------------------------------------------------------------------
# W: the length of the sequence
# A: the transition matrix of the Markov model
# s_0: the state of the first position of the sequence
# Seed: random seed
# ---------------------------------------------------------------------------------------------
gamma_generator = function(W, A, s_0, seed) {
  set.seed(seed);
  gamma <- rep(s_0, W);
  if (length(colnames(A)) == 0) {
    state <- 0:(dim(A)[1] - 1);
    colnames(A) <- state;
  } else {
    state <- colnames(A);
  }
  for (w in 2:W) {
    gamma[w] <- sample(state, 1, replace = TRUE, prob = A[, as.character(gamma[w - 1])]);
  }
  return (gamma);
}
# =============================================================================================



# =============================================================================================
# The function of generating count data from hidden Markov model
# ---------------------------------------------------------------------------------------------
# W: the length of the sequence
# design_matrix: the design matrix
# A: the transition matrix of the Markov model
# mu: the mean of the emission
# sigma: the standard deviation of the emission
# Seed: random seed
# ---------------------------------------------------------------------------------------------
data_generator_hmm = function(W, ip, A, mu, pi, family, phi_lwr, phi_upp, seed) {
  set.seed(seed);
  n <- length(ip);
  
  # Generate the hidden states z
  z <- as.numeric(gamma_generator(W, A, 0, seed));

  # Generate the dispersion phi
  phi <- runif(W, phi_lwr, phi_upp)
  
  # Generate the observed count data Y
  Y <- matrix(0, nrow = n, ncol = W);
  H <- matrix(0, nrow = n, ncol = W);
  index <- which(ip == 1);
  for (w in 1:W) {
    if(z[w] == 1) {
      if (family == "poisson" || family == "zip") {
        Y[index, w] <- rpois(length(index), mu[3]);
      } else if (family == "nb" || family == "zinb") {
        Y[index, w] <- rnbinom(length(index), mu = mu[3], size = phi[w]);
      }
    } else {
      if (family == "poisson" || family == "zip") {
        Y[index, w] <- rpois(length(index), mu[1]);
      } else if (family == "nb" || family == "zinb") {
        Y[index, w] <- rnbinom(length(index), mu = mu[1], size = phi[w]);
      }
    }
  }
  index <- which(ip == 0);
  for (w in 1:W) {
    if (family == "poisson" || family == "zip") {
      Y[index, w] <- rpois(length(index), mu[2]);
    } else if (family == "nb" || family == "zinb") {
      Y[index, w] <- rnbinom(length(index), mu = mu[2], size = phi[w]);
    }
  }
  if (family == "zinb" || family == "zip") {
    H <- matrix(rbinom(n*W, 1, pi), nrow = n)
    Y[which(H == 1)] <- 0;
  }

  return(list(Y = Y, ip = ip, A = A, mu = mu, z = z, phi = phi, H = H, family = family));
}
# =============================================================================================



# =============================================================================================
# The function of estimating size factors from sequencing data
# ---------------------------------------------------------------------------------------------
# Y: a n-by-W count matrix, where n is the number of samples and W is the number of features
# s_mode: "total", "median", "quantile"
# ---------------------------------------------------------------------------------------------
size_factors_estimator = function(Y, s_mode) {
  if (s_mode == "total") {
    s <- rowSums(Y)/sum(Y);
  } else if (s_mode == "median") {
    geoMeans = apply(t(Y), 1, gm_mean)
    s <- estimateSizeFactorsForMatrix(t(Y), geoMeans = geoMeans);
    s <- s/sum(s);
  } else if (s_mode == "quantile") {
    s <- apply(Y, 1, quantile, 0.75)
    s <- s/sum(s);
  }
  g <- colSums(Y);
  if (sum(g == 0) > 0) {
    g[which(g == 0)] <- 1;
  }
  return (s);
}



# =============================================================================================



# =============================================================================================
# Miscellaneous
# ---------------------------------------------------------------------------------------------
tabulate_error = function(gamma_true, gamma) {
  table = matrix(0L, 2, 2);
  p <- length(gamma_true);
  for (i in 1:p) {
    table[gamma[i] + 1, gamma_true[i] + 1] <- table[gamma[i] + 1, gamma_true[i] + 1] + 1;
  }
  return (table);
}

# ---------------------------------------------------------------------------------------------
roc = function(gamma_true, mpv, cutoff, increasing) {
  result <- matrix(NA, nrow = length(cutoff), ncol = 2);
  for (i in 1:length(cutoff)) {
    if (increasing) {
      tab <- tabulate_error(gamma_true, mpv >= cutoff[i]);
    } else {
      tab <- tabulate_error(gamma_true, mpv < cutoff[i]);
    }
    result[i, 1] <- tab[2, 1]/sum(tab[, 1]);
    result[i, 2] <- tab[2, 2]/sum(tab[, 2]);
  }
  if (increasing) {
    result <- rbind(c(1, 1), result);
    result <- rbind(result, c(0, 0));
  } else {
    result <- rbind(c(0, 0), result);
    result <- rbind(result, c(1, 1));
  }
  
  auc <- sum(diff(rev(result[, 1]))*rollmean(rev(result[, 2]), 2));
  return (list(result = result, auc = auc));
}

# ---------------------------------------------------------------------------------------------
BayFDR <- function(PPI, alpha){
  PPI_sorted = sort(PPI,decreasing = TRUE)
  k = 1
  fdr = 0
  while(fdr<alpha){
    fdr = mean(1-PPI_sorted[1:k])
    k = k+1
  }
  return(PPI_sorted[k])
}

# ---------------------------------------------------------------------------------------------
gm_mean = function(x, na.rm = TRUE){
  return(exp(sum(log(x), na.rm=na.rm)/length(x)));
}
# =============================================================================================
