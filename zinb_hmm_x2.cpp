#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double norm_rs(double a, double b);
double half_norm_rs(double a, double b);
double unif_rs(double a, double b);
double exp_rs(double a, double b);
double rnorm_trunc(double mu, double sigma, double lower, double upper);

// [[Rcpp::export]]
Rcpp::List zinb_hmm(IntegerMatrix Y, NumericVector ip, NumericVector s, NumericVector g, double fold, NumericVector z_s, int iter) {
  int n = Y.nrow();
  int W = Y.ncol();
  int Q = 2;
  int total_count_lwr = 10;
  bool store = TRUE;
  bool HMM = TRUE;
  double fold_max = 4.0;
  
  // Hyperparameters
  double a_pi = 1;
  double b_pi = 1;
  double a_d = 0.001;
  double b_d = 0.001;
  double a_phi = 0.001;
  double b_phi = 0.001; 
  double mu_0 = 0;
  double tau_0 = 10;
  double a_0 = 2;
  double b_0 = 1;
  NumericMatrix mu_limit(2, Q);
  mu_limit(0, 0) = 0;
  mu_limit(1, 0) = fold;
  mu_limit(0, 1) = fold;
  mu_limit(1, 1) = 1000;
  
  // Algorithm settings
  int burn = iter/2;
  double tau_phi = 1;
  double tau_d_0 = 0.1;
  double tau_delta = 1;
  LogicalVector flag(W);
  LogicalVector flag_2(W);
  int z_2 = 0;
  
  // Model parameters
  IntegerMatrix H(n, W);
  double pi;
  NumericVector d_0(W);
  NumericVector delta(W);
  NumericVector phi(W);
  arma::vec z(W);
  NumericMatrix A(Q, Q);
  NumericVector mu(Q);
  NumericVector sigma2(Q);
  
  // MCMC samples
  NumericMatrix d_0_store(iter, W);
  NumericMatrix delta_store(iter, W);
  NumericMatrix phi_store(iter, W);
  NumericVector pi_store(iter);
  NumericMatrix mu_store(iter, Q);
  NumericMatrix sigma2_store(iter, Q);
  NumericMatrix A_store(iter, Q*Q);
  
  NumericVector d_0_mean(W);
  NumericVector delta_mean(W);
  NumericVector phi_mean(W);
  double pi_mean = 0;
  NumericVector mu_mean(Q);
  NumericVector sigma2_mean(Q);
  NumericMatrix H_ppi(n, W);
  NumericMatrix A_mean(Q, Q);
  NumericMatrix z_ppi(Q, W);
  IntegerVector z_sum(iter);
  
  // Acceptance rate
  double accept_d_0 = 0;
  double accept_delta = 0;
  double accept_phi = 0;
  double try_d_0 = 0;
  double try_delta = 0;
  double try_phi = 0;
  
  // Temp variables
  int ii, i, w, q, qq, n1, w_temp;
  int count = 0;
  int countcount = 0;
  double hastings, max_temp, sum_temp, phi_temp, d_0_temp, delta_temp, temp, count_temp;
  NumericVector prob_temp_2(2);
  IntegerMatrix C(Q, Q);
  IntegerVector state(Q);
  NumericVector prob_temp(Q);
  NumericVector count_temp_2(Q);
  
  // HMM Initialization
  for(q = 0; q < Q; q++)
  {
    for(qq = 0; qq < Q; qq++)
    {
      C(q, qq) = 0;
      A(q, qq) = 0.5;
    }
  }
  for(w = 0; w < W; w++)
  {
    z(w) = z_s(w);
    count = 0;
    for(i = 0; i < n; i++)
    {
      count = count + Y(i, w);
    }
    if(count <= total_count_lwr)
    {
      flag(w) = 1; 
      z(w) = 0;
    }
    else
    {
      flag(w) = 0;
    }
    if(count == 0)
    {
      flag_2(w) = 1;
    }
    else 
    {
      flag_2(w) = 0;
    }
    if(w > 0 && flag(w) == 0)
    {
      C(z(w), z(w - 1))++;
    }
  }
  temp = Q;
  for(q = 0; q < Q; q++)
  {
    mu(q) = runif(1, mu_limit(0, q), mu_limit(1, q))(0);
    sigma2(q) = 1;
    sigma2_mean(q) = 0;
    mu_mean(q) = 0;
    state(q) = q;
    prob_temp(q) = 1/temp;
    for(w = 0; w < W; w++)
    {
      z_ppi(q, w) = 0;
    }
  }
  for(w = 0; w < W; w++)
  {
    z_2 = z_2 + z(w);
  }
  
  // ZINB Initialization
  n1 = 0;
  for(w = 0; w < W; w++)
  {
    d_0_mean(w) = 0;
    delta_mean(w) = 0;
    phi_mean(w) = 0;
    if(flag_2(w) == 1)
    {
      phi(w) = 0;
      d_0(w) = 0;
      delta(w) = 0;
      for(i = 0; i < n; i++)
      {
        H_ppi(i, w) = 0;
        H(i, w) = 1;
      }
      n1 = n1 + n;
    }
    else
    {
      if(flag(w) == 1)
      {
        delta(w) = 0;
      }
      else
      {
        delta(w) = runif(1, 0, fold_max)(0);
      }
      d_0(w) = rgamma(1, 1, 1)(0);
      phi(w) = 100;
      for(i = 0; i < n; i++)
      {
        H_ppi(i, w) = 0;
        if(Y(i, w) == 0) 
        {
          H(i, w) = rbinom(1, 1, 0.5)(0);
          if (H(i, w) == 1) {
            n1 = n1 + 1;
          }
        } 
        else 
        {
          H(i, w) = 0;
        }
      }
    }
  }
  
  
  
  // MCMC
  count = 0;
  for(ii = 0; ii < iter; ii++)
  {
    // Update pi
    pi = rbeta(1, a_pi + n1, b_pi + n*W - n1)(0);
    
    // Update H
    n1 = 0;
    for(w = 0; w < W; w++)
    {
      if(flag_2(w) == 0)
      {
        for(i = 0; i < n; i++)
        {
          if(Y(i, w) == 0)
          {
            prob_temp_2(0) = phi(w)*(log(phi(w)) - log(s(i)*g(w)*d_0(w) + phi(w))) + log(1 - pi);
            prob_temp_2(1) = log(pi);
            max_temp = max(prob_temp_2);
            prob_temp_2(0) = prob_temp_2(0) - max_temp;
            prob_temp_2(1) = prob_temp_2(1) - max_temp;
            prob_temp_2(0) = exp(prob_temp_2(0));
            prob_temp_2(1) = exp(prob_temp_2(1));
            sum_temp = prob_temp_2(0) + prob_temp_2(1);
            prob_temp_2(0) = prob_temp_2(0)/sum_temp;
            prob_temp_2(1) = prob_temp_2(1)/sum_temp;
            H(i, w) = rbinom(1, 1, prob_temp_2(1))(0);
            n1 = n1 + H(i, w);
          }
        }
      }
      else
      {
        n1 = n1 + n;
      }
    }
    
    // Update phi
    for(w = 0; w < W; w++)
    {
      if(flag_2(w) == 0)
      {
        countcount = 0;
        do {
          phi_temp = rgamma(1, phi(w)*phi(w)/tau_phi, tau_phi/phi(w))(0);
          countcount++;
        } while (phi_temp < 1 && countcount < 1000);
        if (countcount == 1000)
        {
          // phi_temp = runif(1, 1, 100)(0);
          phi_temp = rgamma(1, 100, 1)(0);
        }
        hastings = 0;
        try_phi++;
        for(i = 0; i < n; i++)
        {
          if(H(i, w) == 0) {
            hastings = hastings + phi_temp*log(phi_temp) - lgamma(phi_temp) + lgamma(phi_temp + Y(i, w)) - (phi_temp + Y(i, w))*log(phi_temp + s(i)*g(w)*d_0(w));
            hastings = hastings - (phi(w)*log(phi(w)) - lgamma(phi(w)) + lgamma(phi(w) + Y(i, w)) - (phi(w) + Y(i, w))*log(phi(w) + s(i)*g(w)*d_0(w)));
          }
        }
        hastings = hastings + (a_phi - 1)*log(phi_temp) - b_phi*phi_temp;
        hastings = hastings - ((a_phi - 1)*log(phi(w)) - b_phi*phi(w));
        if(hastings >= log(double(rand()%10001)/10000))
        {
          phi(w) = phi_temp;
          if (ii >= burn) {
            accept_phi++;
          }
        }
      }
    }
    
    // Update d_0
    for(w = 0; w < W; w++)
    {
      if(flag_2(w) == 0)
      {
        d_0_temp = rgamma(1, d_0(w)*d_0(w)/tau_d_0, tau_d_0/d_0(w))(0);
        hastings = 0;
        try_d_0++;
        for(i = 0; i < n; i++)
        {
          if(H(i, w) == 0) {
            if(ip(i) == 0)
            {
              hastings = hastings + Y(i, w)*log(s(i)*g(w)*d_0_temp) - (phi(w) + Y(i, w))*log(phi(w) + s(i)*g(w)*d_0_temp);
              hastings = hastings - (Y(i, w)*log(s(i)*g(w)*d_0(w)) - (phi(w) + Y(i, w))*log(phi(w) + s(i)*g(w)*d_0(w)));
            }
            else
            {
              hastings = hastings + Y(i, w)*log(s(i)*g(w)*(d_0_temp*delta(w))) - (phi(w) + Y(i, w))*log(phi(w) + s(i)*g(w)*(d_0_temp*delta(w)));
              hastings = hastings - (Y(i, w)*log(s(i)*g(w)*(d_0(w)*delta(w))) - (phi(w) + Y(i, w))*log(phi(w) + s(i)*g(w)*(d_0(w)*delta(w))));
            }
          }
        }
        hastings = hastings + (a_d - 1)*log(d_0_temp) - b_d*d_0_temp;
        hastings = hastings - ((a_d - 1)*log(d_0(w)) - b_d*d_0(w));
        if(hastings >= log(double(rand()%10001)/10000))
        {
          d_0(w) = d_0_temp;
          if(ii >= burn) {
            accept_d_0++;
          }
        }
      }
    }
    
    // Update delta
    for(w = 0; w < W; w++)
    {
      if(flag(w) == 0)
      {
        do {
          delta_temp = rnorm_trunc(delta(w), tau_delta, 0, fold_max);
        } while (delta_temp == 0);
        hastings = 0;
        try_delta++;
        for(i = 0; i < n; i++)
        {
          if(H(i, w) == 0 && ip(i) == 1) {
            hastings = hastings + Y(i, w)*log(s(i)*g(w)*(d_0(w)*delta_temp)) - (phi(w) + Y(i, w))*log(phi(w) + s(i)*g(w)*(d_0(w)*delta_temp));
            hastings = hastings - (Y(i, w)*log(s(i)*g(w)*(d_0(w)*delta(w))) - (phi(w) + Y(i, w))*log(phi(w) + s(i)*g(w)*(d_0(w)*delta(w))));
          }
        }
        hastings = hastings + (-(delta_temp - mu(z(w)))*(delta_temp - mu(z(w)))/2/sigma2(z(w)));
        hastings = hastings - (-(delta(w) - mu(z(w)))*(delta(w) - mu(z(w)))/2/sigma2(z(w)));
        if(hastings >= log(double(rand()%10001)/10000))
        {
          delta(w) = delta_temp;
          if(ii >= burn) {
            accept_delta++;
          }
        }
      }
    }
    
    // Update A
    if (HMM)
    {
      for(q = 0; q < Q; q++)
      {
        sum_temp = 0;
        for(qq = 0; qq < Q; qq++)
        {
          prob_temp(qq) = rgamma(1, C(qq, q) + 1, 1)(0);
          sum_temp = sum_temp + prob_temp(qq);
        }
        for(qq = 0; qq < Q; qq++)
        {
          A(qq, q) = prob_temp(qq)/sum_temp;
        }
      }
    }
    else
    {
      temp = rbeta(1, 1 + z_2, 1 + W - z_2)(0);
      A(0, 0) = 1 - temp;
      A(1, 0) = 1 - temp;
      A(0, 1) = temp;
      A(1, 1) = temp;
    }
    
    // Update mu and sigma2
    for(q = 0; q < Q; q++)
    {
      count_temp = 0;
      sum_temp = 0;
      for(w = 0; w < W; w++)
      {
        if(z(w) == q && flag(w) == 0)
        {
          count_temp = count_temp + 1;
          sum_temp = sum_temp + delta(w);
        }
      }
      mu(q) = rnorm_trunc(1/(1/tau_0/tau_0 + count_temp/sigma2(q))*(mu_0/tau_0/tau_0 + sum_temp/sigma2(q)), 1/(1/tau_0/tau_0 + count_temp/sigma2(q)), mu_limit(0, q), mu_limit(1, q)); 
      sum_temp = 0;
      for(w = 0; w < W; w++)
      {
        if(z(w) == q && flag(w) == 0)
        {
          sum_temp = sum_temp + (delta(w) - mu(q))*(delta(w) - mu(q));
        }
      }
      sigma2(q) = 1/rgamma(1, a_0 + count_temp/2, 1/(b_0 + sum_temp/2))(0);
    }
    
    // Update z
    for(q = 0; q < Q; q++)
    {
      for(qq = 0; qq < Q; qq++)
      {
        C(q, qq) = 0;
      }
    }
    w = W - 1;
    while (flag(w) == 1) {
      w = w - 1;
    }
    w_temp = w;
    for (q = 0; q < Q; q++)
    {
      prob_temp(q) = -log(sigma2(q)) - (delta(w) - mu(q))*(delta(w) - mu(q))/2/sigma2(q);
    }
    sum_temp = 0;
    max_temp = max(prob_temp);
    for(q = 0; q < Q; q++)
    {
      prob_temp(q) = prob_temp(q) - max_temp;
      prob_temp(q) = exp(prob_temp(q));
      sum_temp = sum_temp + prob_temp(q);
    }
    for (q = 0; q < Q; q++)
    {
      prob_temp(q) = prob_temp(q)/sum_temp;
    }
    z(w) = RcppArmadillo::sample(state, 1, true, prob_temp)(0);
    for(w = w_temp - 1; w >= 0; w--)
    {
      if(flag(w) == 0)
      {
        for (q = 0; q < Q; q++)
        {
          prob_temp(q) = -log(sigma2(q)) - (delta(w) - mu(q))*(delta(w) - mu(q))/2/sigma2(q) + log(A(z(w + 1), q));
        }
        sum_temp = 0;
        max_temp = max(prob_temp);
        for(q = 0; q < Q; q++)
        {
          prob_temp(q) = prob_temp(q) - max_temp;
          prob_temp(q) = exp(prob_temp(q));
          sum_temp = sum_temp + prob_temp(q);
        }
        for (q = 0; q < Q; q++)
        {
          prob_temp(q) = prob_temp(q)/sum_temp;
        }
        z(w) = RcppArmadillo::sample(state, 1, true, prob_temp)(0);
      }
    }
    for(w = 1; w < W; w++)
    {
      if(flag(w) == 0)
      {
        C(z(w), z(w - 1))++;
      }
    }
    z_2 = 0;
    for(w = 0; w < W; w++)
    {
      z_2 = z_2 + z(w);
    }
    
    // Monitor the process
    if(ii*100/iter == count)
    {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
    if(store)
    {
      pi_store(ii) = pi;
      for(w = 0; w < W; w++)
      {
        d_0_store(ii, w) = d_0(w);
        delta_store(ii, w) = delta(w);
        phi_store(ii, w) = phi(w);
      }
      for(q = 0; q < Q; q++)
      {
        mu_store(ii, q) = mu(q);
        sigma2_store(ii, q) = sigma2(q);
        for(qq = 0; qq < Q; qq++)
        {
          A_store(ii, q*Q + qq) = A(q, qq);
        }
      }
    }
    z_sum(ii) = z_2;
    if(ii >= burn)
    {
      pi_mean = pi_mean + pi;
      for(w = 0; w < W; w++)
      {
        d_0_mean(w) = d_0_mean(w) + d_0(w);
        delta_mean(w) = delta_mean(w) + delta(w);
        phi_mean(w) = phi_mean(w) + phi(w);
        for(i = 0; i < n; i++)
        {
          H_ppi(i, w) = H_ppi(i, w) + H(i, w);
        }
        z_ppi(z(w), w)++;
      }
      for(q = 0; q < Q; q++)
      {
        mu_mean(q) = mu_mean(q) + mu(q);
        sigma2_mean(q) = sigma2_mean(q) + sigma2(q);
        for(qq = 0; qq < Q; qq++)
        {
          A_mean(q, qq) = A_mean(q, qq) + A(q, qq);
        }
      }
    }
  }
  accept_d_0 = accept_d_0/try_d_0;
  accept_delta = accept_delta/try_delta;
  accept_phi = accept_phi/try_phi;
  pi_mean = pi_mean/(iter - burn);
  for(w = 0; w < W; w++)
  {
    d_0_mean(w) = d_0_mean(w)/(iter - burn);
    delta_mean(w) = delta_mean(w)/(iter - burn);
    phi_mean(w) = phi_mean(w) /(iter - burn);
    for(i = 0; i < n; i++)
    {
      H_ppi(i, w) = H_ppi(i, w)/(iter - burn);
    }
  }
  for(q = 0; q < Q; q++)
  {
    mu_mean(q) = mu_mean(q)/(iter - burn);
    sigma2_mean(q) = sigma2_mean(q)/(iter - burn);
    for(qq = 0; qq < Q; qq++)
    {
      A_mean(q, qq) = A_mean(q, qq)/(iter - burn);
    }
    for(w = 0; w < W; w++)
    {
      z_ppi(q, w) = z_ppi(q, w)/(iter - burn);
    }
  }
  if(store)
  {
    return Rcpp::List::create(Rcpp::Named("d_0") = d_0_mean, Rcpp::Named("delta") = delta_mean, Rcpp::Named("phi") = phi_mean, Rcpp::Named("pi") = pi_mean, Rcpp::Named("H") = H_ppi, Rcpp::Named("accept_d_0") = accept_d_0, Rcpp::Named("accept_delta") = accept_delta, Rcpp::Named("accept_phi") = accept_phi, Rcpp::Named("A") = A_mean, Rcpp::Named("mu") = mu_mean, Rcpp::Named("sigma2") = sigma2_mean, Rcpp::Named("z") = z_ppi, Rcpp::Named("d_0_store") = d_0_store, Rcpp::Named("delta_store") = delta_store, Rcpp::Named("phi_store") = phi_store, Rcpp::Named("pi_store") = pi_store, Rcpp::Named("A_store") = A_store, Rcpp::Named("mu_store") = mu_store, Rcpp::Named("sigma2_store") = sigma2_store, Rcpp::Named("z_sum") = z_sum);
  }
  else
  {
    return Rcpp::List::create(Rcpp::Named("d_0") = d_0_mean, Rcpp::Named("delta") = delta_mean, Rcpp::Named("phi") = phi_mean, Rcpp::Named("pi") = pi_mean, Rcpp::Named("H") = H_ppi, Rcpp::Named("accept_d_0") = accept_d_0, Rcpp::Named("accept_delta") = accept_delta, Rcpp::Named("accept_phi") = accept_phi, Rcpp::Named("A") = A_mean, Rcpp::Named("mu") = mu_mean, Rcpp::Named("sigma2") = sigma2_mean, Rcpp::Named("z") = z_ppi, Rcpp::Named("z_sum") = z_sum);
  }
}

// [[Rcpp::export]]
double rnorm_trunc(double mu, double sigma, double lower, double upper)
{
  int change;
  double a, b;
  double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
  double z, tmp, lograt;
  
  change = 0;
  a = (lower - mu)/sigma;
  b = (upper - mu)/sigma;
  
  // First scenario
  if( (a == R_NegInf)||(b == R_PosInf))
  {
    if(a == R_NegInf)
    {
      change = 1;
      a = -b;
      b = R_PosInf;
    }
    // The two possibilities for this scenario
    if(a <= 0.45) z = norm_rs(a, b);
    else z = exp_rs(a, b);
    if(change) z = -z;
  }
  
  // Second scenario
  else if((a*b) <= 0.0)
  {
    // The two possibilities for this scenario
    if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
    {
      z = norm_rs(a, b);
    }
    else z = unif_rs(a,b);
  }
  
  // Third scenario
  else
  {
    if(b < 0)
    {
      tmp = b; b = -a; a = -tmp; change = 1;
    }
    
    lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
    if(lograt <= logt2)
    {
      z = unif_rs(a,b);
    }
    else if((lograt > logt1)&&(a < t3))
    {
      z = half_norm_rs(a,b);
    }
    else
    {
      z = exp_rs(a,b);
    }
    if(change)
    {
      z = -z;
    }
  }
  double output;
  output = sigma*z + mu;
  return (output);
}

// [[Rcpp::export]]
double exp_rs(double a, double b)
{
  double  z, u, rate;
  rate = 1/a;
  
  // Generate a proposal on (0, b-a)
  z = R::rexp(rate);
  while(z > (b-a))
  {
    z = R::rexp(rate);
  }
  u = R::runif(0.0, 1.0);
  
  while( log(u) > (-0.5*z*z))
  {
    z = R::rexp(rate);
    while(z > (b-a))
    {
      z = R::rexp(rate);
    }
    u = R::runif(0.0,1.0);
  }
  return(z+a);
}

// [[Rcpp::export]]
double unif_rs(double a, double b)
{
  double xstar, logphixstar, x, logu;
  
  // Find the argmax (b is always >= 0)
  // This works because we want to sample from N(0,1)
  if(a <= 0.0) 
  {
    xstar = 0.0;
  }
  else 
  {
    xstar = a;
  }
  logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);
  
  x = R::runif(a, b);
  logu = log(R::runif(0.0, 1.0));
  while(logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
  {
    x = R::runif(a, b);
    logu = log(R::runif(0.0, 1.0));
  }
  return x;
}

// [[Rcpp::export]]
double half_norm_rs(double a, double b)
{
  double x;
  x = fabs(norm_rand());
  while((x<a)||(x>b)) 
  {
    x = fabs(norm_rand());
  }
  return x;
}

// [[Rcpp::export]]
double norm_rs(double a, double b)
{
  double x;
  x = Rf_rnorm(0.0, 1.0);
  while((x < a)||(x > b)) 
  {
    x = norm_rand();
  }
  return x;
}