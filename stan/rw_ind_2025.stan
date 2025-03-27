//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> N; // number of annual samples
  int L; // years covered by time-series
  int J; // number of stocks
  array[N] int J_i; //index of stocks
  array[N] int J_ii; // index of stock-year combination
  array[N] real R_S; // log(recruits per spawner)
  matrix[N,J] S; // design matrix of spawners in time T
  vector[N] X1; // vector of covariate 1 in time T
  vector[N] X2; // vector of covariate 2 in time T
  array[J] real pSmax_mean; // priors on smax - based on observed spawner abundance
  array[J] real pSmax_sig; // sigma prior on smax
}
transformed data{
  vector[J] smax_pr;
vector[J] smax_pr_sig;

for(j in 1:J){
smax_pr_sig[j]=sqrt(log(1+((pSmax_sig[j])^2)/((pSmax_mean[j])^2))); //this converts sigma on the untransformed scale to a log scale
smax_pr[j]=log(pSmax_mean[j])-0.5*(smax_pr_sig[j])^2; //convert smax prior to per capita slope - transform to log scale with bias correction
}
}
parameters {
  vector[J] g0; // initial X1 effect (stock-specific)
  vector[J] k0; // initial X2 effect (stock-specific)
  vector[J] log_a; // intrinsic productivity - make hierarchical?
  vector<lower=0>[J] Smax; // spawners where recruitment is maximized
  
  // variance components
  vector<lower = 0>[J] sigma_t; // total variance
  vector<lower = 0, upper = 1>[J] p_rw_g; // proportion gamma process variance
  vector<lower = 0, upper = 1>[J] p_rw_k; // proportion kappa process variance
  
  matrix[L-1,J] z_dev_g; // year-to-year deviations in gamma
  matrix[L-1, J] z_dev_k; // year-to-year deviations in kappa
  
}
transformed parameters{
  
  simplex[3] F; // simplex to bound sigma fractions to sum to 1
  vector<lower = 0>[J] sigma; //annual observation error
  vector<lower = 0>[J] sigma_g; // variance in gamma among years
  vector<lower = 0>[J] sigma_k; // variance in kappa among years

  for(j in 1:J){
    F[1] = 1-(p_rw_g[j] + p_rw_k[j]);
    F[2] = p_rw_g[j];
    F[3] = p_rw_k[j];
    sigma[j] = sigma_t[j].*F[1]; // obs error
    sigma_g[j] = sigma_t[j].*F[2]; // process error = prop total sig attributed to gamma RW
    sigma_k[j] = sigma_t[j].*F[3]; // process error = prop total sig attributed to kappa RW
  } 
  
  
  vector[J] b; // capacity rate
  b = 1.0./Smax; 

  //stock state process
  matrix[L,J] g_t; // gamma at time T 
  matrix[L,J] k_t; // kappa at time T
  
  g_t[1,] = to_row_vector(g0); 
  k_t[1,] = to_row_vector(k0);
  
  for(t in 2:L){
    g_t[t,] = g_t[t-1,] + z_dev_g[t-1,].*to_row_vector(sigma_g);
    k_t[t,] = k_t[t-1,] + z_dev_k[t-1,].*to_row_vector(sigma_k);
  }
  
  
}
model {
  g0 ~ normal(0,5); 
  k0 ~ normal(0,5); 
  log_a ~ normal(1.5,2); 

  
  //capacity for each stock
  for(j in 1:J) Smax[j] ~ lognormal(smax_pr[j], smax_pr_sig[j]);
  
  for(j in 1:J) z_dev_g[,j] ~ std_normal();
  for(j in 1:J) z_dev_k[,j] ~ std_normal();
  
  // variance terms
  sigma_t ~ gamma(2,1);
  p_rw_g ~ beta(2,5); // half normal on variance (lower limit of zero) // try (2,10)
  p_rw_k ~ beta(2,5); // half normal on variance (lower limit of zero)
  
  // Ricker
  R_S ~ normal(log_a - S*b + to_vector(g_t[J_ii]).*X1 + to_vector(k_t[J_ii]).*X2, sigma[J_i]);
  
  
}
generated quantities {
  vector[N] log_lik;
  vector[L] mu_g_t;
  vector[L] mu_k_t;
  vector[L] mu_umsy;
  
  for(l in 1:L){
      mu_g_t[l] = sum(g_t[l,])/J;
      mu_k_t[l] = sum(k_t[l,])/J;
  }
  
  for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|log_a - S[i]*b + to_vector(g_t)[J_ii[i]]*X1 + to_vector(k_t)[J_ii[i]]*X2, sigma[J_i[i]]);

  
}

