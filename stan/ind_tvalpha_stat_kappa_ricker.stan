data{
  int<lower=1> N;//number of annual samples
  int L; //years covered by time-series
  int J;//number of stocks
  int R;//number of ocean regions
  array[N] int J_i;//index of stocks
  array[N] int J_ii;//index of stock-year combination
  array[N] real R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners
  vector[N] X1; // vector of covariate 1 in time T
  vector[N] X2; // vector of covariate 2 in time T
  array[J] int J_or;//index of ocean region by stock identity
  array[J] real pSmax_mean; //priors on smax - based on observed spawner abundance
  array[J] real pSmax_sig;
}
parameters{
  vector[J] log_a;// productivity (stock-specific)
  vector[J] g0; // initial gamma (stock-specific)
  vector[J] kappa; // stock-specific kappa (stationary -- no longer a random-walk initial value) // CHANGED: kappa is now a single stationary coefficient, not a random-walk starting point
  vector<lower = 0>[J] Smax; // spawners where recruitment is maximized
 //variance components
  vector<lower = 1e-6, upper=10>[J] sigma_t; //total variance
  array[J] simplex[2] p_sig; //proportions of sigma_t that are obs error, gamma rw // CHANGED: simplex now has 2 parts (obs error, gamma rw) instead of 3, since kappa no longer has a process-error component
//MVN parameters
  matrix[L-1,J] z_dev_g;  // deviations in gamma from year-to-year
  // CHANGED: z_dev_k (year-to-year kappa deviations) removed -- kappa no longer varies over time
}
transformed parameters{
 //capacity rate (ricker b)
  vector[J] b;
  b=1.0./Smax;
 //observation and process error
 vector<lower = 0>[J] sigma; //annual observation error
 vector<lower = 0>[J] sigma_g; // variance in gamma among years
 // CHANGED: sigma_k removed -- kappa has no process-error variance since it no longer varies among years

  for (j in 1:J) {
  sigma[j]   = sigma_t[j] * p_sig[j][1];// obs error
  sigma_g[j] = sigma_t[j] * p_sig[j][2];// process error = prop total sig attributed to gamma RW
  }

  //stock-specific random walk process
  matrix[L,J] g_t; // gamma coefficient at time t (local + global)
  // CHANGED: k_t (time-varying kappa matrix) removed -- kappa is now a single stationary vector[J], no time dimension

  g_t[1,] = to_row_vector(g0); //initial gamma

  for(t in 2:L){
    g_t[t,] = g_t[t-1,] + z_dev_g[t-1,].*to_row_vector(sigma_g); // stock-specific gamma random walk
  }
  // CHANGED: kappa random-walk loop removed -- kappa is constant over time, so no t-loop is needed for it


}
model{
  //priors
  log_a ~ normal(1.5,5); //stock-specific intrinsic productivity at the initial time-step of all series
  g0 ~ normal(0,5);
  kappa ~ normal(0,5); // CHANGED: prior on stationary kappa (replaces the k0 prior)

 //capacity for each stock
  for(j in 1:J) Smax[j] ~ normal(pSmax_mean[j],pSmax_sig[j]);

  //annual deviations in random walk
  for(j in 1:J){
    z_dev_g[,j] ~ std_normal();
    // CHANGED: z_dev_k prior removed -- kappa no longer has year-to-year deviations
  }

  //variance terms
  sigma_t ~ gamma(2,1);
  for(j in 1:J){
    p_sig[j] ~ dirichlet(rep_vector(1.0,2)); //uniform prior on the simplex // CHANGED: simplex now length 2 (obs error, gamma rw)
  }

  //Modified Ricker
 R_S ~ normal(log_a[J_i] - S*b + to_vector(g_t)[J_ii] .*X1 + kappa[J_i] .*X2, sigma[J_i]); // CHANGED: kappa applied directly by stock (kappa[J_i]), not looked up through a time-varying k_t matrix

}
generated quantities{
vector[N] log_lik;
matrix[L,R] mu_g_r;
vector[R] mu_k_r; // CHANGED: regional kappa mean is now a single vector (no time dimension), not a matrix[L,R]
row_vector[R] R_n; //n stocks per region

for(l in 1:L){

  for(r in 1:R){
  R_n[r] = 0; //start with 0 stocks
  mu_g_r[l,]=rep_row_vector(0.0,R);
  }
    for(j in 1:J){
    int r = J_or[j];
    R_n[r] += 1; //count stocks
    mu_g_r[l,r]+=g_t[l,j]; //summing coefficients

	}
   mu_g_r[l,]=mu_g_r[l,]./R_n; //dividing to get the mean
}

// CHANGED: regional mean kappa computed once (outside the year loop), since kappa is stationary
{
  vector[R] k_n;
  for(r in 1:R){
    k_n[r] = 0;
    mu_k_r[r] = 0.0;
  }
  for(j in 1:J){
    int r = J_or[j];
    k_n[r] += 1;
    mu_k_r[r] += kappa[j];
  }
  mu_k_r = mu_k_r ./ k_n;
}

for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|log_a - S[i]*b + to_vector(g_t)[J_ii[i]]*X1[i] + kappa[J_i[i]]*X2[i], sigma[J_i[i]]); // CHANGED: kappa indexed by stock, not by the (now-removed) k_t matrix

}
