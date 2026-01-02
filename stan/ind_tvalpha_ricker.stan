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
  vector[J] k0; // initial kappa (stock-specific)
  vector<lower = 0>[J] Smax; // spawners where recruitment is maximized
 //variance components
  vector<lower = 1e-6, upper=10>[J] sigma_t; //total variance
  //vector<lower=0,upper=1>[J] p_rw_g; //proportion of variance that is gamma process variance
  //vector<lower=0,upper=1>[J] p_rw_k; //proportion of variance that is kappa process variance
  array[J] simplex[3] p_sig; //proportions of sigma_t that are obs error, gamma rw, kappa rw
//MVN parameters
  matrix[L-1,J] z_dev_g;  // deviations in gamma from year-to-year
  matrix[L-1,J] z_dev_k;  // deviations in kappa from year-to-year
}
transformed parameters{
 //capacity rate (ricker b)
  vector[J] b;
  b=1.0./Smax;
 //observation and process error
 vector<lower = 0>[J] sigma; //annual observation error
 vector<lower = 0>[J] sigma_g; // variance in gamma among years
 vector<lower = 0>[J] sigma_k; // variance in kappa among years

  for (j in 1:J) {
  sigma[j]   = sigma_t[j] * p_sig[j][1];// obs error
  sigma_g[j] = sigma_t[j] * p_sig[j][2];// process error = prop total sig attributed to gamma RW
  sigma_k[j] = sigma_t[j] * p_sig[j][3];// process error = prop total sig attributed to kappa RW
  }

  //stock-specific random walk process
  matrix[L,J] g_t; // gamma coefficient at time t (local + global)
  matrix[L,J] k_t; // kappa coefficient at time t (local + global)

  g_t[1,] = to_row_vector(g0); //initial gamma
  k_t[1,] = to_row_vector(k0); //initial kappa

  for(t in 2:L){
    g_t[t,] = g_t[t-1,] + z_dev_g[t-1,].*to_row_vector(sigma_g); // stock-specific gamma random walk
    k_t[t,] = k_t[t-1,] + z_dev_k[t-1,].*to_row_vector(sigma_k); // stock-specific kappa random walk
  }


}
model{
  //priors
  log_a ~ normal(1.5,5); //stock-specific intrinsic productivity at the initial time-step of all series
  g0 ~ normal(0,5);
  k0 ~ normal(0,5);

 //capacity for each stock
  for(j in 1:J) Smax[j] ~ normal(pSmax_mean[j],pSmax_sig[j]);

  //annual deviations in random walk
  for(j in 1:J){
    z_dev_g[,j] ~ std_normal();
    z_dev_k[,j] ~ std_normal();
  }

  //variance terms
  sigma_t ~ gamma(2,1);
  for(j in 1:J){
    p_sig[j] ~ dirichlet(rep_vector(1.0,3)); //uniform prior on the simplex
  }

  //Modified Ricker
 R_S ~ normal(log_a[J_i] - S*b + to_vector(g_t)[J_ii] .*X1 + to_vector(k_t)[J_ii] .*X2, sigma[J_i]);

}
generated quantities{
vector[N] log_lik;
matrix[L,R] mu_g_r;
matrix[L,R] mu_k_r;
row_vector[R] R_n; //n stocks per region

for(l in 1:L){

  for(r in 1:R){
  R_n[r] = 0; //start with 0 stocks
  mu_g_r[l,]=rep_row_vector(0.0,R);
  mu_k_r[l,]=rep_row_vector(0.0,R);
  }
    for(j in 1:J){
    int r = J_or[j];
    R_n[r] += 1; //count stocks
    mu_g_r[l,r]+=g_t[l,j]; //summing coefficients
    mu_k_r[l,r]+=k_t[l,j];

	}
   mu_g_r[l,]=mu_g_r[l,]./R_n; //dividing to get the mean
   mu_k_r[l,]=mu_k_r[l,]./R_n;
}

for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|log_a - S[i]*b + to_vector(g_t)[J_ii[i]]*X1[i] + to_vector(k_t)[J_ii[i]]*X2[i], sigma[J_i[i]]);

}
