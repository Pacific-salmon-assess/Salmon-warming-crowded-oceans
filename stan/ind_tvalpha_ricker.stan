data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  array[N] int J_i;//index of stocks
  array[N] int J_ii;//index of stock-year combination
  array[N] real R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  array[J] real pSmax_mean; //priors on smax - based on observed spawner abundance
  array[J] real pSmax_sig;
}
transformed data{
vector[J] smax_pr;
vector[J] smax_pr_sig;

for(j in 1:J){
smax_pr_sig[j]=sqrt(log(1+((pSmax_sig[j])^2)/((pSmax_mean[j])^2))); //this converts sigma on the untransformed scale to a log scale
smax_pr[j]=log(pSmax_mean[j])-0.5*(smax_pr_sig[j])^2; //convert smax prior to per capita slope - transform to log scale with bias correction
}
}
parameters{
  vector[J] log_a0_s;// initial productivity (stock-specific)
  vector<lower = 0>[J] Smax; // spawners where recruitment is maximized

 //variance components  
 vector<lower = 0>[J] sigma_t; //total variance
  vector<lower=0,upper=1>[J] p_rw; //proportion of variance that is productivity process variance
//MVN parameters  
  matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  
}
transformed parameters{
 vector<lower = 0>[J] sigma; //annual observation error
 vector<lower = 0>[J] sigma_a; //variance in stock productivity among years
 vector[J] b; //capacity rate
  //global state process
  //stock state process
  matrix[L,J] log_a_t; //stock productivity at time t (local + global)
  
 sigma=sigma_t.*(1-p_rw); //obs error = 1 - prop. attributed to random walk in productivity
 sigma_a=sigma_t.*p_rw; //proc. error = prop. of sigma total attributed to rw in productivity
 b=1.0./Smax;
  
  //initial productivities
  
  log_a_t[1,] = to_row_vector(log_a0_s); //global prod + stock deviation from global prod.
  for(t in 2:L){
    log_a_t[t,] = log_a_t[t-1,] + z_dev_s[t-1,].*to_row_vector(sigma_a); //stock-specific deviation in random walk
  }
    
    
}  
model{
  //priors
  log_a0_s ~ normal(1.5,5); //stock-specific intrinsic productivity at the initial time-step of all series
  
//capacity for each stock
  for(j in 1:J) Smax[j] ~ lognormal(smax_pr[j],smax_pr_sig[j]);
  for(j in 1:J) z_dev_s[,j] ~ std_normal(); //annual deviation in prod.
 
  //variance terms
  sigma_t ~ gamma(2,1);
  p_rw ~ beta(2,5); //half normal on variance (lower limit of zero)
 
 R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  
}
generated quantities{
vector[N] log_lik;
matrix[L,J] Smsy;
matrix[L,J] Umsy;
vector[L] mu_log_a;
vector[L] mu_umsy;

for(l in 1:L){
    mu_log_a[l] = sum(log_a_t[l,])/J;
	for(j in 1:J){
	Umsy[l,j] = 1-lambert_w0(exp(1-log_a_t[l,j]));
    Smsy[l,j] = (1-lambert_w0(exp(1-log_a_t[l,j])))/b[j];
	}
}

for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(log_a_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);

}