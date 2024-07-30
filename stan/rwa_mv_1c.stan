// HH/DG version of MV RW
data{
  int<lower=1> N; //number of annual samples 
  int L; //years covered by time-series
  int J; //number of stocks
  int J_i[N]; //index of stocks
  int J_ii[N]; //index of stock-year combination (eg. stock 1-year 1, stock 1-year 2, etc)
  vector[N] R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  vector[N] X; // vector of covariate X in time T
  vector[J] pSmax_mean; //priors on smax - based on observed spawner abundance
  vector[J] pSmax_sig;
}
transformed data{
vector[J] logbeta_pr;
vector[J] logbeta_pr_sig;

for(t in 1:J){
logbeta_pr_sig[t]=sqrt(log(1+((1/pSmax_sig[t])*(1/pSmax_sig[t]))/((1/pSmax_mean[t])*(1/pSmax_mean[t])))); //this converts sigma on the untransformed scale to a log scale
logbeta_pr[t]=log(1/pSmax_mean[t])-0.5*logbeta_pr_sig[t]*logbeta_pr_sig[t]; //convert smax prior to per capita slope - transform to log scale with bias correction
}
}
parameters{
 //row_vector[J] log_a_s0; // initial productivity (stock-specific)
 row_vector[J] g_s0; // gamma initial (stock-specific) -HH
 vector<upper = 0>[J] log_b; // rate capacity - fixed in this
 vector[J] log_a; // instrinsic growth - also fixed in this version -HH
 
 //variance components
 vector<lower = 0>[J] sigma;

 
//MVN parameters  
  cholesky_factor_corr[J] Lcorr;
  //vector[J] z_dev_s0;
  vector[J] z_dev_g0; // does this parameter do anything? was unlabeled -HH 
  //matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  matrix[L-1,J] z_dev_g; // deviations in gamma from year-to-year
  //vector<lower = 0>[J] sigma_a_s; //variance in stock productivity among years //try sharing?
  real<lower=0> sigma_g_s; // -- SIMPLIFYING CHANGE to shared sigma
  //vector<lower = 0>[J] sigma_g_s; //variance in gammas among years

}
transformed parameters{
  vector[J] b; //capacity rate
  //stock state process
  //matrix[L,J] log_a_t; //stock productivity over time
  matrix[L,J] g_t; //covar effect over time -HH
  //matrix[L-1,J] a_dev; //stock-level deviations in year-to-year productivity
  matrix[L-1,J] g_dev; //stock-level deviations in year-to-year covar effect -HH

  b=exp(log_b);
  
  //initial productivities
  g_t[1,] = g_s0;
  
  for(t in 1:L-1){
   g_dev[t,] = (diag_pre_multiply(rep_vector(sigma_g_s, J), Lcorr) * to_vector(z_dev_g[t,]))';
  }
  
  for(t in 2:L){
    g_t[t,] = g_t[t-1] + g_dev[t-1,]; //final estimate of stock prod. (global + stock)
  }
}  
model{
  
  for(j in 1:J)log_b[j] ~ normal(logbeta_pr[j],logbeta_pr_sig[j]); //capacity with informative prior for each stock
  //log_a_s0 ~ normal(1.5,5);//initial prod. for each stock
  log_a ~ normal(1.5,2); //intercept
  g_s0 ~ normal(0,5); // -HH

  Lcorr ~ lkj_corr_cholesky(1.0); // prior for correlation of process deviances
 
  to_vector(z_dev_g) ~ std_normal(); //global deviation in prod.
  z_dev_g0 ~ std_normal(); // - HH does this parameter do anything?
  
  //variance terms
  sigma ~ gamma(2,2);
 // sigma_a_s ~ gamma(2,2);
  sigma_g_s ~ gamma(1,2); //-HH
  
  // R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  R_S ~ normal(log_a[J_i] - S*b + to_vector(g_t)[J_ii] .*X, sigma[J_i]); // HH
  
}
generated quantities{
corr_matrix[J] Cor_J = multiply_lower_tri_self_transpose(Lcorr);
vector[N] log_lik;
vector[L] gamma_t_m; //average gamma through time

for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(g_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);
for(t in 1:L) gamma_t_m[t]=sum(g_t[t,])/J; //arithmetic mean of log_alpha at each time-step
}
