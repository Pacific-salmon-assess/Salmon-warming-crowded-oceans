// HH/DG version of MV RW - with a second covar added
data {
  int N; // number of observations in data
  int Ng; // gammas/kappas to be estimated: sum( L[r] x J[r] )
  int R; // number of regions, = 4
  int J_tot; // total number of stocks
  int L_tot; // total years estimated, sum of L[r]
  array[R] int J; // number of stocks in each region
  array[R] int L; // years to be estimated in each region
  array[N] int J_i; // index of stocks in data
  array[R] int J_c; // cumulative sum of J[r]
  array[N] int J_ii; //index of stock-year combination in data (eg. stock 1-year 1, stock 1-year 2, etc)
  array[J_tot] int j_start; //position of start of each stock in Ng
  array[J_tot] int j_end; // position of end of each stock in Ng
  array[J_tot] int j_in; // index of stocks within regions, i.e. 1...J[r=1], 1...J[r=2]
  array[Ng] int L_in; // index of year 1...L[r] of Ng length
  vector[N] R_S; // log(recruits per spawner)
  matrix[N, J_tot] S; // design matrix of spawners in time T
  vector[N] X1; // vector of covariate X1 in time T
  vector[N] X2; // vector of covariate X2 in time T
  vector[J_tot] pSmax_mean; // priors on smax - based on observed spawner abundance
  vector[J_tot] pSmax_sig; // priors on smax - based on observed spawner abundance
}
transformed data {
  vector[J_tot] logbeta_pr;
  vector[J_tot] logbeta_pr_sig;
  for (t in 1 : J_tot) {
    logbeta_pr_sig[t] = sqrt(log(1
                                 + ((1 / pSmax_sig[t]) * (1 / pSmax_sig[t]))
                                   / ((1 / pSmax_mean[t])
                                      * (1 / pSmax_mean[t])))); //this converts sigma on the untransformed scale to a log scale
    logbeta_pr[t] = log(1 / pSmax_mean[t])
                    - 0.5 * logbeta_pr_sig[t] * logbeta_pr_sig[t]; //convert smax prior to per capita slope - transform to log scale with bias correction
  }
}
parameters {
  vector[J_tot] g_s0; // gamma initial (stock-specific) 
  vector[J_tot] k_s0; // kappa initial (stock-specific) 
  vector<upper=0>[J_tot] log_b; // rate capacity - fixed in this
  vector[J_tot] log_a; // instrinsic growth - want this to be hierarchical
  //variance components
  vector<lower=0, upper=100>[J_tot] sigma;
  //MVN parameters  
  cholesky_factor_corr[J[1]] Lcorr_g_1;
  cholesky_factor_corr[J[2]] Lcorr_g_2;
  cholesky_factor_corr[J[3]] Lcorr_g_3;
  cholesky_factor_corr[J[4]] Lcorr_g_4;
  cholesky_factor_corr[J[1]] Lcorr_k_1;
  cholesky_factor_corr[J[2]] Lcorr_k_2;
  cholesky_factor_corr[J[3]] Lcorr_k_3;
  cholesky_factor_corr[J[4]] Lcorr_k_4;
  //vector[Ng] z_dev_g;
  //vector[Ng] z_dev_k;
  matrix[L[1]-1, J[1]] z_dev_g_1; // DG's suggestion to make it a matrix
  matrix[L[2]-1, J[2]] z_dev_g_2;
  matrix[L[3]-1, J[3]] z_dev_g_3;
  matrix[L[4]-1, J[4]] z_dev_g_4;
  matrix[L[1]-1, J[1]] z_dev_k_1;
  matrix[L[2]-1, J[2]] z_dev_k_2;
  matrix[L[3]-1, J[3]] z_dev_k_3;
  matrix[L[4]-1, J[4]] z_dev_k_4;
  
  //vector<lower=0>[R] sigma_g_s; //variance in gammas among years, shared among stocks w/in region
  //vector<lower=0>[R] sigma_k_s; //variance in kappas among years, shared among stocks w/in region
  real<lower=0> sigma_g_s_1;
  real<lower=0> sigma_g_s_2;
  real<lower=0> sigma_g_s_3;
  real<lower=0> sigma_g_s_4;
  real<lower=0> sigma_k_s_1;
  real<lower=0> sigma_k_s_2;
  real<lower=0> sigma_k_s_3;
  real<lower=0> sigma_k_s_4;
}

transformed parameters {
  vector[J_tot] b;
  //gamma and kappa process
  vector[Ng] g_t;
  vector[Ng] k_t;
  vector[Ng] g_dev;
  vector[Ng] k_dev;
  vector[N] ypred;
  //
  b = exp(log_b);
 
  for(i in 1:J_tot){ // set inital g and k
    g_t[j_start[i]] = g_s0[i];
    k_t[j_start[i]] = k_s0[i];
    g_dev[j_start[i]] = 0;
    k_dev[j_start[i]] = 0;
  }
  
  for (j in 1:J_c[1]){ // for each stk in region 1
    for(t in (j_start[j]+1):j_end[j]){ // for each year in stk j
      g_dev[t] =  sum((sigma_g_s_1 * Lcorr_g_1[j_in[j]]) * z_dev_g_1[L_in[t-1],j_in[j]]); // no good?
      k_dev[t] = sum((sigma_k_s_1 * Lcorr_k_1[j_in[j]]) * z_dev_k_1[L_in[t-1],j_in[j]]);
    }
  }
  for (j in (J_c[1]+1):J_c[2]){ // for each stk in region 2
    for(t in (j_start[j]+1):j_end[j]){ // for each year in stk j
      g_dev[t] = sum((sigma_g_s_2 * Lcorr_g_2[j_in[j]]) * z_dev_g_2[L_in[t-1],j_in[j]]);
      k_dev[t] = sum((sigma_k_s_2 * Lcorr_g_2[j_in[j]]) * z_dev_k_2[L_in[t-1],j_in[j]]);
    }
  }
  for (j in (J_c[2]+1):J_c[3]){ // for each stk in region 3
    for(t in (j_start[j]+1):j_end[j]){ // for each year in stk j
      g_dev[t] = sum((sigma_g_s_3 * Lcorr_g_3[j_in[j]]) * z_dev_g_3[L_in[t-1],j_in[j]]);
      k_dev[t] = sum((sigma_k_s_3 * Lcorr_g_3[j_in[j]]) * z_dev_k_3[L_in[t-1],j_in[j]]);
    }
  }
  for (j in (J_c[3]+1):J_c[4]){ // for each stk in region 4
    for(t in (j_start[j]+1):j_end[j]){ // for each year in stk j
      g_dev[t] = sum((sigma_g_s_4 * Lcorr_g_4[j_in[j]]) * z_dev_g_4[L_in[t-1],j_in[j]]);
      k_dev[t] = sum((sigma_k_s_4 * Lcorr_g_4[j_in[j]]) * z_dev_k_4[L_in[t-1],j_in[j]]);
    }
  } 
  
  for(j in 1:J_tot){
    for(t in (j_start[j]+1):j_end[j]){
      g_t[t] = g_t[t-1] + g_dev[t];
      k_t[t] = k_t[t-1] + k_dev[t];
    }
  }
  ypred = log_a[J_i] - S * b + g_t[J_ii] .* X1 + k_t[J_ii] .* X2;
  //print(ypred);
}
model {
  for (j in 1 : J_tot) {
    log_b[j] ~ normal(logbeta_pr[j], logbeta_pr_sig[j]);
  } // informative priors for b, by stock
  log_a ~ normal(1.5, 2);
  g_s0 ~ normal(0, 5);
  k_s0 ~ normal(0, 5);
  // priors for correlations of g/k deviances
  Lcorr_g_1 ~ lkj_corr_cholesky(1.0);
  Lcorr_g_2 ~ lkj_corr_cholesky(1.0);
  Lcorr_g_3 ~ lkj_corr_cholesky(1.0);
  Lcorr_g_4 ~ lkj_corr_cholesky(1.0);
  Lcorr_k_1 ~ lkj_corr_cholesky(1.0);
  Lcorr_k_2 ~ lkj_corr_cholesky(1.0);
  Lcorr_k_3 ~ lkj_corr_cholesky(1.0);
  Lcorr_k_4 ~ lkj_corr_cholesky(1.0);
  // Z deviations
  //z_dev_g ~ std_normal();
  //z_dev_k ~ std_normal();
  to_vector(z_dev_g_1) ~ std_normal(); // separate by region
  to_vector(z_dev_g_2) ~ std_normal();
  to_vector(z_dev_g_3) ~ std_normal();
  to_vector(z_dev_g_4) ~ std_normal();
  to_vector(z_dev_k_1) ~ std_normal();
  to_vector(z_dev_k_2) ~ std_normal();
  to_vector(z_dev_k_3) ~ std_normal();  
  to_vector(z_dev_k_4) ~ std_normal();
// variance terms
  sigma ~ gamma(2, 2);
  //sigma_g_s ~ gamma(1, 2); 
  //sigma_k_s ~ gamma(1, 2);
  sigma_g_s_1 ~ gamma(1,2); // separate by region
  sigma_g_s_2 ~ gamma(1,2);
  sigma_g_s_3 ~ gamma(1,2);
  sigma_g_s_4 ~ gamma(1,2);
  sigma_k_s_1 ~ gamma(1,2);
  sigma_k_s_2 ~ gamma(1,2);
  sigma_k_s_3 ~ gamma(1,2);
  sigma_k_s_4 ~ gamma(1,2);
  
  // Likelihood
  R_S ~ normal(ypred,
               sigma[J_i]);
}
generated quantities {
  // Output correlation matrices
  corr_matrix[J[1]] Cor_J_g_1 = multiply_lower_tri_self_transpose(Lcorr_g_1);
  corr_matrix[J[2]] Cor_J_g_2 = multiply_lower_tri_self_transpose(Lcorr_g_2);
  corr_matrix[J[3]] Cor_J_g_3 = multiply_lower_tri_self_transpose(Lcorr_g_3);
  corr_matrix[J[4]] Cor_J_g_4 = multiply_lower_tri_self_transpose(Lcorr_g_4);
  corr_matrix[J[1]] Cor_J_k_1 = multiply_lower_tri_self_transpose(Lcorr_k_1);
  corr_matrix[J[2]] Cor_J_k_2 = multiply_lower_tri_self_transpose(Lcorr_k_2);
  corr_matrix[J[3]] Cor_J_k_3 = multiply_lower_tri_self_transpose(Lcorr_k_3);
  corr_matrix[J[4]] Cor_J_k_4 = multiply_lower_tri_self_transpose(Lcorr_k_4);
  // Average gamma and kappa through time
  // Need new method since removed J_l, Lstart and Lend indices
  //vector[L_tot] gamma_t_m; 
  //vector[L_tot] kappa_t_m;
  //int pos; // temporary position integer
  //pos=1;
  //for (r in 1 : R) {
    //for (t in 1:) { 
    //gamma_t_m[t] = sum(g_t[]) / J[r]; //arithmetic mean of g_t at each time-step
      //kappa_t_m[t] = sum(k_t[J_l[Lstart[t] : Lend[t]]]) / J[r]; //arithmetic mean of k_t at each time-step
    //}
    //pos = L_c[r] + 1;
  //}
}
