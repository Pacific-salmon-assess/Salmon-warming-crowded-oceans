// HH/DG version of MV RW - with a second covar added
data {
  int N; // number of observations in data
  int Ng; // gammas/kappas to be estimated: sum( L[r] x J[r] )
  int R; // number of regions, = 4
  int J_tot; // total number of stocks
  int L_tot; // total years estimated, sum of L[r]
  array[R] int J; // number of stocks in each region
  array[R] int L; // years to be estimated in each region
  array[R] int L_c; // cumulative sum of L[r]
  array[N] int J_i; // index of stocks in data
  array[N] int J_ii; //index of stock-year combination in data (eg. stock 1-year 1, stock 1-year 2, etc)
  array[Ng] int J_l; // index mapping 1:L order of Ng to 1:J order
  array[L_tot] int Lstart; // start position of each year in J_l
  array[L_tot] int Lend; // end position of each year in J_l
  array[J_tot] int j_init; // position of first year in J_l for each stock
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
  vector[Ng] z_dev_g;
  vector[Ng] z_dev_k;
  vector<lower=0>[R] sigma_g_s; //variance in gammas among years, shared among stocks w/in region
  vector<lower=0>[R] sigma_k_s; //variance in kappas among years, shared among stocks w/in region
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
  // calculate first-year gammas for all stocks
  g_t[j_init] = g_s0;
  k_t[j_init] = k_s0;
  // calculate annual deviates separately for each region
  // region 1
  for (i in 2 : L_c[1]) {
    g_dev[J_l[Lstart[i] : Lend[i]]] = diag_pre_multiply(rep_vector(sigma_g_s[1],
                                                                   J[1]),
                                                        Lcorr_g_1)
                                      * z_dev_g[J_l[Lstart[i] : Lend[i]]];
    k_dev[J_l[Lstart[i] : Lend[i]]] = diag_pre_multiply(rep_vector(sigma_k_s[1],
                                                                   J[1]),
                                                        Lcorr_k_1)
                                      * z_dev_k[J_l[Lstart[i] : Lend[i]]];
    g_t[J_l[Lstart[i] : Lend[i]]] = g_t[J_l[Lstart[i - 1] : Lend[i - 1]]]
                                    + g_dev[J_l[Lstart[i] : Lend[i]]];
    k_t[J_l[Lstart[i] : Lend[i]]] = k_t[J_l[Lstart[i - 1] : Lend[i - 1]]]
                                    + k_dev[J_l[Lstart[i] : Lend[i]]];
  }
  // region 2
  for (i in (L_c[1] + 2) : L_c[2]) {
    g_dev[J_l[Lstart[i] : Lend[i]]] = diag_pre_multiply(rep_vector(sigma_g_s[2],
                                                                   J[2]),
                                                        Lcorr_g_2)
                                      * z_dev_g[J_l[Lstart[i] : Lend[i]]];
    k_dev[J_l[Lstart[i] : Lend[i]]] = diag_pre_multiply(rep_vector(sigma_k_s[2],
                                                                   J[2]),
                                                        Lcorr_k_2)
                                      * z_dev_k[J_l[Lstart[i] : Lend[i]]];
    g_t[J_l[Lstart[i] : Lend[i]]] = g_t[J_l[Lstart[i - 1] : Lend[i - 1]]]
                                    + g_dev[J_l[Lstart[i] : Lend[i]]];
    k_t[J_l[Lstart[i] : Lend[i]]] = k_t[J_l[Lstart[i - 1] : Lend[i - 1]]]
                                    + k_dev[J_l[Lstart[i] : Lend[i]]];
  }
  // region 3 
  for (i in (L_c[2] + 2) : L_c[3]) {
    g_dev[J_l[Lstart[i] : Lend[i]]] = diag_pre_multiply(rep_vector(sigma_g_s[3],
                                                                   J[3]),
                                                        Lcorr_g_3)
                                      * z_dev_g[J_l[Lstart[i] : Lend[i]]];
    k_dev[J_l[Lstart[i] : Lend[i]]] = diag_pre_multiply(rep_vector(sigma_k_s[3],
                                                                   J[3]),
                                                        Lcorr_k_3)
                                      * z_dev_k[J_l[Lstart[i] : Lend[i]]];
    g_t[J_l[Lstart[i] : Lend[i]]] = g_t[J_l[Lstart[i - 1] : Lend[i - 1]]]
                                    + g_dev[J_l[Lstart[i] : Lend[i]]];
    k_t[J_l[Lstart[i] : Lend[i]]] = k_t[J_l[Lstart[i - 1] : Lend[i - 1]]]
                                    + k_dev[J_l[Lstart[i] : Lend[i]]];
  }
  // region 4 
  for (i in (L_c[3] + 2) : L_c[4]) {
    g_dev[J_l[Lstart[i] : Lend[i]]] = diag_pre_multiply(rep_vector(sigma_g_s[4],
                                                                   J[4]),
                                                        Lcorr_g_4)
                                      * z_dev_g[J_l[Lstart[i] : Lend[i]]];
    k_dev[J_l[Lstart[i] : Lend[i]]] = diag_pre_multiply(rep_vector(sigma_k_s[4],
                                                                   J[4]),
                                                        Lcorr_k_4)
                                      * z_dev_k[J_l[Lstart[i] : Lend[i]]];
    g_t[J_l[Lstart[i] : Lend[i]]] = g_t[J_l[Lstart[i - 1] : Lend[i - 1]]]
                                    + g_dev[J_l[Lstart[i] : Lend[i]]];
    k_t[J_l[Lstart[i] : Lend[i]]] = k_t[J_l[Lstart[i - 1] : Lend[i - 1]]]
                                    + k_dev[J_l[Lstart[i] : Lend[i]]];
  }
  ypred = log_a[J_i] - S * b + g_t[J_ii] .* X1 + k_t[J_ii] .* X2;
  print(ypred);
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
  z_dev_g ~ std_normal();
  z_dev_k ~ std_normal();
  // variance terms
  sigma ~ gamma(2, 2);
  sigma_g_s ~ gamma(1, 2);
  sigma_k_s ~ gamma(1, 2);
  // Likelihoo
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
  vector[L_tot] gamma_t_m; 
  vector[L_tot] kappa_t_m;
  int pos; // temporary position integer
  for (r in 1 : R) {
    for (t in pos : L_c[r]) {
      gamma_t_m[t] = sum(g_t[J_l[Lstart[t] : Lend[t]]]) / J[r]; //arithmetic mean of g_t at each time-step
      kappa_t_m[t] = sum(k_t[J_l[Lstart[t] : Lend[t]]]) / J[r]; //arithmetic mean of k_t at each time-step
    }
    pos = L_c[r] + 1;
  }
}
