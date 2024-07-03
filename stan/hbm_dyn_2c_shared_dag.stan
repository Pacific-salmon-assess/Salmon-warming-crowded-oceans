// Hierarchical Ricker model + AR1 error term + time-varying covariate (group)
//
// Parameters:
//   alpha: exchangeable -> Normal(mu_alpha, sigma_alpha)
//      - non-centered parameterization
//   beta:  series-specific
//   gamma: time-varying  g[t] = g[t-1] + N(0, sigma_gamma)
//      - non-centered parameteriztion
//      - gamma's shared among series in a group (same as kappa groups)
//   kappa: time-varying  k[t] = k[t-1] + N(0, sigma_kappa)
//      - non-centered parameteriztion
//      - kappa's shared among series in a group (same as gamma groups)
//   sigma: series-specific
//   phi:   shared

data {
    int<lower=0> N;                      // total number of years
    int<lower=0> n_series;               // number of time series
    int<lower=0> Ng_groups;              // number of gamma groups
    int<lower=0> Na_groups;              // number of alpha groups
    int<lower=0> a_group[n_series];      // alpha grouping factor
    int<lower=0> g_group[n_series];      // gamma/kappa grouping factor
    int<lower=0> year[N];                // year of response (maps gamma -> y)
    int<lower=0> y_start[n_series];      // ragged start points for series
    int<lower=0> y_end[n_series];        // ragged end points for series
    int<lower=0> g_start[Ng_groups];     // ragged start points for gamma/kappa
    int<lower=0> g_end[Ng_groups];       // ragged end points for gamma/kappa
    int<lower=0> Ng;                     // total number of gammas/kappas
    real x1[N];                          // 1st covariate
    real x2[N];                          // 2nd covariate
    real x3[N];                          // 3rd covariate
    real y[N];                           // response
    int priors_only;                     // should likelihood be ignored?
}
parameters {
    real mu_alpha[Na_groups];             // population-level intercept
    real<lower=0> sigma_alpha[Na_groups]; // population-level intercept SD
    real beta[n_series];                  // slopes for 1st covariate
    real<lower=-1, upper=1> phi;          // autocorrelation parameter
    real<lower=0> sigmaNC[n_series];      // residual SD (not corrected for AR1)
    real g0[Ng_groups];                   // gamma at t = 1
    real k0[Ng_groups];                   // kappa at t=1  
    real d_alpha[n_series];               // random alpha deviate
    matrix[Ng_groups,Ng-1] d_gamma;         // random gamma deviate matrix
    matrix[Ng_groups,Ng-1] d_kappa;         // random kappa deviate matrix
    real<lower=0> sigma_gamma[Ng_groups]; // time-varying gamma SD
    real<lower=0> sigma_kappa[Ng_groups]; // time-varying kappa SD
    vector[n_series] d_gamma_i;	                // z-score for stock-specific time-invariant gamma
    vector<lower=0>[Ng_groups] sigma_gamma_i;   //stock-specific gamma SD (indexed by group)
    vector[n_series] d_kappa_i;	                // z-score for stock-specific time-invariant kappa
    vector<lower=0>[Ng_groups] sigma_kappa_i;   //stock-specific gamma SD (indexed by group)
}
transformed parameters {
    real alpha[n_series];                 // series-specific intercept
    real yhat[N];                         // predicted values
    real epsilon[N];                      // residuals
    real<lower=0> sigma[n_series];        // residual SD (corrected for AR1)
    matrix[Ng_groups,Ng] gamma;                       // 2nd covariate effect - stratified by gamma group
    matrix[Ng_groups,Ng] kappa;                       // 3rd covariate effect - stratified by gamma group
    real tmp_epsilon;                     // temporary var to avoid deep copy
    real tmp_gamma;                       // temporary var to avoid deep copy
    real tmp_kappa;                       // temporary var to avoid deep copy
    vector[n_series] gamma_i;                 // series-specific time-invariant covariate
    vector[n_series] kappa_i;                 // series-specific time-invariant covariate
    
	gamma_i = sigma_gamma_i[a_group].*d_gamma_i;
	kappa_i = sigma_kappa_i[a_group].*d_kappa_i;
	
	for(j in 1:Ng_groups) {
        gamma[j,1] = g0[j];
        kappa[j,1] = k0[j];
        for(t in 2:n_years) {
            tmp_gamma = gamma[j,t-1];
            gamma[j,t] = tmp_gamma + sigma_gamma[j] * d_gamma[j,t-1];
            tmp_kappa = kappa[j,t-1];
            kappa[j,t] = tmp_kappa + sigma_kappa[j] * d_kappa[j,t-1];
        }
    }

    for(i in 1:n_series) {
        alpha[i] = mu_alpha[a_group[i]] + sigma_alpha[a_group[i]] * d_alpha[i];

        // first data point in series
        yhat[y_start[i]] = alpha[i] + beta[i] * x1[y_start[i]]+ gamma_i[i]*x2[y_start[i]] + gamma[g_group[i],year[y_start[i]]] * x2[y_start[i]]+ kappa_i[i]*x3[y_start[i]] +  kappa[g_group[i],year[y_start[i]]] * x3[y_start[i]];
        epsilon[y_start[i]] = y[y_start[i]] - yhat[y_start[i]];

        for(t in (y_start[i]+1):y_end[i]) {
            tmp_epsilon = epsilon[t-1];

            yhat[t] = alpha[i] + beta[i] * x1[t]+ gamma_i[i]*x2[t] + gamma[year[t]] * x2[t]+ kappa_i[i]*x3[t] + kappa[year[t]] * x3[t] + phi * tmp_epsilon;
            epsilon[t] = y[t] - (yhat[t] - (phi * tmp_epsilon));

            yhat[t] = alpha[i] + beta[i] * x1[t]+ gamma_i[i]*x2[t] + gamma[g_group[i],year[t]] * x2[t]+ kappa_i[i]*x3[t] + kappa[g_group[i],year[t]] * x3[t] + (phi^(year[t]-year[t-1]))*tmp_epsilon;
            epsilon[t] = y[t] - (yhat[t] - ((phi^(year[t]-year[t-1])) * tmp_epsilon));
        }

        // AR1 impacts variance: sigma^2 = sigmaNC^2 * (1 - phi^2)
        sigma[i] = sqrt(sigmaNC[i] * sigmaNC[i] * (1 - phi * phi));
    }
}
model {
    // priors: population-level
    mu_alpha ~ normal(1, 5);
    sigma_alpha ~ student_t(3, 0, 3);
    phi ~ normal(0, 1);
    g0 ~ normal(0, 3);
    to_vector(d_gamma) ~ normal(0, 1);
    to_vector(d_kappa) ~ normal(0,1);
	d_gamma_i ~ normal(0, 1);
	d_kappa_i ~ normal(0, 1);
	
    sigma_gamma ~ student_t(3, 0, 0.5);
    sigma_kappa ~ student_t(3, 0, 0.5);
	 sigma_kappa_i ~ student_t(3, 0, 0.5);
	  sigma_gamma_i ~ student_t(3, 0, 0.5);

    // priors: series-specific
    sigmaNC ~ student_t(3, 0, 3);
    beta ~ normal(0, 5);
    d_alpha ~ normal(0, 1);

    // likelihood
    if(!priors_only) {
        for(i in 1:n_series) {
            y[y_start[i]:y_end[i]] ~ normal(yhat[y_start[i]:y_end[i]], sigma[i]);
        }
    }
}
generated quantities {
    real log_lik[N];                      // log-likelihood
    real yrep[N];                         // replicated datasets
    real signal_noise_g[n_series];          // signal-to-noise ratio
    real signal_noise_k[n_series];
    for(i in 1:n_series) {
        for(t in y_start[i]:y_end[i]) {
            log_lik[t] = normal_lpdf(y[t] | yhat[t], sigma[i]);
            yrep[t] = normal_rng(yhat[t], sigma[i]);
        }
    }

    for(j in 1:Ng_groups) {
        for(i in 1:n_series) {
            signal_noise_g[i] = sigma_gamma[j] / sigma[i];
            signal_noise_k[i] = sigma_kappa[j] / sigma[i];
        }
    }
}

