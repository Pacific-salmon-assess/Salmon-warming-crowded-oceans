// Hierarchical Ricker model + AR1 error term + time-varying covariate
//
// Parameters:
//   alpha: exchangeable -> Normal(mu_alpha, sigma_alpha)
//      - non-centered parameterization
//   beta:  series-specific
//   gamma: time-varying  g[t] = g[t-1] + N(0, sigma_gamma)
//      - non-centered parameteriztion
//      - gamma's are series-specific
//      - sigma_gamma is series specific
//   sigma: series-specific
//   phi:   shared

data {
    int<lower=0> N;                       // total number of years
    int<lower=0> n_series;                // number of time series
    int<lower=0> y_start[n_series];       // ragged start points for series
    int<lower=0> y_end[n_series];         // ragged end points for series
    int<lower=0> Na_groups;               // number of alpha groups
    int<lower=0> a_group[n_series];       // alpha grouping factor
    real x1[N];                           // 1st covariate
    real x2[N];                           // 2nd covariate
    real y[N];                            // response
    int priors_only;                      // should likelihood be ignored?
    real sigma_gamma_df;                  // prior df for sigma_gamma
    real sigma_gamma_mu;                  // prior mu for sigma_gamma
    real sigma_gamma_sd;                  // prior sd for sigma_gamma
}
parameters {
    real mu_alpha[Na_groups];             // population-level intercept
    real<lower=0> sigma_alpha[Na_groups]; // population-level intercept SD
    real beta[n_series];                  // slopes for 1st covariate
    real<lower=-1, upper=1> phi;          // autocorrelation parameter
    real<lower=0> sigmaNC[n_series];      // residual SD (not corrected for AR1)
    real g0[n_series];                    // gamma at t = 1
    real d_alpha[n_series];               // random alpha deviate
    real d_gamma[N];                      // random gamma deviate
    real<lower=0> sigma_gamma[n_series];  // gamma SD
}
transformed parameters {
    real alpha[n_series];                 // series-specific intercept
    real yhat[N];                         // predicted values
    real epsilon[N];                      // residuals
    real<lower=0> sigma[n_series];        // residual SD (corrected for AR1)
    real gamma[N];                        // 2nd covariate effect
    real tmp_epsilon;                     // temporary var to avoid deep copy
    real tmp_gamma;                       // temporary var to avoid deep copy

    for(i in 1:n_series) {
        alpha[i] = mu_alpha[a_group[i]] + sigma_alpha[a_group[i]] * d_alpha[i];

        // first data point in series
        gamma[y_start[i]] = g0[i];
        yhat[y_start[i]] = alpha[i] + beta[i] * x1[y_start[i]] + gamma[y_start[i]] * x2[y_start[i]];
        epsilon[y_start[i]] = y[y_start[i]] - yhat[y_start[i]];

        for(t in (y_start[i]+1):y_end[i]) {
            tmp_gamma = gamma[t-1];
            tmp_epsilon = epsilon[t-1];
            gamma[t] = tmp_gamma + sigma_gamma[i] * d_gamma[t-1];
            yhat[t] = alpha[i] + beta[i] * x1[t] + gamma[t] * x2[t] + phi * tmp_epsilon;
            epsilon[t] = y[t] - (yhat[t] - (phi * tmp_epsilon));
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

    // priors: series-specific
    sigma_gamma ~ student_t(sigma_gamma_df, sigma_gamma_mu, sigma_gamma_sd);
    sigmaNC ~ student_t(3, 0, 3);
    beta ~ normal(0, 5);
    g0 ~ normal(0, 3);
    d_alpha ~ normal(0, 1);
    d_gamma ~ normal(0, 1);

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
    real signal_noise[n_series];          // signal-to-noise ratio
    for(i in 1:n_series) {
        for(t in y_start[i]:y_end[i]) {
            log_lik[t] = normal_lpdf(y[t] | yhat[t], sigma[i]);
            yrep[t] = normal_rng(yhat[t], sigma[i]);
        }
    }

    for(i in 1:n_series) {
        signal_noise[i] = sigma_gamma[i] / sigma[i];
    }
}

