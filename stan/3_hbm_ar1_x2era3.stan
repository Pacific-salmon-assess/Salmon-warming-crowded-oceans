// Hierarchical Ricker model + AR1 error term + covariate (group)
//
// Parameters:
//   alpha: exchangeable -> Normal(mu_alpha, sigma_alpha)
//      - non-centered parameterization
//   beta:  series-specific
//   gamma: exchangeable within a group and era -> Normal(mu_gamma, sigma_gamma)
//      - non-centered parameterization
//      - mu_gamma shared among series within a group and era
//      - sigma_gamma shared among series within a group (shared arcross eras)
//      - three eras -> early (mu_gamma2); middle (mu_gamma2); late (mu_gamma3)
//   sigma: series-specific
//   phi:   shared

data {
    int<lower=0> N;                       // total number of years
    int<lower=0> n_series;                // number of time series
    int<lower=0> Na_groups;               // number of alpha groups
    int<lower=0> Ng_groups;               // number of gamma groups
    int<lower=0> a_group[n_series];       // alpha grouping factor
    int<lower=0> g_group[n_series];       // gamma grouping factor
    int<lower=0> y_start[n_series];       // ragged start points for series
    int<lower=0> y_end[n_series];         // ragged end points for series
    int<lower=0, upper=1> era1[N];        // 1st era logical
    int<lower=0, upper=1> era2[N];        // 2nd era logical
    int<lower=0, upper=1> era3[N];        // 3rd era logical
    real x1[N];                           // 1st covariate
    real x2[N];                           // 2nd covariate
    real y[N];                            // response
    int priors_only;                      // should likelihood be ignored?
}
parameters {
    real mu_alpha[Na_groups];             // population-level intercept
    real mu_gamma1[Ng_groups];            // population-level x2 mean -- early era
    real mu_gamma2[Ng_groups];            // population-level x2 mean -- middle era
    real mu_gamma3[Ng_groups];            // population-level x2 mean -- late era
    real d_alpha[n_series];               // random alpha deviate
    real d_gamma1[n_series];              // random gamma1 deviate
    real d_gamma2[n_series];              // random gamma2 deviate
    real d_gamma3[n_series];              // random gamma3 deviate
    real beta[n_series];                  // slopes for 1st covariate
    real<lower=0> sigmaNC[n_series];      // residual SD (not corrected for AR1)
    real<lower=0> sigma_alpha[Na_groups]; // population-level intercept SD
    real<lower=-1, upper=1> phi;          // autocorrelation parameter
    real<lower=0> sigma_gamma[Ng_groups]; // population-level x2 SD
}
transformed parameters {
    real alpha[n_series];                 // series-specific intercept
    real gamma1[n_series];                // series-specific x2 effect (early)
    real gamma2[n_series];                // series-specific x2 effect (middle)
    real gamma3[n_series];                // series-specific x2 effect (late)
    real yhat[N];                         // predicted values
    real epsilon[N];                      // residuals
    real<lower=0> sigma[n_series];        // residual SD (corrected for AR1)
    real tmp_epsilon;                     // temporary var to avoid deep copy

    for(i in 1:n_series) {
        alpha[i] = mu_alpha[a_group[i]] + sigma_alpha[a_group[i]] * d_alpha[i];

        for(t in (y_start[i]):y_end[i]) {
            gamma1[i] = (mu_gamma1[g_group[i]] + sigma_gamma[g_group[i]] * d_gamma1[i]);
            gamma2[i] = (mu_gamma2[g_group[i]] + sigma_gamma[g_group[i]] * d_gamma2[i]);
            gamma3[i] = (mu_gamma3[g_group[i]] + sigma_gamma[g_group[i]] * d_gamma3[i]);
            if(t == y_start[i]) { // first data point in series
                yhat[y_start[i]] = alpha[i] + beta[i] * x1[y_start[i]] +
                    (gamma1[i] * x2[y_start[i]] * era1[y_start[i]]) +
                    (gamma2[i] * x2[y_start[i]] * era2[y_start[i]]) +
                    (gamma3[i] * x2[y_start[i]] * era3[y_start[i]]);
                epsilon[y_start[i]] = y[y_start[i]] - yhat[y_start[i]];
            } else {
                tmp_epsilon = epsilon[t-1];
                yhat[t] = alpha[i] + beta[i] * x1[t] +
                    (gamma1[i] * x2[t] * era1[t]) +
                    (gamma2[i] * x2[t] * era2[t]) +
                    (gamma3[i] * x2[t] * era3[t]) +
                    phi * tmp_epsilon;
                epsilon[t] = y[t] - (yhat[t] - (phi * tmp_epsilon));
            }
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

    // priors: group-specific
    mu_gamma1 ~ normal(0, 5);
    mu_gamma2 ~ normal(0, 5);
    mu_gamma3 ~ normal(0, 5);
    sigma_gamma ~ student_t(3, 0, 3);

    // priors: series-specific
    d_alpha ~ normal(0, 1);
    d_gamma1 ~ normal(0, 1);
    d_gamma2 ~ normal(0, 1);
    d_gamma3 ~ normal(0, 1);
    beta ~ normal(0, 5);
    sigmaNC ~ student_t(3, 0, 3);

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
    for(i in 1:n_series) {
        for(t in y_start[i]:y_end[i]) {
            log_lik[t] = normal_lpdf(y[t] | yhat[t], sigma[i]);
            yrep[t] = normal_rng(yhat[t], sigma[i]);
        }
    }
}

