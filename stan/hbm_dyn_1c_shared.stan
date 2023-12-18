// Hierarchical Ricker model + AR1 error term + time-varying covariate (group)
//
// Parameters:
//   alpha: exchangeable -> Normal(mu_alpha, sigma_alpha)
//      - non-centered parameterization
//   beta:  series-specific
//   gamma: time-varying  g[t] = g[t-1] + N(0, sigma_gamma)
//      - non-centered parameteriztion
//      - gamma's shared among series in a group
//   sigma: series-specific
//   phi:   shared

data {
    int<lower=0> N;                      // total number of years
    int<lower=0> n_series;               // number of time series
    int<lower=0> Ng_groups;              // number of gamma groups
    int<lower=0> Na_groups;              // number of alpha groups
    int<lower=0> a_group[n_series];      // alpha grouping factor
    int<lower=0> g_group[n_series];      // gamma grouping factor
    int<lower=0> year[N];                // year of response (maps gamma -> y)
    int<lower=0> y_start[n_series];      // ragged start points for series
    int<lower=0> y_end[n_series];        // ragged end points for series
    int<lower=0> g_start[Ng_groups];     // ragged start points for gamma
    int<lower=0> g_end[Ng_groups];       // ragged end points for gamma
    int<lower=0> Ng;                     // total number of gammas
    real x1[N];                          // 1st covariate
    real x2[N];                          // 2nd covariate
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
    real d_alpha[n_series];               // random alpha deviate
    real d_gamma[Ng-1];                   // random gamma deviate
    real<lower=0> sigma_gamma[Ng_groups]; // gamma SD
}
transformed parameters {
    real alpha[n_series];                 // series-specific intercept
    real yhat[N];                         // predicted values
    real epsilon[N];                      // residuals
    real<lower=0> sigma[n_series];        // residual SD (corrected for AR1)
    real gamma[Ng];                       // 2nd covariate effect
    real tmp_epsilon;                     // temporary var to avoid deep copy
    real tmp_gamma;                       // temporary var to avoid deep copy

    for(j in 1:Ng_groups) {
        gamma[g_start[j]] = g0[j];
        for(t in (g_start[j]+1):g_end[j]) {
            tmp_gamma = gamma[t-1];
            gamma[t] = tmp_gamma + sigma_gamma[j] * d_gamma[t-1];
        }
    }

    for(i in 1:n_series) {
        alpha[i] = mu_alpha[a_group[i]] + sigma_alpha[a_group[i]] * d_alpha[i];

        // first data point in series
        yhat[y_start[i]] = alpha[i] + beta[i] * x1[y_start[i]] + gamma[year[y_start[i]]] * x2[y_start[i]];
        epsilon[y_start[i]] = y[y_start[i]] - yhat[y_start[i]];

        for(t in (y_start[i]+1):y_end[i]) {
            tmp_epsilon = epsilon[t-1];
            yhat[t] = alpha[i] + beta[i] * x1[t] + gamma[year[t]] * x2[t] + phi * tmp_epsilon;
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
    g0 ~ normal(0, 3);
    d_gamma ~ normal(0, 1);
    sigma_gamma ~ student_t(3, 0, 0.5);

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
    real signal_noise[n_series];          // signal-to-noise ratio
    for(i in 1:n_series) {
        for(t in y_start[i]:y_end[i]) {
            log_lik[t] = normal_lpdf(y[t] | yhat[t], sigma[i]);
            yrep[t] = normal_rng(yhat[t], sigma[i]);
        }
    }

    for(j in 1:Ng_groups) {
        for(i in 1:n_series) {
            signal_noise[i] = sigma_gamma[j] / sigma[i];
        }
    }
}

