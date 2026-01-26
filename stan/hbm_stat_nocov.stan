// Hierarchical Ricker model + AR1 error term + NO covariates
//
// Parameters:
//   alpha: exchangeable -> Normal(mu_alpha, sigma_alpha)
//      - non-centered parameterization
//      - sigma_alpha shared among series within a group
//   beta:  series-specific
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
    real x1[N];                           // 1st covariate
    real y[N];                            // response
    int priors_only;                      // should likelihood be ignored?
}
parameters {
    real mu_alpha[Na_groups];             // population-level intercept
    real d_alpha[n_series];               // random alpha deviate
    real beta[n_series];                  // slopes for 1st covariate
    real<lower=0> sigmaNC[n_series];      // residual SD (not corrected for AR1)
    real<lower=0> sigma_alpha[Na_groups]; // population-level intercept SD
    real<lower=-1, upper=1> phi;          // autocorrelation parameter
}
transformed parameters {
    real alpha[n_series];                 // series-specific intercept
    real yhat[N];                         // predicted values
    real epsilon[N];                      // residuals
    real<lower=0> sigma[n_series];        // residual SD (corrected for AR1)
    real tmp_epsilon;                     // temporary var to avoid deep copy

    for(i in 1:n_series) {
        alpha[i] = mu_alpha[a_group[i]] + sigma_alpha[a_group[i]] * d_alpha[i];

        // first data point in series
        yhat[y_start[i]] = alpha[i] + beta[i] * x1[y_start[i]];
        epsilon[y_start[i]] = y[y_start[i]] - yhat[y_start[i]];

        for(t in (y_start[i]+1):y_end[i]) {
            tmp_epsilon = epsilon[t-1];
            yhat[t] = alpha[i] + beta[i] * x1[t] + phi * tmp_epsilon;
            epsilon[t] = y[t] - (yhat[t] - (phi * tmp_epsilon));
        }

        // AR1 impacts variance: sigma^2 = sigmaNC^2 * (1 - phi^2)
        sigma[i] = sqrt(sigmaNC[i] * sigmaNC[i] * (1 - phi * phi));
    }
}
model {
    // priors: population-level
    mu_alpha ~ normal(0, 5);
    sigma_alpha ~ cauchy(0, 2.5);
    phi ~ normal(0, 1);

    // priors: series-specific
    d_alpha ~ normal(0, 1);
    beta ~ normal(0, 5);
    sigmaNC ~ cauchy(0, 2.5);

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
    real yresid[N];                       // realized residuals
    for(i in 1:n_series) {
        for(t in y_start[i]:y_end[i]) {
            log_lik[t] = normal_lpdf(y[t] | yhat[t], sigma[i]);
            yrep[t] = normal_rng(yhat[t], sigma[i]);
            yresid[t] = y[t] - yhat[t];
        }
    }
}

