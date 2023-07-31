functions {
  vector normalize(vector x) {
  return x / sum(x);
}
}
data {
  int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[N] X1; // vector for covariate 1
  matrix[K,K] alpha_dirichlet; //prior inputs for dirichlet 
 }
parameters {
  // Discrete state model
  simplex[K] A[K]; // transition probabilities
  simplex[K] pi1; // initial state probabilities

  // A[i][j] = p(z_t = j | z_{t-1} = i)
  // Continuous observation model
  real<lower=0> log_a; // max. productivity
  real log_b; // rate capacity - fixed in this
  vector<lower=0>[K] sigma; // observation standard deviations
  ordered[K] beta1; //covariate 1 - K states
  
  vector<lower = -1, upper = 1>[K] rho; //autocorrelation in prod. residuals
}

transformed parameters {
  vector[K] logalpha[N];
  real b; //
  vector[K] mu[N]; //expectation for each state at each time
  vector[K] epsilon[N]; //residuals
  vector<lower=0>[K] sigmaAR;

  b=exp(log_b);

 
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
  real accumulator1[K];
  mu[1][1] = log_a - b*S[1] + beta1[1]*X1[1]; //initial expectation - state 1
  mu[1][2] = log_a - b*S[1] + beta1[2]*X1[1]; //initial expectation - state 2
  epsilon[1][1] = R_S[1] - mu[1][1]; //initial residual - state 1
  epsilon[1][2] = R_S[1] - mu[1][2]; //initial residual - state 2
  
  //first-step of log likelihoods
  for(j in 1:2){logalpha[1][j] = log(pi1[j]) + normal_lpdf(R_S[1] |mu[1][j], sigma[j]);
  
   sigmaAR[j] = sigma[j]*sqrt(1-rho[j]^2);
  }
  
  //for every other observation:
  for (t in 2:N) {
  for (j in 1:K) { // j = current (t)
  for (i in 1:K) { // i = previous (t-1)
  mu[t][j]= log_a - b*S[t] + beta1[j]*X1[t]+(rho[j]*epsilon[t-1][j]);
  epsilon[t][j]= R_S[t] - mu[t][j];
  
  // Murphy (2012) p. 609 eq. 17.48
  // belief state + transition prob + local evidence at t
  accumulator1[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |mu[t][j], sigmaAR[j]);
  }
  logalpha[t, j] = log_sum_exp(accumulator1);
  }
  }
  } // Forward
  
   
}
model{
  log_a ~ normal(1.5,2.5);
  log_b ~ normal(-12,3);
  beta1 ~ normal(0,1); //set priors for both sets at 1 std
  sigma ~ normal(0,1); //half normal on variance (lower limit of zero)
  
  //autocorrelation term
  rho ~ uniform(-1,1);
  
  //initial state probs.
  pi1 ~ dirichlet(rep_vector(1,K));
  
  //transition matrix
  for(k in 1:K){
  A[k,] ~ dirichlet(alpha_dirichlet[k,]);
  }
  
  target += log_sum_exp(logalpha[N]);
}
generated quantities {
  vector[N] log_lik;
  int<lower=1, upper=K> zstar[N];
  real logp_zstar;
  vector[K] alpha[N];
  vector[K] logbeta[N];
  vector[K] loggamma[N];
  vector[K] beta[N];
  vector[K] gamma[N];
  
  real S_max; 
  
  { // Forward algortihm
  for (t in 1:N)
  alpha[t] = softmax(logalpha[t]);
  } // Forward
  
  { // Backward algorithm log p(y_{t+1:T} | z_t = j)
  real accumulator2[K];
  for (j in 1:K)
  logbeta[N, j] = 1;
  for (tforward in 0:(N-2)) {
  int t;
  t = N - tforward;
  for (j in 1:K) { // j = previous (t-1)
  for (i in 1:K) { // i = next (t)
  // Murphy (2012) Eq. 17.58
  // backwards t + transition prob + local evidence at t
  accumulator2[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] | mu[t][i], sigmaAR[i]);
  }
  logbeta[t-1, j] = log_sum_exp(accumulator2);
  }
  }
  for (t in 1:N)
  beta[t] = softmax(logbeta[t]);
  } // Backward

  { // Forward-backward algorithm log p(z_t = j | y_{1:N})
  for(t in 1:N) {
  loggamma[t] = alpha[t] .* beta[t];
  }
  for(t in 1:N)
  gamma[t] = normalize(loggamma[t]);
  } // Forward-backward
  
  { // Viterbi algorithm
  int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
  real delta[N, K]; // max prob for the sequence up to t
  // that ends with an emission from state k
  for (j in 1:K)
  delta[1, K] = normal_lpdf(R_S[1] | mu[1][j], sigma);
  for (t in 2:N) {
    for (j in 1:K) { // j = current (t)
      delta[t, j] = negative_infinity();
      for (i in 1:K) { // i = previous (t-1)
        real logp;
        logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | mu[t][j], sigmaAR[j]);
        if (logp > delta[t, j]) {
          bpointer[t, j] = i;
          delta[t, j] = logp;
        }
      }
    }
  }
  logp_zstar = max(delta[N]);
  for (j in 1:K)
    if (delta[N, j] == logp_zstar)
      zstar[N] = j;
  for (t in 1:(N - 1)) {
    zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
  }
  }

log_lik[1] = normal_lpdf(R_S[1]|mu[1][zstar[1]], sigma);
for(n in 2:N)log_lik[n] = normal_lpdf(R_S[n]|mu[n][zstar[n]], sigmaAR);
   
S_max = 1/b;

}
