functions {
  vector normalize(vector x) {
  return x / sum(x);
}
}
data {
  int<lower=1> N;//number of osbservations
  int<lower=1> L;//Time-series length (ie. time span)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  matrix[K,K] alpha_dirichlet; //prior inputs for dirichlet
  int ii[N]; //index of brood years by observation
  int ll[L]; //index of observations by brood year - missing = 0
 }
parameters {
  // Discrete state model
  simplex[K] A[K]; // transition probabilities
  simplex[K] pi1; // initial state probabilities
  
  // A[i][j] = p(z_t = j | z_{t-1} = i)
  
  // Continuous observation model
  real<lower=0> log_a; // max. productivity
  real log_b; // rate capacity - fixed in this
  real<lower=0> sigma; // observation standard deviations
  
  //state variables
  ordered[K] beta1; //covariate 1 - K states
}

transformed parameters {
  vector[K] logalpha[L];
  real b; //

  b=exp(log_b);

{ // Forward algorithm log p(z_t = j | y_{1:t})
  real accumulator1[K];

  logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1] + beta1*X1[1], sigma);
+ normal_lpdf(R_S[t] |log_a - b*S[t]+ beta1[i]*X1[t], sigma)
  for (t in 2:L) {
  for (j in 1:K) { // j = current (t)
  for (i in 1:K) { // i = previous (t-1)
  // Murphy (2012) p. 609 eq. 17.48
  // belief state + transition prob + local evidence at t
 
  if(ll[t]==0){
  accumulator1[i] = logalpha[t-1, i] + log(A[i, j]);
  }else{
    accumulator1[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[ll[t]] |log_a - b*S[ll[t]]+ beta1[i]*X1[ll[t]], sigma);
}
  }
  logalpha[t, j] = log_sum_exp(accumulator1);
  }
  } // Forward
}
}
model{
  log_a ~ normal(1.5,2.5);
  log_b ~ normal(-12,3);
  
  sigma ~ normal(0,1); //half normal on variance (lower limit of zero)
  
  pi1~ dirichlet(rep_vector(1,K));

  for(k in 1:K){
  A[k,] ~ dirichlet(alpha_dirichlet[k,]);
  }
  
  target += log_sum_exp(logalpha[L]);
}
generated quantities {
  vector[N] log_lik;
  int<lower=1, upper=K> zstar[L];
  real logp_zstar;
  vector[K] alpha[L];
  vector[K] logbeta[L];
  vector[K] loggamma[L];
  vector[K] beta[L];
  vector[K] gamma[L];
  
  real S_max;
  vector[K] U_msy;
  vector[K] S_msy;
  
  
  { // Forward algortihm
  for (t in 1:L)
  alpha[t] = softmax(logalpha[t]);
  } // Forward
  
  { // Backward algorithm log p(y_{t+1:T} | z_t = j)
  real accumulator2[K];
  for (j in 1:K)
  logbeta[L, j] = 1;
  for (tforward in 0:(L-2)) {
  int t;
  t = L - tforward;
  for (j in 1:K) { // j = previous (t-1)
  for (i in 1:K) { // i = next (t)
  // Murphy (2012) Eq. 17.58
  // backwards t + transition prob + local evidence at t
  
  if(ll[t]==0){
  accumulator2[i] = logbeta[t, i] + log(A[j, i]); //if observation missing - propagate transition probabilities
  }
  else{
  accumulator2[i] = logbeta[t, i] + log(A[j, i])+normal_lpdf(R_S[ll[t]] | log_a - b*S[ll[t]]+beta1[i]*X1[ll[t]], sigma); //if observed - add local evidence
  }
  }
  logbeta[t-1, j] = log_sum_exp(accumulator2);
  }
  }
  for (t in 1:L)
  beta[t] = softmax(logbeta[t]);
  } // Backward

  { // Forward-backward algorithm log p(z_t = j | y_{1:N})
  for(t in 1:L) {
  loggamma[t] = alpha[t] .* beta[t];
  }
  for(t in 1:L)
  gamma[t] = normalize(loggamma[t]);
  } // Forward-backward
  
  { // Viterbi algorithm
  int bpointer[L, K]; // backpointer to the most likely previous state on the most probable path
  real delta[L, K]; // max prob for the sequence up to t
  // that ends with an emission from state k
  for (j in 1:K)
  delta[1, K] = normal_lpdf(R_S[1] | log_a - b*S[1]+beta1[j]*X1[1], sigma);
  for (t in 2:L) {
    for (j in 1:K) { // j = current (t)
      delta[t, j] = negative_infinity();
      for (i in 1:K) { // i = previous (t-1)
        real logp;
		
	if(ll[t]==0){
	logp = delta[t-1, i] + log(A[i, j]);
	}
	else{
        logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[ll[t]] | log_a - b*S[ll[t]]+beta1[j]*X1[ll[t]], sigma);
	}
	
        if (logp > delta[t, j]) {
          bpointer[t, j] = i;
          delta[t, j] = logp;
        }
      }
    }
  }
  logp_zstar = max(delta[L]);
  for (j in 1:K)
    if (delta[L, j] == logp_zstar)
      zstar[L] = j;
  for (t in 1:(L - 1)) {
    zstar[L - t] = bpointer[L - t + 1, zstar[L - t + 1]];
  }
  }


for(n in 1:N) log_lik[n] = normal_lpdf(R_S[ii[n]]|log_a - S[n]*b+beta1[zstar[ii[n]]]*X1[ii[n]], sigma);
normal_lpdf(R_S[n]|log_a[zstar[ii[n]]] - S[n]*b, sigma);
   
S_max = 1/b;

}