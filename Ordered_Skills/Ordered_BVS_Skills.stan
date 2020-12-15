

data { 
  int N;  // total number of matches
  int Y[N];  // response variable 
  int<lower=2> ncat;  // number of response categories 
  int<lower=1> K;  // number of candidate variables
  matrix[N, K]   X;  // design matrix (our design matrix was defined without intercepts, only skills we have )
  vector[K] gammas; // binary indicators of candidate variables
  vector[K] post_mean_betas; //posterior means of betas from a MCMC pilot run of full model
  vector[K] post_sd_betas; //posterior stan. deviat. of betas from a MCMC pilot run of full model
  } 

parameters { 
  vector[K] betas;  // parameters of candidate variables 
  real first_temp_Intercept;        // fake first thresholds
  vector<lower=0>[ncat-2] delta;   // delta parameters in the threshold prior;
} 



transformed parameters { 
  ordered[ncat-1] temp_Intercept;  // temporary thresholds
   vector[K] gb;
   
  temp_Intercept[1] = first_temp_Intercept;
  for (k in 2:(ncat-1)){
    temp_Intercept[k] = temp_Intercept[k-1] + delta[k-1]; // threshold transformation
  }
   
  for (j in 1:K){
    gb[j]=gammas[j]*betas[j];
  }
  // sum to zero constraint for the general ability parameters
}

model { 
  // Linear predictor
  vector[N] mu = X * gb;
  
  // Priors of all parameters except for the candidate variables' parameters
  target += normal_lpdf(first_temp_Intercept | 0, 10); // first threshold prior 
  for (k in 1:(ncat-2)){
    target+= lognormal_lpdf(delta[k]|0, 10);           // delta prior
    } 

  
  // Prior+Pseudoprior of the candidate variables' parameters
  for (j in 1:K){
   target+=normal_lpdf(betas[j]|gammas[j]*0+(1-gammas[j])*post_mean_betas[j],(sqrt(N)^(gammas[j]))*post_sd_betas[j]);//mixture normal in mean
  }
    
  // Likelihood
  
  for (i in 1:N) {
    target += ordered_logistic_lpmf(Y[i] | mu[i], temp_Intercept);
  }
} 

generated quantities{
  matrix[N,K] log_lik_one; // matrix with log likelihoods when gamma[j]=1
  matrix[N,K] log_lik_zero; // matrix with log likelihoods when gamma[j]=0
  
  vector[K] gb_new; // the gammas*betas when we want to estimate the log likelihood
  // for the opposite value of binary indicator within the if-loop.
  
  vector[N] mu_new;// the linear predictor when we want to estimate the 
  // log likelihood for the opposite value of binary indicator within the if-loop.
  
  vector[N] mu = X * gb;
   for (i in 1:N) {
     for (j in 1:K){
      // In each case of if-loop, we estimate the log likelihood for both cases
      // of binary indicators gamma
      if (gb[j]==0){
         log_lik_zero[i,j] = ordered_logistic_lpmf(Y[i] |mu[i], temp_Intercept);
         gb_new=gb;
         gb_new[j]=1*betas[j];// gamma[j]=1
         mu_new[i]=X[i,] * gb_new;
         log_lik_one[i,j]=  ordered_logistic_lpmf(Y[i]| mu_new[i], temp_Intercept);
      } else {
         log_lik_one[i,j]=ordered_logistic_lpmf(Y[i] |mu[i], temp_Intercept);
         gb_new=gb;
         gb_new[j]=0;
         mu_new[i]=X[i,] * gb_new;
         log_lik_zero[i,j]=  ordered_logistic_lpmf(Y[i]| mu_new[i], temp_Intercept);
      }
    }
  }
}



