// generated with brms 2.6.0
functions { 
  /* cumulative-logit log-PDF for a single response
  * Args:
    *   y: response category
  *   mu: linear predictor
  *   thres: ordinal thresholds
  *   disc: discrimination parameter
  * Returns:
    *   a scalar to be added to the log posterior
  */
    real cumulative_logit_lpmf(int y, real mu, vector thres, real disc) {
      int ncat = num_elements(thres) + 1;
      real p;
      if (y == 1) {
        p = inv_logit(disc * (thres[1] - mu));
      } else if (y == ncat) {
        p = 1 - inv_logit(disc * (thres[ncat - 1] - mu));
      } else {
        p = inv_logit(disc * (thres[y] - mu)) -
          inv_logit(disc * (thres[y - 1] - mu));
      }
      return log(p);
    } 
  /* cumulative-logit log-PDF for a single response
  * including category specific effects
  * Args:
    *   y: response category
  *   mu: linear predictor
  *   mucs: predictor for category specific effects
  *   thres: ordinal thresholds
  *   disc: discrimination parameter
  * Returns:
    *   a scalar to be added to the log posterior
  */
    real cumulative_logit_cs_lpmf(int y, real mu, row_vector mucs, vector thres, real disc) {
      int ncat = num_elements(thres) + 1;
      real p;
      if (y == 1) {
        p = inv_logit(disc * (thres[1] - mucs[1] - mu));
      } else if (y == ncat) {
        p = 1 - inv_logit(disc * (thres[ncat - 1] - mucs[ncat - 1] - mu));
      } else {
        p = inv_logit(disc * (thres[y] - mucs[y] - mu)) -
          inv_logit(disc * (thres[y - 1] - mucs[y - 1] - mu));
      }
      return log(p);
    } 
} 

data { 
  int N;  // total number of observations 
  int Y[N];  // response variable 
  int<lower=2> ncat;  // number of categories 
  int<lower=1> K;  // number of population-level effects 
  matrix[N, K]   X;  // population-level design matrix (our design matrix was defined without intercepts, only skills we have )
} 
transformed data { 
  int Kc = K; 
  matrix[N, K] Xc;  // 
  Xc=X;
} 

parameters { 
  vector[Kc] beta;  // population-level effects 
  real first_temp_Intercept;        // fake first thresholds
  vector<lower=0>[ncat-2] delta;   // delta parameters in the threshold prior;
}
transformed parameters{
  ordered[ncat-1] temp_Intercept;  // temporary thresholds
   temp_Intercept[1] = first_temp_Intercept;
  for (k in 2:(ncat-1)){
    temp_Intercept[k] = temp_Intercept[k-1] + delta[k-1]; // threshold transformation
  }
}

model { 
  vector[N] mu = Xc * beta;
  
  // priors including all constants 
  target += normal_lpdf(first_temp_Intercept | 0, 10);  // first threshold prior 
  for (k in 1:(ncat-2)){
    target+= lognormal_lpdf(delta[k]|0, 10);         // delta prior
    }
  target+=normal_lpdf(beta|0,10);
  // likelihood including all constants 
  
  for (n in 1:N) {
    target += ordered_logistic_lpmf(Y[n] | mu[n], temp_Intercept);
  }
  
} 

generated quantities {

  vector[N] log_lik;
  real dev;
  int y_pred[N];
  vector[N] mu = Xc * beta;
  dev=0;
  for (i in 1:N) {
    log_lik[i] = ordered_logistic_lpmf(Y[i] | mu[i], temp_Intercept);
    dev=dev-2*log_lik[i];
    y_pred[i]=ordered_logistic_rng(mu[i], temp_Intercept);
  }



}
