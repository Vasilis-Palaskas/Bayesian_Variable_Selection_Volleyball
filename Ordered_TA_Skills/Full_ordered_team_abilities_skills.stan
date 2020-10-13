data { 
  int N;  // total number of matches
  int Y[N];  // response variable 
  int<lower=2> ncat;  // number of response categories 
  int<lower=1> K;  // number of candidate variables
  matrix[N, K]   X;  // population-level design matrix (our design matrix was defined without intercepts, only skills we have )
  int  home_team[N]; // vector of home ground teams in the whole league
  int  away_team[N]; // vector of away ground teams in the whole league
} 


parameters { 
  vector[K] beta;  // parameters of candidate variables 
  ordered[ncat-1] temp_Intercept;  // temporary thresholds 
  real gen_abil_raw[11]; // general ability parameters (12 teams in total)
} 

transformed parameters { 
  // sum to zero constraint parameterisation for the general ability parameters
  vector[12]   gen_abil; 
  for (t in 1:(12-1)) {
    gen_abil[t] = gen_abil_raw[t];
  }
  gen_abil[12] = -sum(gen_abil_raw);
}

model { 
  vector[N] mu = X * beta+gen_abil[home_team]-gen_abil[away_team];
  
  
  // priors 
  target += normal_lpdf(temp_Intercept | 0, 10); 
  target += normal_lpdf(gen_abil |0,10); 
  target+=normal_lpdf(beta|0,10);
  // likelihood
  
  for (n in 1:N) {
    target += ordered_logistic_lpmf(Y[n] | mu[n], temp_Intercept);
  }
  
} 
