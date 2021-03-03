// Define the pmf of zdts distribution
functions {
  
  real skellam_lpmf(int k, real mu1, real mu2) {
    real total;
    real log_prob;
    
    total = (- mu1 - mu2) + (log(mu1) - log(mu2)) * k / 2;
    log_prob = total + log(modified_bessel_first_kind(abs(k), 2 * sqrt(mu1*mu2)));
    
    return log_prob;
  }
  
  real skellam_without_lpmf(int k, real mu1, real mu2) {
    real log_prob_new;
    vector[6] lpmfs;
    real normalization;
    
    for(i in 1:3) {
      lpmfs[i] = skellam_lpmf(i - 4 | mu1, mu2);
      lpmfs[i + 3] = skellam_lpmf(i | mu1, mu2);
    }
    normalization = log_sum_exp(lpmfs);
    if (k > 3 )  {
      log_prob_new=-700;
      return log_prob_new;
    } else if (k <(-3)){
      log_prob_new=-700;
      return log_prob_new;
    } else if (k ==0){
      log_prob_new=-700;
      return log_prob_new;
    }  else {
      log_prob_new=skellam_lpmf(k|mu1,mu2)-normalization;
      return log_prob_new;
    }
    
  }
}


data {
  int <lower=1> n_games; // total number of matches
  int <lower=1> n_teams; // number of teams 12
  int<lower=1> K;       // number of candidate variables
  matrix[n_games, K] X_home;   // design matrix for (home team's) skills
  matrix[n_games,K] X_away;    // design matrix for (away team's) skills
  int <lower=0,upper=3> home_sets[n_games]; // (home team's) winning sets
  int <lower=0,upper=3> away_sets[n_games]; // away team winning sets
  int gammas_home[K]; // binary indicators of (home team's) candidate variables
  int gammas_away[K]; // binary indicators of (away team's) candidate variables
  vector[K] post_mean_beta_home; //posterior means of (home team's) betas from a MCMC pilot run of full model
  vector[K] post_mean_beta_away; //posterior means of (away team's) betas from a MCMC pilot run of full model
  vector[K] post_sd_beta_home; //posterior stan. deviat. of (home team's) betas from a MCMC pilot run of full model
  vector[K] post_sd_beta_away; //posterior stan. deviat. of (away team's) betas from a MCMC pilot run of full model
  real<lower=0> c_thres;//c: upper threshold multiplicator for lambdas parameters
  real<lower=0> c_std;//c: upper threshold multiplicator for lambdas parameters
}


parameters {
  
  vector[K] beta_home; // parameters of (home team's) candidate variables 
  vector[K] beta_away; // parameters of (away team's) candidate variables  
  
  real mu; // constant parameter
  real home; // common home effect

}

transformed parameters {


  vector[n_games]   lambda1;
  vector[n_games]   lambda2;
  vector[n_games]   lambda1_star;
  vector[n_games]   lambda2_star; 
  vector[K] gb_home;
  vector[K] gb_away;
  
  for (j in 1:K){
    gb_home[j]=gammas_home[j]*beta_home[j];
    gb_away[j]=gammas_away[j]*beta_away[j];
  }

  // sum to zero constraint for both attacking and defensive parameters


  // Linear predictor
  lambda1= exp(mu+home+X_home * gb_home);          
  lambda2= exp(mu+X_away* gb_away); 
  
  // We specified an upper bound in order to avoid numeric overflow problem for model parameters
   for (g in 1:n_games) {
     if (lambda1[g]>(100*c_thres)){
      lambda1_star[g]=(100*c_thres);
    } else {
       lambda1_star[g]=lambda1[g];
    }
    if (lambda2[g]>(100*c_thres)){
       lambda2_star[g]=(100*c_thres);
    } else {
       lambda2_star[g]=lambda2[g];
     }
   }
}

model {

  // Priors of all parameters except for the candidate variables' parameters
  
  target+=normal_lpdf(mu|0,0.37);
  target+=normal_lpdf(home|0,0.37);


  // Prior+Pseudoprior of the candidate variables' parameters
    for (j in 1:K){
     target+=normal_lpdf(beta_home[j]|gammas_home[j] *0+(1-gammas_home[j])*post_mean_beta_home[j],(sqrt(n_games)^gammas_home[j])*post_sd_beta_home[j]);
     target+=normal_lpdf(beta_away[j]|gammas_away[j] *0+(1-gammas_away[j])*post_mean_beta_away[j],(sqrt(n_games)^gammas_away[j])*post_sd_beta_away[j]);
  }
  
 // Likelihood
  
  for (g in 1:n_games){
    target+=skellam_without_lpmf((home_sets[g]-away_sets[g])|lambda1_star[g],lambda2_star[g]);
    
  }
}

generated quantities{
  matrix[n_games,2*K] log_lik_one; // matrix with log likelihoods when gamma[j]=1
  matrix[n_games,2*K] log_lik_zero; // matrix with log likelihoods when gamma[j]=0

  // the gammas*betas when we want to estimate the log likelihood
  // for the opposite value of binary indicator within the if-loop.
  vector[K] gb_new_home;
  vector[K] gb_new_away;

  // the linear predictor when we want to estimate the 
  // log likelihood for the opposite value of binary indicator within the if-loop.

  vector[n_games] lambda1_new;
  vector[n_games] lambda2_new;
    vector[n_games] lambda1_star_new;
  vector[n_games] lambda2_star_new;
  
   for (i in 1:n_games) {
     for (j in 1:K){

      // In each case of if-loop, we estimate the log likelihood for both cases
      // of binary indicators gamma
      if (gb_home[j]==0){
         log_lik_zero[i,j] = skellam_without_lpmf(home_sets[i]-away_sets[i] |lambda1_star[i],lambda2_star[i]);
         gb_new_home=gb_home;
         gb_new_home[j]=1*beta_home[j];
         lambda1_new[i]=exp(mu+home+X_home[i,]*gb_new_home);

	// We specified an upper bound in order to avoid numeric overflow problem for model parameters
         if (lambda1_new[i]>c_thres*100){
              lambda1_star_new[i]=c_thres*100;
          } else {
               lambda1_star_new[i]=lambda1_new[i];
          }
         
         log_lik_one[i,j]=  skellam_without_lpmf(home_sets[i]-away_sets[i] |lambda1_star_new[i],lambda2_star[i]);
      } else {
         log_lik_one[i,j]=skellam_without_lpmf(home_sets[i]-away_sets[i] |lambda1_star[i],lambda2_star[i]);
         gb_new_home=gb_home;
         gb_new_home[j]=0;
         lambda1_new[i]=exp(mu+home+X_home[i,]*gb_new_home);
         
	// We specified an upper bound in order to avoid numeric overflow problem for model parameters
          if (lambda1_new[i]>c_thres*100){
              lambda1_star_new[i]=c_thres*100;
          } else {
              lambda1_star_new[i]=lambda1_new[i];
          }
          
         log_lik_zero[i,j]=skellam_without_lpmf(home_sets[i]-away_sets[i] |lambda1_star_new[i],lambda2_star[i]);
      }
      
      if (gb_away[j]==0){
         log_lik_zero[i,j+K] = skellam_without_lpmf(home_sets[i]-away_sets[i] |lambda1_star[i],lambda2_star[i]);
         gb_new_away=gb_away;
         gb_new_away[j]=1*beta_away[j];
         lambda2_new[i]=exp(mu+X_away[i,]*gb_new_away);
        
	// We specified an upper bound in order to avoid numeric overflow problem for model parameters
        if (lambda2_new[i]>c_thres*100){
             lambda2_star_new[i]=c_thres*100;
        } else {
             lambda2_star_new[i]=lambda2_new[i];
        }
          
        log_lik_one[i,j+K]=  skellam_without_lpmf(home_sets[i]-away_sets[i] |lambda1_star[i],lambda2_star_new[i]);
        
      } else {
         log_lik_one[i,j+K]=skellam_without_lpmf(home_sets[i]-away_sets[i] |lambda1_star[i],lambda2_star[i]);
         gb_new_away=gb_away;
         gb_new_away[j]=0;
         lambda2_new[i]=exp(mu+X_away[i,]*gb_new_away);
              
	// We specified an upper bound in order to avoid numeric overflow problem for model parameters
         if (lambda2_new[i]>c_thres*100){
             lambda2_star_new [i]=c_thres*100;
         } else {
            lambda2_star_new[i]=lambda2_new[i];
         }
        
        log_lik_zero[i,j+K]=skellam_without_lpmf(home_sets[i]-away_sets[i] |lambda1_star[i],lambda2_star_new[i]);

      }
    }
  }
}
