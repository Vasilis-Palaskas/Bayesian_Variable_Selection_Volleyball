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

      log_prob_new=skellam_lpmf(k|mu1,mu2)-normalization;
      return log_prob_new;
    
  }
}



data {
  int <lower=1> n_games; //number of games 132
  int <lower=1> n_teams; //number of teams 12
  int<lower=1> K;       // number of candidate variables
  matrix[n_games, K] X_home;   // design matrix for (home team's) skills
  matrix[n_games,K] X_away;    // design matrix for (away team's) skills
  int <lower=0,upper=3> home_sets[n_games];//0-3 sets can have each team
  int <lower=0,upper=3> away_sets[n_games];//0-3 sets can have each team
  real<lower=0> c_thres;//c: upper threshold multiplicator for lambdas parameters
  real<lower=0> c_std;//c: upper threshold multiplicator for lambdas parameters
 
}

parameters {
  
  vector[K] beta_home;
  vector[K] beta_away;
  real mu;
  real home;
 
}

transformed parameters {


  vector[n_games]   lambda1;
  vector[n_games]   lambda2;
  vector[n_games]   lambda1_star;
  vector[n_games]   lambda2_star; 
 

  // Creation of linear predictor
  lambda1= exp(mu+home+X_home * beta_home);          
  lambda2= exp(mu+X_away* beta_away);    
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
  
  int    sets_diff[n_games];
  for (g in 1:n_games){
    sets_diff[g]=home_sets[g]-away_sets[g];
    
  }
  

  //Priors
  target+=normal_lpdf(beta_home|0,1*c_std);
  target+=normal_lpdf(beta_away|0,1*c_std);
  target+=normal_lpdf(mu|0,0.37);
  target+=normal_lpdf(home|0,0.37);

  
  
  //likelihood-systematic component
  for (g in 1:n_games) {
    target+=skellam_without_lpmf(sets_diff[g]|lambda1_star[g],lambda2_star[g]) ;
  }
  
}
generated quantities{
  vector[n_games] log_lik;
 // vector[n_games] log_lik_star;
  real dev;

  //dev=0;
  dev=0;
    for (g in 1:n_games) {
        log_lik[g] =skellam_without_lpmf(home_sets[g]-away_sets[g]|lambda1_star[g],lambda2_star[g]) ;
        dev=dev-2*log_lik[g];

    }

  //overall=attack-defense;
  //DIC=mean(dev)+0.5*variance(dev);
}
