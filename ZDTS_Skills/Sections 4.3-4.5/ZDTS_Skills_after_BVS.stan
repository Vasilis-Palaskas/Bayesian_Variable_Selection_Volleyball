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
  int <lower=1> n_games; //number of games 132
  int <lower=1> n_teams; //number of teams 12
  int<lower=1> K_home;       // number of home skill variables
    int<lower=1> K_away;       // number of away skill variables

  matrix[n_games, K_home] X_home;   // design matrix for (home team's) skills
  matrix[n_games,K_away] X_away;    // design matrix for (away team's) skills
  int <lower=0,upper=3> home_sets[n_games];//0-3 sets can have each team
  int <lower=0,upper=3> away_sets[n_games];//0-3 sets can have each team
}

parameters {
  
  vector[K_home] beta_home;
  vector[K_away] beta_away;
  real mu;
  real home;
}

transformed parameters {
  // Enforce sum-to-zero constraints

  vector[n_games]   lambda1_star;
  vector[n_games]   lambda2_star; 
  vector[n_games]   lambda1;
  vector[n_games]   lambda2;
  
  // Creation of linear predictor
  lambda1_star= exp(mu+X_home * beta_home+home);          
  lambda2_star= exp(mu+X_away * beta_away);  
      lambda1=lambda1_star;
    lambda2=lambda2_star;
  // for (g in 1:n_games) {
  //   if (lambda1_star[g]>150.0){
  //     lambda1[g]=150.0;
  //   } else {
  //     lambda1[g]=lambda1_star[g];
  //   }
  //   if (lambda2_star[g]>150.0){
  //     lambda2[g]=150.0;
  //   } else {
  //     lambda2[g]=lambda2_star[g];
  //   }
  // }
}

model {

  
  //Priors
  
  //Priors
  target+=normal_lpdf(beta_home|0,1);
  target+=normal_lpdf(beta_away|0,1);
  target+=normal_lpdf(mu|0,0.37);
  target+=normal_lpdf(home|0,0.37);


  
  //likelihood-systematic component
  for (i in 1:n_games) {
    target+=skellam_without_lpmf(home_sets[i]-away_sets[i]|lambda1[i],lambda2[i]) ;
  }
  
}


generated quantities{
  vector[n_games] log_lik;
  real dev;
  // vector[n_teams]   overall;// overall ability

  //real DIC;

  dev=0;
  for (i in 1:n_games) {
    log_lik[i] = skellam_without_lpmf(home_sets[i]-away_sets[i] |lambda1[i],lambda2[i]);
    dev=dev-2*log_lik[i];
  }
  // overall=attack-defense;
  //DIC=mean(dev)+0.5*variance(dev);
}
