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
  int<lower=1> K;       // Dimension of model matrix (without intercept because mu is our intercept)
  matrix[n_games, K] X_home;   //Model matrix for home teams
  matrix[n_games,K] X_away;    // Model matrix for away teams
  int <lower=0,upper=3> home_sets[n_games];//0-3 sets can have each team
  int <lower=0,upper=3> away_sets[n_games];//0-3 sets can have each team
  int home_team[n_games];
  int away_team[n_games];
}


parameters {
  
  vector[K] beta_home;
  vector[K] beta_away;
  real mu;
  real home;
  real attack_raw[n_teams - 1];
  real defense_raw[n_teams - 1];
}

transformed parameters {
  // Enforce sum-to-zero constraints
  vector[n_teams]   attack;
  vector[n_teams]   defense;
  vector[n_games]   lambda1_star;
  vector[n_games]   lambda2_star; 
  vector[n_games]   lambda1;
  vector[n_games]   lambda2;
  
  for (t in 1:(n_teams-1)) {
    attack[t] = attack_raw[t];
    defense[t] = defense_raw[t];
  }
  
  attack[n_teams] = -sum(attack_raw);
  defense[n_teams] = -sum(defense_raw);
  // Creation of linear predictor
  lambda1_star= exp(mu+home+attack[home_team]+defense[away_team]+X_home * beta_home);          
  lambda2_star= exp(mu+attack[away_team]+defense[home_team]+X_away* beta_away);  
 
  for (g in 1:n_games){
    if (lambda1_star[g]>100.0){
      lambda1[g]=100.0;
    } else {
      lambda1[g]=lambda1_star[g];
    }
    if (lambda2_star[g]>100.0){
      lambda2[g]=100.0;
    } else {
      lambda2[g]=lambda2_star[g];
    }
    
  }
}

model {
  
  
  //Priors
  target+=normal_lpdf(beta_home|0,1);
  target+=normal_lpdf(beta_away|0,1);
  target+=normal_lpdf(mu|0,0.37);
  target+=normal_lpdf(home|0,0.37);
  target+=normal_lpdf(attack|0,1);
  target+=normal_lpdf(defense|0,1);
  
  
  
  //likelihood-systematic component
  
  for (g in 1:n_games){
    target+=skellam_without_lpmf((home_sets[g]-away_sets[g])|lambda1[g],lambda2[g]);
    
  }
}

