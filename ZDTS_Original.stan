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
  real nomin;
  real denom;
  
  for (t in 1:(n_teams-1)) {
    attack[t] = attack_raw[t];
    defense[t] = defense_raw[t];
  }
  
  attack[n_teams] = -sum(attack_raw);
  defense[n_teams] = -sum(defense_raw);
  // Creation of linear predictor
  lambda1_star= exp(mu+home+attack[home_team]+defense[away_team]+X_home * beta_home);          
  lambda2_star= exp(mu+attack[away_team]+defense[home_team]+X_away* beta_away);    
  nomin= (
    modified_bessel_first_kind(abs(1), 2 * sqrt(lambda1_star*lambda2_star))+
    +(2*modified_bessel_first_kind(abs(2), 2 * sqrt(lambda1_star*lambda2_star))*
    (lambda1_star+lambda2_star)/
    sqrt(lambda1_star*lambda2_star)
    )+(3*modified_bessel_first_kind(abs(3), 2 * sqrt(lambda1_star*lambda2_star))*
    (lambda1_star^2+lambda1_star*lambda2_star+lambda2_star^2)/
        (lambda1_star*lambda2_star)
    ));
    
    denom=modified_bessel_first_kind(abs(1), 2 * sqrt(lambda1_star*lambda2_star))*
    (lambda1_star+lambda2_star)+
   ( (modified_bessel_first_kind(abs(2), 2 * sqrt(lambda1_star*lambda2_star))*
    (lambda1_star^2+lambda2_star^2))/
     sqrt(lambda1_star*lambda2_star)   )+ 
     (modified_bessel_first_kind(abs(3), 2 * sqrt(lambda1_star*lambda2_star))*
    (lambda1_star^3+lambda2_star^3) )/
     (lambda1_star*lambda2_star)   );
  lambda1=lambda1_star*(nomin/denom)
     lambda2=lambda2_star*(nomin/denom)

  //  for (g in 1:n_games) {
  //    if (lambda1_star[g]>(100*c_thres)){
  //     lambda1[g]=(100*c_thres);
  //   } else {
  //      lambda1[g]=lambda1_star[g];
  //   }
  //   if (lambda2_star[g]>(100*c_thres)){
  //     lambda2[g]=(100*c_thres);
  //  } else {
  //      lambda2[g]=lambda2_star[g];
  //   }
  // }
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
  target+=normal_lpdf(attack|0,1);
  target+=normal_lpdf(defense|0,1);
  
  
  
  //likelihood-systematic component
  for (g in 1:n_games) {
    target+=skellam_without_lpmf(sets_diff[g]|lambda1[g],lambda2[g]) ;
  }
  
}
generated quantities{
  vector[n_games] log_lik;
 // vector[n_games] log_lik_star;
  real dev;

  //dev=0;
  dev=0;
    for (g in 1:n_games) {
        log_lik[g] =skellam_without_lpmf(home_sets[g]-away_sets[g]|lambda1[g],lambda2[g]) ;
        dev=dev-2*log_lik[g];
        //log_lik_star[g] =skellam_without_lpmf(home_sets[g]-away_sets[g]|lambda1_star[g],lambda2_star[g]) ;
        //dev=dev-2*log_lik_star[g];
    }

  //overall=attack-defense;
  //DIC=mean(dev)+0.5*variance(dev);
}
