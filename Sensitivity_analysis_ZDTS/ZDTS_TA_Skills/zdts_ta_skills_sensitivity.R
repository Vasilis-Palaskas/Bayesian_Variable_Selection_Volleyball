# Load the proper libraries.
library(rstan)
library(coda)
library(bayesplot)
library(skellam)
options(mc.cores = parallel::detectCores())# Activate multiple cores for stan models

#---Data Preparation
source(file.choose())#-Data_Preparation.R

# Choose the working directory of this file (...\\Submitted_Appendix\\ZDTS\\)
# setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/ZDTS_TA_Skills")

#------Skills for both Home and Away Teams
X_home<-data_by_sets[c(
  "Home_perfect_serves","Home_very_good_serves",
  "Home_failed_serves","Home_perfect_passes","Home_very_good_passes",
  "Home_poor_passes","Home_failed_passes","Home_perfect_att1",
  "Home_blocked_att1","Home_failed_att1","Home_perfect_att2",
  "Home_blocked_att2","Home_failed_att2","Home_perfect_blocks",
  "Home_net_violation_blocks","Home_failed_blocks","Home_failed_settings")
]

X_away<-data_by_sets[c(
  "Away_perfect_serves","Away_very_good_serves",
  "Away_failed_serves","Away_perfect_passes","Away_very_good_passes",
  "Away_poor_passes","Away_failed_passes","Away_perfect_att1",
  "Away_blocked_att1","Away_failed_att1","Away_perfect_att2",
  "Away_blocked_att2","Away_failed_att2","Away_perfect_blocks",
  "Away_net_violation_blocks","Away_failed_blocks","Away_failed_settings")
]



#----Rename properly the skill variables
##----Skill events selected via the BVS process based on PSI Median Threshold
#### Standardization of the Model Matrices for numerical convenience

X_home_std<-data.frame(scale(X_home,center=T,scale=T) )
X_away_std<-data.frame(scale(X_away,center=T,scale=T) )
## Vector of teams names along with
## their ranking positions, points, abilities
teams <- levels(data_by_sets$home_Team)
observed_positions<-c("(7)","(6)","(9)","(8)","(5)","(11)","(1)","(12)","(4)","(10)","(3)","(2)")
observed_points<-c("(36)","(37)","(16)","(28)","(38)","(14)","(62)","(7)","(39)","(16)","(50)","(53)")


teams_attack<-paste0(teams," ","Attack")
teams_defense<-paste0(teams," ","Defense")
teams_over<-paste0(teams," ","Overall")

teams_pos<-paste0(teams," ",observed_positions)
teams_points<-paste0(teams," ",observed_points)
#-------------- Sensitivity Analysis
#--------
###-----Datalists required for the Bayesian model fitting across several values of c

c_thres<-c(1/20,1/10,1/2,1,2,5,8,9,10)#- c_thres: l1,l2 threshold multiplicator for betas parameters
c_std<-c(1/2,1,2,5,8,9,10)#- c_std: prior standard deviation multiplicator for betas parameters
zdts_support<-c(-3,-2,-1,1,2,3)#- Support of the random variable of set-difference

##---Initializing Matrix for minimum deviances-divergent transitions

min_dev_vector<-mean_dev_vector<-median_dev_vector<-sd_dev_vector<-
  div_trans_vector<-deviance_median_vector<-NULL

 for (i in c_thres){
     for (j in c_std){
      data_std_thres_zdts_skills<-list(c_thres=i,c_std=j,
                                       n_games=dim(data_by_sets)[1],
                                       away_team=as.numeric(data_by_sets$away_Team),
                                       home_team=as.numeric(data_by_sets$home_Team),
                                       n_teams=length(levels(data_by_sets$home_Team)),
                                       X_home=X_home_std,X_away=X_away_std,
                                       K=ncol(X_home_std),
                                       home_sets=data_by_sets$home_sets,
                                       away_sets=data_by_sets$away_sets)
      
      # Define the model implemented for sensitivity analysis
      sensit_betas_zdts_skills_std_thres.stan=
        "functions {
  
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

  vector[n_games]   lambda1;
  vector[n_games]   lambda2;
  vector[n_games]   lambda1_star;
  vector[n_games]   lambda2_star; 
  
  for (t in 1:(n_teams-1)) {
    attack[t] = attack_raw[t];
    defense[t] = defense_raw[t];
  }
  
  attack[n_teams] = -sum(attack_raw);
  defense[n_teams] = -sum(defense_raw);
  // Creation of linear predictor
  lambda1= exp(mu+home+attack[home_team]+defense[away_team]+X_home * beta_home);          
  lambda2= exp(mu+attack[away_team]+defense[home_team]+X_away* beta_away);    
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
  

  //Priors
  target+=normal_lpdf(beta_home|0,1*c_std);
  target+=normal_lpdf(beta_away|0,1*c_std);
  target+=normal_lpdf(mu|0,0.37);
  target+=normal_lpdf(home|0,0.37);
  target+=normal_lpdf(attack|0,1);
  target+=normal_lpdf(defense|0,1);
  
  
  
  //likelihood-systematic component
  for (g in 1:n_games) {
    target+=skellam_without_lpmf(home_sets[g]-away_sets[g]|lambda1_star[g],lambda2_star[g]) ;
  }
  
}
generated quantities{
  vector[n_games] log_lik;
  real dev;

  dev=0;
    for (g in 1:n_games) {
        log_lik[g] =skellam_without_lpmf(home_sets[g]-away_sets[g]|lambda1_star[g],lambda2_star[g]) ;
        dev=dev-2*log_lik[g];

    }

  overall=attack-defense;

}

"

#---Bayesian ZDTS Models Running
warmup_iters<-1000
total_iters<-6000
full_zdts_ta_skills<-stan(model_code=sensit_betas_zdts_skills_std_thres.stan,
                       data= data_std_thres_zdts_skills,thin=1,
                       chains=1,cores=1,
                       iter=total_iters,warmup=warmup_iters,
                      init_r=1,seed="1234")

# Extraction of the candidate models' deviances (Table 2)
dev_full_zdts_ta_skills<-extract(full_zdts_ta_skills,pars="dev")
#----Deviance estimated based on the median of l1, l2 parameters (more stable estimation of the deviance)
l1_star<-extract(full_zdts_ta_skills,pars="lambda1_star")
l2_star<-extract(full_zdts_ta_skills,pars="lambda2_star")
l1_star_median<-apply(l1_star$lambda1_star,2,median)
l2_star_median<-apply(l2_star$lambda2_star,2,median)
###--------------Deviance calculation based on median parameters-----------###
#-Initialisation
deviance_median<-0#-Deviance estimated based on the median of l1, l2
log_lik_median<-NULL#-Log likelihood required estimated based on the median of l1, l2

for (l in 1:data_std_thres_zdts_skills$n_games){
  log_lik_median<-log(dskellam(data_std_thres_zdts_skills$home_sets[l]-
                                 data_std_thres_zdts_skills$away_sets[l],
                               l1_star_median[l],l2_star_median[l])/sum(
                                 dskellam(zdts_support,l1_star_median[l],
                                          l2_star_median[l])
                               )
  )
  deviance_median=deviance_median-2*log_lik_median
}
##--Summary Statistics of Model Deviances-divergent transitions
min_dev_vector<-c(min_dev_vector,min(dev_full_zdts_ta_skills$dev))#- Posterior Minimum Deviance
mean_dev_vector<-c(mean_dev_vector,mean(dev_full_zdts_ta_skills$dev))#-Posterior Mean Deviance
sd_dev_vector<-c(sd_dev_vector,sd(dev_full_zdts_ta_skills$dev))#-Posterior Std of Deviance
median_dev_vector<-c(median_dev_vector,median(dev_full_zdts_ta_skills$dev))#-Posterior Median of Deviance
divergent <- get_sampler_params(full_zdts_ta_skills, inc_warmup=FALSE)[[1]][,'divergent__']#-Posterior Divergent transitions
div_trans_vector<-c(div_trans_vector,sum(divergent)/(2*(total_iters-warmup_iters)))#-%Posterior Divergent transitions
deviance_median_vector<-c(deviance_median_vector,deviance_median)#-Deviance estimated based on the median of l1, l2
     }
}
#---Tables with summary statistics of deviances and divergent transitions
table_min_dev<-matrix(min_dev_vector,ncol=7,nrow=9)
table_mean_dev<-matrix(mean_dev_vector,ncol=7,nrow=9)
table_median_dev<-matrix(median_dev_vector,ncol=7,nrow=9)
table_sd_dev<-matrix(sd_dev_vector,ncol=7,nrow=9)
table_div_trans<-matrix(div_trans_vector,ncol=7,nrow=9)
table_deviance_median<-matrix(deviance_median_vector,ncol=7,nrow=9)
# Rounding
table_min_dev<-round(table_min_dev,1)
table_div_trans<-round(table_div_trans,2)
table_median_dev<-round(table_median_dev,1)
table_mean_dev<-round(table_mean_dev,1)
table_sd_dev<-round(table_sd_dev,2)
table_deviance_median<-round(table_deviance_median,1)
#----Save results
write.csv(table_min_dev,file="table_deviances_zdts_ta_skills.csv",row.names=F)
write.csv(table_div_trans,file="table_div_ta_trans.csv",row.names=F)
write.csv(table_median_dev,file="table_median_deviances_zdts_ta_skills.csv",row.names=F)
write.csv(table_mean_dev,file="table_mean_deviances_zdts_ta_skills.csv",row.names=F)
write.csv(table_sd_dev,file="table_sd_deviances_zdts_ta_skills.csv",row.names=F)
write.csv(table_deviance_median,file="table_deviance_median_zdts_ta_skills.csv",row.names=F)




