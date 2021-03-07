# Load the proper libraries.
library(rstan)
library(coda)
library(bayesplot)

# Activate multiple cores for stan models
options(mc.cores = parallel::detectCores())# Choose the working directory of this file (...\\Submitted_Appendix\\ZDTS\\)
setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/ZDTS_TA_Skills")
# Choose the working directory of this file (...\\Submitted_Appendix\\ZDTS\\)
# Load the properly prepared data for both home and away skill events as well as
# both home and away teams in each match
load("X_home")
load("X_away")
load("data_zdts_skills")




#### Standardization of the Model Matrices for numerical convenience
X_home_std<-X_away_std<-matrix(NA,nrow=132,ncol=17)

for (i in 1:dim(X_home)[2]){
  X_home_std[,i]<-scale(X_home[,i])
  X_away_std[,i]<-scale(X_away[,i])
}

colnames(X_home_std)<-c("(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass",
                        "(Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
                        "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
                        "(Home) block net violation","(Home) failed block","(Home) failed setting")

colnames(X_away_std)<-c("(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass",
                        "(Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
                        "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
                        "(Away) block net violation","(Away) failed block","(Away) failed setting")


#------------------------------------------
########-------------- Sensitivity Analysis


###-----Datalists required for the Bayesian model fitting across several values of c
###----- c: prior standard deviation multiplicator for betas parameters

data_zdts_skills_c_1<-list(c=1,n_games=data_zdts_skills$N,
                           away_team=as.numeric(data_zdts_skills$away_team),
                           home_team=as.numeric(data_zdts_skills$home_team),
                           n_teams=data_zdts_skills$n_teams,
                           X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
                           home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)

data_zdts_skills_c_2<-list(c=2,n_games=data_zdts_skills$N,
                           away_team=as.numeric(data_zdts_skills$away_team),
                           home_team=as.numeric(data_zdts_skills$home_team),
                           n_teams=data_zdts_skills$n_teams,
                           X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
                           home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)


data_zdts_skills_c_5<-list(c=5,n_games=data_zdts_skills$N,
                           away_team=as.numeric(data_zdts_skills$away_team),
                           home_team=as.numeric(data_zdts_skills$home_team),
                           n_teams=data_zdts_skills$n_teams,
                           X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
                           home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)

data_zdts_skills_c_10<-list(c=10,n_games=data_zdts_skills$N,
                            away_team=as.numeric(data_zdts_skills$away_team),
                            home_team=as.numeric(data_zdts_skills$home_team),
                            n_teams=data_zdts_skills$n_teams,
                            X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
                            home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)





# Define the model implemented for sensitivity analysis
sensit_betas_zdts_skills_different_thresholds.stan=
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
  int<lower=1> K;       // number of candidate variables
  matrix[n_games, K] X_home;   // design matrix for (home team's) skills
  matrix[n_games,K] X_away;    // design matrix for (away team's) skills
  int <lower=0,upper=3> home_sets[n_games];//0-3 sets can have each team
  int <lower=0,upper=3> away_sets[n_games];//0-3 sets can have each team
  real<lower=0> c;//c: upper threshold multiplicator for lambdas parameters
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
   for (g in 1:n_games) {
     if (lambda1_star[g]>(100*c)){
      lambda1[g]=(100*c);
    } else {
       lambda1[g]=lambda1_star[g];
    }
    if (lambda2_star[g]>(100*c)){
       lambda2[g]=(100*c);
    } else {
       lambda2[g]=lambda2_star[g];
     }
   }
}

model {
  
  int    sets_diff[n_games];
  for (g in 1:n_games){
    sets_diff[g]=home_sets[g]-away_sets[g];
    
  }
  

  //Priors
  target+=normal_lpdf(beta_home|0,10);
  target+=normal_lpdf(beta_away|0,10);
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
 // real dev;

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

"


# Extraction of the candidate models' deviances (Table 2)

# c=1
full_zdts_skills_c_1_different_thresholds<-stan(model_code=sensit_betas_zdts_skills_different_thresholds.stan,data=data_zdts_skills_c_1,thin=1,chains=2,cores=2,
                                                                          iter=10000,warmup=2000,seed="12345",init_r=1)#15948/16000=99.6% divergent transitions after warmup

save(full_zdts_skills_c_1_different_thresholds,file="full_zdts_skills_c_1_different_thresholds")


# c=2
full_zdts_skills_c_2_different_thresholds<-stan(model_code=sensit_betas_zdts_skills_different_thresholds.stan,data=data_zdts_skills_c_2,thin=1,chains=2,cores=2,
                                                                          iter=10000,warmup=2000,seed="12345",init_r=1)#15948/16000=99.6% divergent transitions after warmup

save(full_zdts_skills_c_2_different_thresholds,file="full_zdts_skills_c_2_different_thresholds")

# c=5
full_zdts_skills_c_5_different_thresholds<-stan(model_code=sensit_betas_zdts_skills_different_thresholds.stan,data=data_zdts_skills_c_5,thin=1,chains=2,cores=2,
                                                                          iter=10000,warmup=2000,seed="12345",init_r=1)#15980/16000=99.8%divergent transitions after warmup
save(full_zdts_skills_c_5_different_thresholds,file="full_zdts_skills_c_5_different_thresholds")

# c=10
full_zdts_skills_c_10_different_thresholds<-stan(model_code=sensit_betas_zdts_skills_different_thresholds.stan,data=data_zdts_skills_c_10,thin=1,chains=2,cores=2,
                                                      iter=10000,warmup=2000,seed="12345",init_r=1)#Effective Samples Size (ESS) is too low, 

save(full_zdts_skills_c_10_different_thresholds,file="full_zdts_skills_c_10_different_thresholds")

##---Save several fitted models (c=1,2,5,10)

# launch_shinystan(full_zdts_skills_c_1_different_thresholds)# Not converged
# launch_shinystan(full_zdts_skills_c_2_different_thresholds)# Not converged
# launch_shinystan(full_zdts_skills_c_5_different_thresholds)# Not converged
# launch_shinystan(full_zdts_skills_c_10_different_thresholds)# Not converged


####----- Posterior distributions of betas parameters across several c values
####-----  c: prior standard deviation multiplicator for betas parameters
beta_home_full_zdts_skills_c_1_different_thresholds<-extract(full_zdts_skills_c_1_different_thresholds,pars="beta_home")
beta_away_full_zdts_skills_c_1_different_thresholds<-extract(full_zdts_skills_c_1_different_thresholds,pars="beta_away")

beta_home_full_zdts_skills_c_2_different_thresholds<-extract(full_zdts_skills_c_2_different_thresholds,pars="beta_home")
beta_away_full_zdts_skills_c_2_different_thresholds<-extract(full_zdts_skills_c_2_different_thresholds,pars="beta_away")

beta_home_full_zdts_skills_c_5_different_thresholds<-extract(full_zdts_skills_c_5_different_thresholds,pars="beta_home")
beta_away_full_zdts_skills_c_5_different_thresholds<-extract(full_zdts_skills_c_5_different_thresholds,pars="beta_away")

beta_home_full_zdts_skills_c_10_different_thresholds<-extract(full_zdts_skills_c_10_different_thresholds,pars="beta_home")
beta_away_full_zdts_skills_c_10_different_thresholds<-extract(full_zdts_skills_c_10_different_thresholds,pars="beta_away")

beta_home_full_zdts_skills_c_1_different_thresholds<-as.data.frame(beta_home_full_zdts_skills_c_1_different_thresholds)
beta_away_full_zdts_skills_c_1_different_thresholds<-as.data.frame(beta_away_full_zdts_skills_c_1_different_thresholds)
beta_home_full_zdts_skills_c_2_different_thresholds<-as.data.frame(beta_home_full_zdts_skills_c_2_different_thresholds)
beta_away_full_zdts_skills_c_2_different_thresholds<-as.data.frame(beta_away_full_zdts_skills_c_2_different_thresholds)
beta_home_full_zdts_skills_c_5_different_thresholds<-as.data.frame(beta_home_full_zdts_skills_c_5_different_thresholds)
beta_away_full_zdts_skills_c_5_different_thresholds<-as.data.frame(beta_away_full_zdts_skills_c_5_different_thresholds)
beta_home_full_zdts_skills_c_10_different_thresholds<-as.data.frame(beta_home_full_zdts_skills_c_10_different_thresholds)
beta_away_full_zdts_skills_c_10_different_thresholds<-as.data.frame(beta_away_full_zdts_skills_c_10_different_thresholds)

###-----------Rename the coefficients of stan fit objects
colnames(beta_home_full_zdts_skills_c_1_different_thresholds)<-colnames(beta_home_full_zdts_skills_c_2_different_thresholds)<-colnames(beta_home_full_zdts_skills_c_5_different_thresholds)<-
  colnames(beta_home_full_zdts_skills_c_10_different_thresholds)<-c(
    "(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass",
    "(Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
    "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
    "(Home) block net violation","(Home) failed block","(Home) failed setting")

colnames(beta_away_full_zdts_skills_c_1_different_thresholds)<-colnames(beta_away_full_zdts_skills_c_2_different_thresholds)<-colnames(beta_away_full_zdts_skills_c_5_different_thresholds)<-
  colnames(beta_away_full_zdts_skills_c_10_different_thresholds)<-c(
    "(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass",
    "(Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
    "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
    "(Away) block net violation","(Away) failed block","(Away) failed setting")




##---------------------------------------------------------
##-------1) Deviances Comparison across several values of c

#----c=1
dev_full_zdts_skills_c_1_different_thresholds<-extract(full_zdts_skills_c_1_different_thresholds,pars="dev")
dev_full_zdts_skills_c_1_different_thresholds$dev

#----c=2
dev_full_zdts_skills_c_2_different_thresholds<-extract(full_zdts_skills_c_2_different_thresholds,pars="dev")
dev_full_zdts_skills_c_2_different_thresholds$dev

#----c=5

dev_full_zdts_skills_c_5_different_thresholds<-extract(full_zdts_skills_c_5_different_thresholds,pars="dev")
dev_full_zdts_skills_c_5_different_thresholds$dev


#----c=10

dev_full_zdts_skills_c_10_different_thresholds<-extract(full_zdts_skills_c_10_different_thresholds,pars="dev")
dev_full_zdts_skills_c_10_different_thresholds$dev


#---- 1)Sensitivity Analysis: Deviances Comparison Table (c=2,5,10)
mean(dev_full_zdts_skills_c_1_different_thresholds$dev)#215.3
sd(dev_full_zdts_skills_c_1_different_thresholds$dev)#11.2
mean(dev_full_zdts_skills_c_2_different_thresholds$dev)#212.4
sd(dev_full_zdts_skills_c_2_different_thresholds$dev)#9.2
mean(dev_full_zdts_skills_c_5_different_thresholds$dev)#210.0
sd(dev_full_zdts_skills_c_5_different_thresholds$dev)# 9.1
mean(dev_full_zdts_skills_c_10_different_thresholds$dev)#208.9
sd(dev_full_zdts_skills_c_10_different_thresholds$dev)# 9.9


min(dev_full_zdts_skills_c_1_different_thresholds$dev)#
min(dev_full_zdts_skills_c_2_different_thresholds$dev)
min(dev_full_zdts_skills_c_5_different_thresholds$dev)
min(dev_full_zdts_skills_c_10_different_thresholds$dev)
#---- 2)Sensitivity Analysis:  
rstan::check_divergences(full_zdts_skills_c_1_different_thresholds)#0%
rstan::check_divergences(full_zdts_skills_c_2_different_thresholds)#0%
rstan::check_divergences(full_zdts_skills_c_5_different_thresholds)#50.2%
rstan::check_divergences(full_zdts_skills_c_10_different_thresholds)#86.2%



#---- 3) Posterior summary statistics of betas parameters
###---
names(full_zdts_skills_c_1_different_thresholds)[c(1:17)]<-names(full_zdts_skills_c_2_different_thresholds)[c(1:17)]<-
  names(full_zdts_skills_c_5_different_thresholds)[c(1:17)]<-names(full_zdts_skills_c_10_different_thresholds)[c(1:17)]<-c(
    "(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass",
    "(Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
    "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
    "(Home) block net violation","(Home) failed block","(Home) failed setting")

names(full_zdts_skills_c_1_different_thresholds)[c(18:34)]<-names(full_zdts_skills_c_2_different_thresholds)[c(18:34)]<-
  names(full_zdts_skills_c_5_different_thresholds)[c(18:34)]<-names(full_zdts_skills_c_10_different_thresholds)[c(18:34)]<-c(
    "(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass",
    "(Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
    "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
    "(Away) block net violation","(Away) failed block","(Away) failed setting")


full_zdts_skills_c_1_different_thresholds_summary<-summary(full_zdts_skills_c_1_different_thresholds,
                                                                pars=c("beta_home","beta_away"))
full_zdts_skills_c_2_different_thresholds_summary<-summary(full_zdts_skills_c_2_different_thresholds,
                                                                pars=c("beta_home","beta_away"))
full_zdts_skills_c_5_different_thresholds_summary<-summary(full_zdts_skills_c_5_different_thresholds,
                                                                pars=c("beta_home","beta_away"))
full_zdts_skills_c_10_different_thresholds_summary<-summary(full_zdts_skills_c_10_different_thresholds,
                                                                 pars=c("beta_home","beta_away"))
posterior_summary_betas_zdts_skills<-round(
  cbind(full_zdts_skills_c_1_different_thresholds_summary$summary[,c(1,3)],full_zdts_skills_c_2_different_thresholds_summary$summary[,c(1,3)],
        full_zdts_skills_c_5_different_thresholds_summary$summary[,c(1,3)],full_zdts_skills_c_10_different_thresholds_summary$summary[,c(1,3)]),2)
xtable(posterior_summary_betas_zdts_skills)



##4) Count how many times exceeded these values
lambda1_star_full_zdts_skills_c_1_different_thresholds<-extract(full_zdts_skills_c_1_different_thresholds,pars="lambda1_star")
lambda1_star_full_zdts_skills_c_1_different_thresholds$lambda1_star

lambda2_star_full_zdts_skills_c_1_different_thresholds<-extract(full_zdts_skills_c_1_different_thresholds,pars="lambda2_star")
lambda2_star_full_zdts_skills_c_1_different_thresholds$lambda2_star

lambda1_star_full_zdts_skills_c_2_different_thresholds<-extract(full_zdts_skills_c_2_different_thresholds,pars="lambda1_star")
lambda1_star_full_zdts_skills_c_2_different_thresholds$lambda1_star

lambda2_star_full_zdts_skills_c_2_different_thresholds<-extract(full_zdts_skills_c_2_different_thresholds,pars="lambda2_star")
lambda2_star_full_zdts_skills_c_2_different_thresholds$lambda2_star

lambda1_star_full_zdts_skills_c_5_different_thresholds<-extract(full_zdts_skills_c_5_different_thresholds,pars="lambda1_star")
lambda1_star_full_zdts_skills_c_5_different_thresholds$lambda1_star

lambda2_star_full_zdts_skills_c_5_different_thresholds<-extract(full_zdts_skills_c_5_different_thresholds,pars="lambda2_star")
lambda2_star_full_zdts_skills_c_5_different_thresholds$lambda2_star

lambda1_star_full_zdts_skills_c_10_different_thresholds<-extract(full_zdts_skills_c_10_different_thresholds,pars="lambda1_star")
lambda1_star_full_zdts_skills_c_10_different_thresholds$lambda1_star

lambda2_star_full_zdts_skills_c_10_different_thresholds<-extract(full_zdts_skills_c_10_different_thresholds,pars="lambda2_star")
lambda2_star_full_zdts_skills_c_10_different_thresholds$lambda2_star


# c=100 threshold
sum(lambda1_star_full_zdts_skills_c_1_different_thresholds$lambda1_star>100)/
  length(lambda1_star_full_zdts_skills_c_1_different_thresholds$lambda1_star)#21.9%
sum(lambda2_star_full_zdts_skills_c_1_different_thresholds$lambda2_star>100)/
  length(lambda2_star_full_zdts_skills_c_1_different_thresholds$lambda2_star)#14.2%
# c=200 threshold
sum(lambda1_star_full_zdts_skills_c_2_different_thresholds$lambda1_star>2*100)/
  length(lambda1_star_full_zdts_skills_c_2_different_thresholds$lambda1_star)#19.5%
sum(lambda2_star_full_zdts_skills_c_2_different_thresholds$lambda2_star>2*100)/
  length(lambda2_star_full_zdts_skills_c_2_different_thresholds$lambda2_star)#12.2
# c=500 threshold
sum(lambda1_star_full_zdts_skills_c_5_different_thresholds$lambda1_star>5*100)/
  length(lambda1_star_full_zdts_skills_c_5_different_thresholds$lambda1_star)#16.5%
sum(lambda2_star_full_zdts_skills_c_5_different_thresholds$lambda2_star>5*100)/
  length(lambda2_star_full_zdts_skills_c_5_different_thresholds$lambda2_star)#9.9
# c=1000 threshold
sum(lambda1_star_full_zdts_skills_c_10_different_thresholds$lambda1_star>10*100)/
  length(lambda1_star_full_zdts_skills_c_10_different_thresholds$lambda1_star)#13.7%
sum(lambda2_star_full_zdts_skills_c_10_different_thresholds$lambda2_star>10*100)/
  length(lambda2_star_full_zdts_skills_c_10_different_thresholds$lambda2_star)#8.3
##---- Convertion to array (necessary for the summaries)
array_posterior_full_zdts_skills_c_2_different_thresholds<-as.array(full_zdts_skills_c_2_different_thresholds)
array_posterior_full_zdts_skills_c_5_different_thresholds<-as.array(full_zdts_skills_c_5_different_thresholds)
array_posterior_full_zdts_skills_c_10_different_thresholds<-as.array(full_zdts_skills_c_10_different_thresholds)





### 95% Posterior Density Areas for skills parameters across all candidate models (c=4,25,100 in variances)
plot_beta_home_posterior_full_zdts_skills_c_2_different_thresholds<-mcmc_areas(beta_home_full_zdts_skills_c_2_different_thresholds,
                                                                                    prob = 0.95,prob_outer=0.95,
                                                                                    point_est = c( "mean"))+ggtitle("betas_home_c_2")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))


plot_beta_away_posterior_full_zdts_skills_c_2_different_thresholds<-mcmc_areas(beta_away_full_zdts_skills_c_2_different_thresholds,
                                                                                    prob = 0.95,prob_outer=0.95,
                                                                                    point_est = c( "mean"))+ggtitle("betas_away_c_2")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))

plot_beta_home_posterior_full_zdts_skills_c_5_different_thresholds<-mcmc_areas(beta_home_full_zdts_skills_c_5_different_thresholds,
                                                                                    prob = 0.95,prob_outer=0.95,
                                                                                    point_est = c( "mean"))+ggtitle("betas_home_c_5")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))


plot_beta_away_posterior_full_zdts_skills_c_5_different_thresholds<-mcmc_areas(beta_away_full_zdts_skills_c_5_different_thresholds,
                                                                                    prob = 0.95,prob_outer=0.95,
                                                                                    point_est = c( "mean"))+ggtitle("betas_away_c_5")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))



plot_beta_home_posterior_full_zdts_skills_c_10_different_thresholds<-mcmc_areas(beta_home_full_zdts_skills_c_10_different_thresholds,
                                                                                     prob = 0.95,prob_outer=0.95,
                                                                                     point_est = c( "mean"))+ggtitle("betas_home_c_10")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))


plot_beta_away_posterior_full_zdts_skills_c_10_different_thresholds<-mcmc_areas(beta_away_full_zdts_skills_c_10_different_thresholds,
                                                                                     prob = 0.95,prob_outer=0.95,
                                                                                     point_est = c( "mean"))+ggtitle("betas_away_c_10")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))

#----Plots Combinations
full_zdts_skills_c_2_different_thresholds_5_areas<-ggarrange(plot_beta_home_posterior_full_zdts_skills_c_2_different_thresholds, plot_beta_away_posterior_full_zdts_skills_c_2_different_thresholds,
                                                                  plot_beta_home_posterior_full_zdts_skills_c_5_different_thresholds, plot_beta_away_posterior_full_zdts_skills_c_5_different_thresholds,
                                                                  ncol = 2, nrow = 2)
full_zdts_skills_c_5_different_thresholds_10_areas<-ggarrange(plot_beta_home_posterior_full_zdts_skills_c_5_different_thresholds, plot_beta_away_posterior_full_zdts_skills_c_5_different_thresholds,
                                                                   plot_beta_home_posterior_full_zdts_skills_c_10_different_thresholds, plot_beta_away_posterior_full_zdts_skills_c_10_different_thresholds,
                                                                   ncol = 2, nrow = 2)












###-------------------------Leonardo Plots




#-------------------------------------------------------
#---Autocorrelation, Trace, Cumsum plots (Not ready yet)
##-------------------------------------------------------

# MCMC Convergence diagnostics
# a) Firstly, for gammas and betas indicators
# 
# convert them to a mcmc object in terms of our convenience


mcmc_beta_home_full_zdts_skills_c_2_different_thresholds<-as.mcmc(beta_home_full_zdts_skills_c_2_different_thresholds)
mcmc_beta_away_full_zdts_skills_c_2_different_thresholds<-as.mcmc(beta_away_full_zdts_skills_c_2_different_thresholds)


##--- MCMC Convergence Diagnostics
##--- Autocorrelation, Trace and Cumsum Plots
autocorr.plot(mcmc_beta_home_full_zdts_skills_c_2_different_thresholds)
autocorr.plot(mcmc_beta_away_full_zdts_skills_c_2_different_thresholds)


traceplot(mcmc_final_posterior_values_betas_home)
traceplot(mcmc_beta_away_full_zdts_skills_c_2_different_thresholds)


cumsumplot(mcmc_beta_home_full_zdts_skills_c_2_different_thresholds)
cumsumplot(mcmc_beta_away_full_zdts_skills_c_2_different_thresholds)


