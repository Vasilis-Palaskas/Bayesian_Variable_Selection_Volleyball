# Load the proper libraries.
library(rstan)
library(coda)
library(bayesplot)

# Activate multiple cores for stan models
options(mc.cores = parallel::detectCores())# Choose the working directory of this file (...\\Submitted_Appendix\\ZDTS\\)
# setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/ZDTS_TA_Skills")
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
########-------------- Sensitivity Analysis (Future Meeting 9)


###-----Datalists required for the Bayesian model fitting across several values of c
###----- c: prior standard deviation multiplicator for betas parameters

data_zdts_skills_std_thres<-list(c_thres=1,c_std=1/20,
                       n_games=data_zdts_skills$N,
                       away_team=as.numeric(data_zdts_skills$away_team),
                       home_team=as.numeric(data_zdts_skills$home_team),
                       n_teams=data_zdts_skills$n_teams,
                       X_home=X_home_std,X_away=X_away_std,
                       K=ncol(X_home_std),
                       home_sets=data_zdts_skills$home_sets,
                       away_sets=data_zdts_skills$away_sets)

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
     if (lambda1_star[g]>(100*c_thres)){
      lambda1[g]=(100*c_thres);
    } else {
       lambda1[g]=lambda1_star[g];
    }
    if (lambda2_star[g]>(100*c_thres)){
       lambda2[g]=(100*c_thres);
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

"


# Extraction of the candidate models' deviances (Table 2)


# c_thres=1/20,c_std=10
full_zdts_skills<-stan(model_code=sensit_betas_zdts_skills_std_thres.stan,
                                 data=data_zdts_skills_std_thres,thin=1,chains=1,cores=1,
                                 iter=12000,warmup=2000,seed="12345",init_r=1)#15948/16000=99.6% divergent transitions after warmup

save(full_zdts_skills,file="full_zdts_skills")

# 
# # c=1/10
# full_zdts_skills_c_1_20_1_10<-stan(model_code=sensit_betas_zdts_skills.stan,
#                                                    data=data_zdts_skills_c_1_10,thin=1,chains=2,cores=2,
#                                                    iter=10000,warmup=2000,seed="12345",init_r=1)#15948/16000=99.6% divergent transitions after warmup
# 
# save(full_zdts_skills_c_1_20_1_10,file="full_zdts_skills_c_1_20_1_10")
# 
# # c=1/2
# full_zdts_skills_c_1_20_1_2<-stan(model_code=sensit_betas_zdts_skills.stan,
#                                                   data=data_zdts_skills_c_1_2,thin=1,chains=2,cores=2,
#                                                   iter=10000,warmup=2000,seed="12345",init_r=1)#15948/16000=99.6% divergent transitions after warmup
# 
# save(full_zdts_skills_c_1_20_1_2,file="full_zdts_skills_c_1_20_1_2")
# 
# # c=1
# full_zdts_skills_c_1_20_1<-stan(model_code=sensit_betas_zdts_skills.stan,
#                                                 data=data_zdts_skills_c_1,thin=1,chains=2,cores=2,
#                                                 iter=10000,warmup=2000,seed="12345",init_r=1)#15948/16000=99.6% divergent transitions after warmup
# 
# save(full_zdts_skills_c_1_20_1,file="full_zdts_skills_c_1_20_1")
# 
# # c=2
# full_zdts_skills_c_2<-stan(model_code=sensit_betas_zdts_skills.stan,data=data_zdts_skills_c_2,thin=1,chains=2,cores=2,
#                                                 iter=10000,warmup=2000,seed="12345",init_r=1)#15948/16000=99.6% divergent transitions after warmup
# 
# save(full_zdts_skills_c_2,file="full_zdts_skills_c_2")
# 
# # c=5
# full_zdts_skills_c_5<-stan(model_code=sensit_betas_zdts_skills.stan,data=data_zdts_skills_c_5,thin=1,chains=2,cores=2,
#                                                 iter=10000,warmup=2000,seed="12345",init_r=1)#15980/16000=99.8%divergent transitions after warmup
# save(full_zdts_skills_c_5,file="full_zdts_skills_c_5")
# 
# # c=10
# full_zdts_skills_c_1_20_10<-stan(model_code=sensit_betas_zdts_skills.stan,data=data_zdts_skills_c_10,thin=1,chains=2,cores=2,
#                                                  iter=10000,warmup=2000,seed="12345",init_r=1)#Effective Samples Size (ESS) is too low, 
# 
# save(full_zdts_skills_c_1_20_10,file="full_zdts_skills_c_1_20_10")

##---Save several fitted models (c=1,2,5,10)

# launch_shinystan(full_zdts_skills_c_1_20_1)# Not converged
# launch_shinystan(full_zdts_skills_c_2)# Not converged
# launch_shinystan(full_zdts_skills_c_5)# Not converged
# launch_shinystan(full_zdts_skills_c_1_20_10)# Not converged


####----- Posterior distributions of betas parameters across several c values
####-----  c: prior standard deviation multiplicator for betas parameters
beta_home_full_zdts_skills_c_1_20_1<-extract(full_zdts_skills_c_1_20_1,pars="beta_home")
beta_away_full_zdts_skills_c_1_20_1<-extract(full_zdts_skills_c_1_20_1,pars="beta_away")
# 
# beta_home_full_zdts_skills_c_2<-extract(full_zdts_skills_c_2,pars="beta_home")
# beta_away_full_zdts_skills_c_2<-extract(full_zdts_skills_c_2,pars="beta_away")
# 
# beta_home_full_zdts_skills_c_5<-extract(full_zdts_skills_c_5,pars="beta_home")
# beta_away_full_zdts_skills_c_5<-extract(full_zdts_skills_c_5,pars="beta_away")
# 
# beta_home_full_zdts_skills_c_1_20_10<-extract(full_zdts_skills_c_1_20_10,pars="beta_home")
# beta_away_full_zdts_skills_c_1_20_10<-extract(full_zdts_skills_c_1_20_10,pars="beta_away")
# 
# beta_home_full_zdts_skills_c_1_20_1<-as.data.frame(beta_home_full_zdts_skills_c_1_20_1)
# beta_away_full_zdts_skills_c_1_20_1<-as.data.frame(beta_away_full_zdts_skills_c_1_20_1)
# # beta_home_full_zdts_skills_c_2<-as.data.frame(beta_home_full_zdts_skills_c_2)
# # beta_away_full_zdts_skills_c_2<-as.data.frame(beta_away_full_zdts_skills_c_2)
# # beta_home_full_zdts_skills_c_5<-as.data.frame(beta_home_full_zdts_skills_c_5)
# # beta_away_full_zdts_skills_c_5<-as.data.frame(beta_away_full_zdts_skills_c_5)
# # beta_home_full_zdts_skills_c_1_20_10<-as.data.frame(beta_home_full_zdts_skills_c_1_20_10)
# # beta_away_full_zdts_skills_c_1_20_10<-as.data.frame(beta_away_full_zdts_skills_c_1_20_10)
# 
# ###-----------Rename the coefficients of stan fit objects
# colnames(beta_home_full_zdts_skills_c_1_20_1)<-colnames(beta_home_full_zdts_skills_c_2)<-colnames(beta_home_full_zdts_skills_c_5)<-
#   colnames(beta_home_full_zdts_skills_c_1_20_10)<-c(
#     "(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass",
#     "(Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
#     "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
#     "(Home) block net violation","(Home) failed block","(Home) failed setting")
# 
# colnames(beta_away_full_zdts_skills_c_1_20_1)<-colnames(beta_away_full_zdts_skills_c_2)<-colnames(beta_away_full_zdts_skills_c_5)<-
#   colnames(beta_away_full_zdts_skills_c_1_20_10)<-c(
#     "(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass",
#     "(Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
#     "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
#     "(Away) block net violation","(Away) failed block","(Away) failed setting")
# 



##---------------------------------------------------------
##-------1) Deviances Comparison across several values of c


#----c=1/20, c=1/20
dev_full_zdts_skills_c_1_20_1_20<-extract(full_zdts_skills_c_1_20_1_20,pars="dev")
dev_full_zdts_skills_c_1_20_1_20$dev

#----c=1/20, c=1/10
dev_full_zdts_skills_c_1_20_1_10<-extract(full_zdts_skills_c_1_20_1_10,pars="dev")
dev_full_zdts_skills_c_1_20_1_10$dev

#----c=1/20, c=1/2
dev_full_zdts_skills_c_1_20_1_2<-extract(full_zdts_skills_c_1_20_1_2,pars="dev")
dev_full_zdts_skills_c_1_20_1_2$dev
#----c=1/20, c=1
dev_full_zdts_skills_c_1_20_1<-extract(full_zdts_skills_c_1_20_1,pars="dev")
dev_full_zdts_skills_c_1_20_1$dev

#----c=1/20, c=2
dev_full_zdts_skills_c_1_20_2<-extract(full_zdts_skills_c_1_20_2,pars="dev")
dev_full_zdts_skills_c_1_20_2$dev
#----c=1/20, c=5
dev_full_zdts_skills_c_1_20_5<-extract(full_zdts_skills_c_1_20_5,pars="dev")
dev_full_zdts_skills_c_1_20_5$dev
#----c=1/20, c=10
dev_full_zdts_skills_c_1_20_10<-extract(full_zdts_skills_c_1_20_10,pars="dev")
dev_full_zdts_skills_c_1_20_10$dev
####------Minimum Deviances
min(dev_full_zdts_skills_c_1_20_1_20$dev)#
min(dev_full_zdts_skills_c_1_20_1_10$dev)#
min(dev_full_zdts_skills_c_1_20_1_2$dev)#
min(dev_full_zdts_skills_c_1_20_1$dev)#
min(dev_full_zdts_skills_c_1_20_2$dev)#
min(dev_full_zdts_skills_c_1_20_5$dev)#
min(dev_full_zdts_skills_c_1_20_10$dev)#

#---- 2)Sensitivity Analysis:  
rstan::check_divergences(full_zdts_skills_c_1_20_1)#0%
rstan::check_divergences(full_zdts_skills_c_1_20_1_20)#0%
rstan::check_divergences(full_zdts_skills_c_1_20_1_2)#0%

rstan::check_divergences(full_zdts_skills_c_1_20_1)#0%
rstan::check_divergences(full_zdts_skills_c_2)#0%
rstan::check_divergences(full_zdts_skills_c_5)#50.2%
rstan::check_divergences(full_zdts_skills_c_1_20_10)#86.2%



#---- 3) Posterior summary statistics of betas parameters
###---
names(full_zdts_skills_c_1_20_1)[c(1:17)]<-names(full_zdts_skills_c_2)[c(1:17)]<-
  names(full_zdts_skills_c_5)[c(1:17)]<-names(full_zdts_skills_c_1_20_10)[c(1:17)]<-c(
    "(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass",
    "(Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
    "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
    "(Home) block net violation","(Home) failed block","(Home) failed setting")

names(full_zdts_skills_c_1_20_1)[c(18:34)]<-names(full_zdts_skills_c_2)[c(18:34)]<-
  names(full_zdts_skills_c_5)[c(18:34)]<-names(full_zdts_skills_c_1_20_10)[c(18:34)]<-c(
    "(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass",
    "(Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
    "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
    "(Away) block net violation","(Away) failed block","(Away) failed setting")


full_zdts_skills_c_1_20_1_summary<-summary(full_zdts_skills_c_1_20_1,
                                                           pars=c("beta_home","beta_away"))
full_zdts_skills_c_2_summary<-summary(full_zdts_skills_c_2,
                                                           pars=c("beta_home","beta_away"))
full_zdts_skills_c_5_summary<-summary(full_zdts_skills_c_5,
                                                           pars=c("beta_home","beta_away"))
full_zdts_skills_c_1_20_10_summary<-summary(full_zdts_skills_c_1_20_10,
                                                            pars=c("beta_home","beta_away"))
posterior_summary_betas_zdts_skills<-round(
  cbind(full_zdts_skills_c_1_20_1_summary$summary[,c(1,3)],full_zdts_skills_c_2_summary$summary[,c(1,3)],
        full_zdts_skills_c_5_summary$summary[,c(1,3)],full_zdts_skills_c_1_20_10_summary$summary[,c(1,3)]),2)
xtable(posterior_summary_betas_zdts_skills)



##4) Count how many times exceeded these values
lambda1_star_full_zdts_skills_c_1_20_1<-extract(full_zdts_skills_c_1_20_1,pars="lambda1_star")
lambda1_star_full_zdts_skills_c_1_20_1$lambda1_star

lambda2_star_full_zdts_skills_c_1_20_1<-extract(full_zdts_skills_c_1_20_1,pars="lambda2_star")
lambda2_star_full_zdts_skills_c_1_20_1$lambda2_star

lambda1_star_full_zdts_skills_c_2<-extract(full_zdts_skills_c_2,pars="lambda1_star")
lambda1_star_full_zdts_skills_c_2$lambda1_star

lambda2_star_full_zdts_skills_c_2<-extract(full_zdts_skills_c_2,pars="lambda2_star")
lambda2_star_full_zdts_skills_c_2$lambda2_star

lambda1_star_full_zdts_skills_c_5<-extract(full_zdts_skills_c_5,pars="lambda1_star")
lambda1_star_full_zdts_skills_c_5$lambda1_star

lambda2_star_full_zdts_skills_c_5<-extract(full_zdts_skills_c_5,pars="lambda2_star")
lambda2_star_full_zdts_skills_c_5$lambda2_star

lambda1_star_full_zdts_skills_c_1_20_10<-extract(full_zdts_skills_c_1_20_10,pars="lambda1_star")
lambda1_star_full_zdts_skills_c_1_20_10$lambda1_star

lambda2_star_full_zdts_skills_c_1_20_10<-extract(full_zdts_skills_c_1_20_10,pars="lambda2_star")
lambda2_star_full_zdts_skills_c_1_20_10$lambda2_star


# c=100 threshold
sum(lambda1_star_full_zdts_skills_c_1_20_1$lambda1_star>100)/
  length(lambda1_star_full_zdts_skills_c_1_20_1$lambda1_star)#21.9%
sum(lambda2_star_full_zdts_skills_c_1_20_1$lambda2_star>100)/
  length(lambda2_star_full_zdts_skills_c_1_20_1$lambda2_star)#14.2%
# c=200 threshold
sum(lambda1_star_full_zdts_skills_c_2$lambda1_star>2*100)/
  length(lambda1_star_full_zdts_skills_c_2$lambda1_star)#19.5%
sum(lambda2_star_full_zdts_skills_c_2$lambda2_star>2*100)/
  length(lambda2_star_full_zdts_skills_c_2$lambda2_star)#12.2
# c=500 threshold
sum(lambda1_star_full_zdts_skills_c_5$lambda1_star>5*100)/
  length(lambda1_star_full_zdts_skills_c_5$lambda1_star)#16.5%
sum(lambda2_star_full_zdts_skills_c_5$lambda2_star>5*100)/
  length(lambda2_star_full_zdts_skills_c_5$lambda2_star)#9.9
# c=1000 threshold
sum(lambda1_star_full_zdts_skills_c_1_20_10$lambda1_star>10*100)/
  length(lambda1_star_full_zdts_skills_c_1_20_10$lambda1_star)#13.7%
sum(lambda2_star_full_zdts_skills_c_1_20_10$lambda2_star>10*100)/
  length(lambda2_star_full_zdts_skills_c_1_20_10$lambda2_star)#8.3
##---- Convertion to array (necessary for the summaries)
array_posterior_full_zdts_skills_c_2<-as.array(full_zdts_skills_c_2)
array_posterior_full_zdts_skills_c_5<-as.array(full_zdts_skills_c_5)
array_posterior_full_zdts_skills_c_1_20_10<-as.array(full_zdts_skills_c_1_20_10)





### 95% Posterior Density Areas for skills parameters across all candidate models (c=4,25,100 in variances)
plot_beta_home_posterior_full_zdts_skills_c_2<-mcmc_areas(beta_home_full_zdts_skills_c_2,
                                                                               prob = 0.95,prob_outer=0.95,
                                                                               point_est = c( "mean"))+ggtitle("betas_home_c_2")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))


plot_beta_away_posterior_full_zdts_skills_c_2<-mcmc_areas(beta_away_full_zdts_skills_c_2,
                                                                               prob = 0.95,prob_outer=0.95,
                                                                               point_est = c( "mean"))+ggtitle("betas_away_c_2")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))

plot_beta_home_posterior_full_zdts_skills_c_5<-mcmc_areas(beta_home_full_zdts_skills_c_5,
                                                                               prob = 0.95,prob_outer=0.95,
                                                                               point_est = c( "mean"))+ggtitle("betas_home_c_5")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))


plot_beta_away_posterior_full_zdts_skills_c_5<-mcmc_areas(beta_away_full_zdts_skills_c_5,
                                                                               prob = 0.95,prob_outer=0.95,
                                                                               point_est = c( "mean"))+ggtitle("betas_away_c_5")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))



plot_beta_home_posterior_full_zdts_skills_c_1_20_10<-mcmc_areas(beta_home_full_zdts_skills_c_1_20_10,
                                                                                prob = 0.95,prob_outer=0.95,
                                                                                point_est = c( "mean"))+ggtitle("betas_home_c_10")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))


plot_beta_away_posterior_full_zdts_skills_c_1_20_10<-mcmc_areas(beta_away_full_zdts_skills_c_1_20_10,
                                                                                prob = 0.95,prob_outer=0.95,
                                                                                point_est = c( "mean"))+ggtitle("betas_away_c_10")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))

#----Plots Combinations
full_zdts_skills_c_2_5_areas<-ggarrange(plot_beta_home_posterior_full_zdts_skills_c_2, plot_beta_away_posterior_full_zdts_skills_c_2,
                                                             plot_beta_home_posterior_full_zdts_skills_c_5, plot_beta_away_posterior_full_zdts_skills_c_5,
                                                             ncol = 2, nrow = 2)
full_zdts_skills_c_5_10_areas<-ggarrange(plot_beta_home_posterior_full_zdts_skills_c_5, plot_beta_away_posterior_full_zdts_skills_c_5,
                                                              plot_beta_home_posterior_full_zdts_skills_c_1_20_10, plot_beta_away_posterior_full_zdts_skills_c_1_20_10,
                                                              ncol = 2, nrow = 2)












###-------------------------Leonardo Plots




#-------------------------------------------------------
#---Autocorrelation, Trace, Cumsum plots (Not ready yet)
##-------------------------------------------------------

# MCMC Convergence diagnostics
# a) Firstly, for gammas and betas indicators
# 
# convert them to a mcmc object in terms of our convenience


mcmc_beta_home_full_zdts_skills_c_2<-as.mcmc(beta_home_full_zdts_skills_c_2)
mcmc_beta_away_full_zdts_skills_c_2<-as.mcmc(beta_away_full_zdts_skills_c_2)


##--- MCMC Convergence Diagnostics
##--- Autocorrelation, Trace and Cumsum Plots
autocorr.plot(mcmc_beta_home_full_zdts_skills_c_2)
autocorr.plot(mcmc_beta_away_full_zdts_skills_c_2)


traceplot(mcmc_final_posterior_values_betas_home)
traceplot(mcmc_beta_away_full_zdts_skills_c_2)


cumsumplot(mcmc_beta_home_full_zdts_skills_c_2)
cumsumplot(mcmc_beta_away_full_zdts_skills_c_2)


