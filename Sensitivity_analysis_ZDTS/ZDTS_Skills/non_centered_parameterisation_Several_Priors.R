# Load the proper libraries.
library(rstan)
library(coda)
library(bayesplot)
# Choose the working directory of this file (...\\Submitted_Appendix\\ZDTS\\)
setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/ZDTS_Skills")

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
" (Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
                        "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
                        "(Home) block net violation","(Home) failed block","(Home) failed setting")

colnames(X_away_std)<-c("(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass",
" (Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
                        "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
                        "(Away) block net violation","(Away) failed block","(Away) failed setting")


#------------------------------------------
########-------------- Sensitivity Analysis


###-----Datalists required for the Bayesian model fitting across several values of c
###----- c: prior standard deviation multiplicator for betas parameters

data_zdts_only_skills_c_1<-list(c=1,n_games=data_zdts_skills$N,
                                n_teams=data_zdts_skills$n_teams,
                                X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
                                home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)

data_zdts_only_skills_c_2<-list(c=2,n_games=data_zdts_skills$N,
                                n_teams=data_zdts_skills$n_teams,
                                X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
                                home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)


data_zdts_only_skills_c_5<-list(c=5,n_games=data_zdts_skills$N,
                                n_teams=data_zdts_skills$n_teams,
                                X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
                                home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)

data_zdts_only_skills_c_10<-list(c=10,n_games=data_zdts_skills$N,
                                 n_teams=data_zdts_skills$n_teams,
                                 X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
                                 home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)


###------Future Meeting 8 Action a)


data_zdts_only_skills_c_1_20<-list(c=1/20,n_games=data_zdts_skills$N,
                                   n_teams=data_zdts_skills$n_teams,
                                   X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
                                   home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)

data_zdts_only_skills_c_1_10<-list(c=1/10,n_games=data_zdts_skills$N,
                                   n_teams=data_zdts_skills$n_teams,
                                   X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
                                   home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)

data_zdts_only_skills_c_1_2<-list(c=1/2,n_games=data_zdts_skills$N,
                                  n_teams=data_zdts_skills$n_teams,
                                  X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
                                  home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)

#------- Define the model implemented for sensitivity analysis
sensit_betas_zdts_skills.stan=
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
  int<lower=0> c;
}

parameters {

  vector[K] beta_home_raw;
  vector[K] beta_away_raw;
  real mu;
  real home;
}

transformed parameters {
 // Non-centered parameterizations
  vector[K] beta_home =   beta_home_raw*c;
  vector[K] beta_away = beta_away_raw*c;

  vector[n_games]   lambda1_star;
  vector[n_games]   lambda2_star; 
  //vector[n_games]   lambda1;
  //vector[n_games]   lambda2;
  
  // Creation of linear predictor
  lambda1_star= exp(mu+home+X_home * beta_home);          
  lambda2_star= exp(mu+X_away * beta_away);  
  //for (g in 1:n_games) {
  //  if (lambda1_star[g]>100.0){
     // lambda1[g]=100.0;
   // } else {
    //  lambda1[g]=lambda1_star[g];
    //}
    //if (lambda2_star[g]>100.0){
    //  lambda2[g]=100.0;
    //} else {
    //  lambda2[g]=lambda2_star[g];
    //}
  //}
  
}

model {
  
  int    sets_diff[n_games];
  for (g in 1:n_games){
    sets_diff[g]=home_sets[g]-away_sets[g];
  }
  
  //Priors
  
  //Priors
  target +=normal_lpdf(beta_home_raw | 0, 1);
  target += normal_lpdf(beta_away_raw | 0, 1);
  target+=normal_lpdf(mu|0,0.37);
  target+=normal_lpdf(home|0,0.37);
  target+=normal_lpdf(home|0,0.37);

  //likelihood-systematic component
  for (g in 1:n_games) {
    target+=skellam_without_lpmf(sets_diff[g]|lambda1_star[g],lambda2_star[g]) ;
  }
  
}


generated quantities{
  vector[n_games] log_lik;
  real dev;
  // vector[n_teams]   overall;// overall ability

  //real DIC;

  dev=0;
  for (i in 1:n_games) {
    log_lik[i] = skellam_without_lpmf(home_sets[i]-away_sets[i] |lambda1_star[i],lambda2_star[i]);
    dev=dev-2*log_lik[i];
  }
   //overall=attack-defense;
  //DIC=mean(dev)+0.5*variance(dev);
}



"


# Extraction of the candidate models' deviances (Table 2)

# c=1
full_zdts_only_skills_c_1<-stan(model_code=sensit_betas_zdts_skills.stan,data=data_zdts_only_skills_c_1,
                                thin=1,chains=1,
                                iter=10000,warmup=2000,seed="12345",init_r=1)#15948/16000=99.6% divergent transitions after warmup



# c=2
full_zdts_only_skills_c_2<-stan(model_code=sensit_betas_zdts_skills.stan,
                                data=data_zdts_only_skills_c_2,thin=1,chains=1,
                                iter=10000,warmup=2000,seed="12345",init_r=1)#15948/16000=99.6% divergent transitions after warmup




# c=5
full_zdts_only_skills_c_1_2<-stan(model_code=sensit_betas_zdts_skills.stan,
                                  data=data_zdts_only_skills_c_1_2,thin=1,chains=1,
                                  iter=10000,warmup=2000,seed="12345",init_r=1)#15980/16000=99.8%divergent transitions after warmup

# c=10
full_zdts_only_skills_c_10<-stan(model_code=sensit_betas_zdts_skills.stan,data=data_zdts_only_skills_c_10,thin=1,chains=2,
                                 iter=10000,warmup=2000,seed="12345",init_r=1)#Effective Samples Size (ESS) is too low, 

##----- Future Meeting 7 actions
# c=1/20
full_zdts_only_skills_c_1_20<-stan(model_code=sensit_betas_zdts_skills.stan,
                                   data=data_zdts_only_skills_c_1_20,thin=1,chains=2,
                                   iter=10000,warmup=2000,seed="12345",init_r=1)#Effective Samples Size (ESS) is too low, 
# c=1/10
full_zdts_only_skills_c_1_10<-stan(model_code=sensit_betas_zdts_skills.stan,
                                   data=data_zdts_only_skills_c_1_10,thin=1,chains=2,
                                   iter=10000,warmup=2000,seed="12345",init_r=1)#Effective Samples Size (ESS) is too low, 
# c=1/2
full_zdts_only_skills_c_1_2<-stan(model_code=sensit_betas_zdts_skills.stan,
                                  data=data_zdts_only_skills_c_1_2,thin=1,chains=2,
                                  iter=10000,warmup=2000,seed="12345",init_r=1)#Effect

##---Save several fitted models (c=1,2,5,10)
save(full_zdts_only_skills_c_1,file="full_zdts_only_skills_c_1")
save(full_zdts_only_skills_c_2,file="full_zdts_only_skills_c_2")
save(full_zdts_only_skills_c_5,file="full_zdts_only_skills_c_5")
save(full_zdts_only_skills_c_10,file="full_zdts_only_skills_c_10")

save(full_zdts_only_skills_c_1_20,file="full_zdts_only_skills_c_1_20")
save(full_zdts_only_skills_c_1_10,file="full_zdts_only_skills_c_1_10")
save(full_zdts_only_skills_c_1_2,file="full_zdts_only_skills_c_1_2")
# launch_shinystan(full_zdts_only_skills_c_1)# Not converged
# launch_shinystan(full_zdts_only_skills_c_2)# Not converged
# launch_shinystan(full_zdts_only_skills_c_5)# Not converged
# launch_shinystan(full_zdts_only_skills_c_10)# Not converged

####----- Posterior distributions of betas parameters across several c values
####-----  c: prior standard deviation multiplicator for betas parameters
beta_home_full_zdts_only_skills_c_1<-extract(full_zdts_only_skills_c_1,pars="beta_home")
beta_away_full_zdts_only_skills_c_1<-extract(full_zdts_only_skills_c_1,pars="beta_away")

beta_home_full_zdts_only_skills_c_2<-extract(full_zdts_only_skills_c_2,pars="beta_home")
beta_away_full_zdts_only_skills_c_2<-extract(full_zdts_only_skills_c_2,pars="beta_away")

beta_home_full_zdts_only_skills_c_5<-extract(full_zdts_only_skills_c_5,pars="beta_home")
beta_away_full_zdts_only_skills_c_5<-extract(full_zdts_only_skills_c_5,pars="beta_away")

beta_home_full_zdts_only_skills_c_10<-extract(full_zdts_only_skills_c_10,pars="beta_home")
beta_away_full_zdts_only_skills_c_10<-extract(full_zdts_only_skills_c_10,pars="beta_away")

beta_home_full_zdts_only_skills_c_1_20<-extract(full_zdts_only_skills_c_1_20,pars="beta_home")
beta_away_full_zdts_only_skills_c_1_20<-extract(full_zdts_only_skills_c_1_20,pars="beta_away")

beta_home_full_zdts_only_skills_c_1_10<-extract(full_zdts_only_skills_c_1_10,pars="beta_home")
beta_away_full_zdts_only_skills_c_1_10<-extract(full_zdts_only_skills_c_1_10,pars="beta_away")

beta_home_full_zdts_only_skills_c_1_2<-extract(full_zdts_only_skills_c_1_2,pars="beta_home")
beta_away_full_zdts_only_skills_c_1_2<-extract(full_zdts_only_skills_c_1_2,pars="beta_away")




beta_home_full_zdts_only_skills_c_1<-as.data.frame(beta_home_full_zdts_only_skills_c_1)
beta_away_full_zdts_only_skills_c_1<-as.data.frame(beta_away_full_zdts_only_skills_c_1)
beta_home_full_zdts_only_skills_c_2<-as.data.frame(beta_home_full_zdts_only_skills_c_2)
beta_away_full_zdts_only_skills_c_2<-as.data.frame(beta_away_full_zdts_only_skills_c_2)
beta_home_full_zdts_only_skills_c_5<-as.data.frame(beta_home_full_zdts_only_skills_c_5)
beta_away_full_zdts_only_skills_c_5<-as.data.frame(beta_away_full_zdts_only_skills_c_5)
beta_home_full_zdts_only_skills_c_10<-as.data.frame(beta_home_full_zdts_only_skills_c_10)
beta_away_full_zdts_only_skills_c_10<-as.data.frame(beta_away_full_zdts_only_skills_c_10)



beta_home_full_zdts_only_skills_c_1_20<-as.data.frame(beta_home_full_zdts_only_skills_c_1_20)
beta_away_full_zdts_only_skills_c_1_20<-as.data.frame(beta_away_full_zdts_only_skills_c_1_20)

beta_home_full_zdts_only_skills_c_1_10<-as.data.frame(beta_home_full_zdts_only_skills_c_1_10)
beta_away_full_zdts_only_skills_c_1_10<-as.data.frame(beta_away_full_zdts_only_skills_c_1_10)

beta_home_full_zdts_only_skills_c_1_2<-as.data.frame(beta_home_full_zdts_only_skills_c_1_2)
beta_away_full_zdts_only_skills_c_1_2<-as.data.frame(beta_away_full_zdts_only_skills_c_1_2)
###-----------Rename the coefficients of stan fit objects
colnames(beta_home_full_zdts_only_skills_c_1_20)<-
  colnames(beta_home_full_zdts_only_skills_c_1_10)<-colnames(beta_home_full_zdts_only_skills_c_1_2)<-
  colnames(beta_home_full_zdts_only_skills_c_1)<-colnames(beta_home_full_zdts_only_skills_c_2)<-
  colnames(beta_home_full_zdts_only_skills_c_5)<-
  colnames(beta_home_full_zdts_only_skills_c_10)<-c(
    "(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass",
    "(Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
    "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
    "(Home) block net violation","(Home) failed block","(Home) failed setting")

colnames(beta_away_full_zdts_only_skills_c_1)<-colnames(beta_away_full_zdts_only_skills_c_2)<-colnames(beta_away_full_zdts_only_skills_c_5)<-
  colnames(beta_away_full_zdts_only_skills_c_10)<-c(
    "(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass",
    "(Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
    "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
    "(Away) block net violation","(Away) failed block","(Away) failed setting")





##---------------------------------------------------------
##-------1) Deviances Comparison across several values of c

#----c=1/20
dev_full_zdts_only_skills_c_1_20<-extract(full_zdts_only_skills_c_1_20,pars="dev_star")
dev_full_zdts_only_skills_c_1_20$dev_star
#----c=1/10
dev_full_zdts_only_skills_c_1_10<-extract(full_zdts_only_skills_c_1_10,pars="dev_star")
dev_full_zdts_only_skills_c_1_10$dev_star
#----c=1/2
dev_full_zdts_only_skills_c_1_2<-extract(full_zdts_only_skills_c_1_2,pars="dev_star")
dev_full_zdts_only_skills_c_1_2$dev_star
#----c=1
dev_full_zdts_only_skills_c_1<-extract(full_zdts_only_skills_c_1,pars="dev_star")
dev_full_zdts_only_skills_c_1$dev_star

#----c=2
dev_full_zdts_only_skills_c_2<-extract(full_zdts_only_skills_c_2,pars="dev_star")
dev_full_zdts_only_skills_c_2$dev_star

#----c=5

dev_full_zdts_only_skills_c_5<-extract(full_zdts_only_skills_c_5,pars="dev_star")
dev_full_zdts_only_skills_c_5$dev_star


#----c=10

dev_full_zdts_only_skills_c_10<-extract(full_zdts_only_skills_c_10,pars="dev_star")
dev_full_zdts_only_skills_c_10$dev_star


#---- 1)Sensitivity Analysis: Deviances Comparison Table (c=2,5,10)
mean(dev_full_zdts_only_skills_c_1$dev_star)#204.6
min(dev_full_zdts_only_skills_c_1$dev_star)#204.6
sd(dev_full_zdts_only_skills_c_1$dev_star)#7.7

mean(dev_full_zdts_only_skills_c_2$dev_star)#206.2
sd(dev_full_zdts_only_skills_c_2$dev_star)#8.2
min(dev_full_zdts_only_skills_c_2$dev_star)#204.6

mean(dev_full_zdts_only_skills_c_5$dev_star)#207.5
sd(dev_full_zdts_only_skills_c_5$dev_star)# 8.6
mean(dev_full_zdts_only_skills_c_10$dev_star)#207.3
sd(dev_full_zdts_only_skills_c_10$dev_star)# 8.1


###-----Minimum deviances
min(dev_full_zdts_only_skills_c_1_20$dev_star)#204.6
min(dev_full_zdts_only_skills_c_1_10$dev_star)#204.6
min(dev_full_zdts_only_skills_c_1_2$dev_star)#204.6

min(dev_full_zdts_only_skills_c_1$dev_star)#204.6
min(dev_full_zdts_only_skills_c_2$dev_star)#204.6
min(dev_full_zdts_only_skills_c_5$dev_star)#204.6
min(dev_full_zdts_only_skills_c_10$dev_star)#204.6

#---- 2)Sensitivity Analysis:
rstan::check_divergences(full_zdts_only_skills_c_1_20)#99.9%
rstan::check_divergences(full_zdts_only_skills_c_1_10)#99.9%
rstan::check_divergences(full_zdts_only_skills_c_1_2)#28.6%

rstan::check_divergences(full_zdts_only_skills_c_1)#97.0%
rstan::check_divergences(full_zdts_only_skills_c_2)#99.7%
rstan::check_divergences(full_zdts_only_skills_c_5)#99.9%
rstan::check_divergences(full_zdts_only_skills_c_10)#99.9%



#---- 3) Posterior summary statistics of betas parameters
###---
names(full_zdts_only_skills_c_1)[c(1:17)]<-names(full_zdts_only_skills_c_2)[c(1:17)]<-
  names(full_zdts_only_skills_c_5)[c(1:17)]<-names(full_zdts_only_skills_c_10)[c(1:17)]<-c(
    "(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass",
    "(Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
    "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
    "(Home) block net violation","(Home) failed block","(Home) failed setting")

names(full_zdts_only_skills_c_1)[c(18:34)]<-names(full_zdts_only_skills_c_2)[c(18:34)]<-
  names(full_zdts_only_skills_c_5)[c(18:34)]<-names(full_zdts_only_skills_c_10)[c(18:34)]<-c(
    "(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass",
    "(Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
    "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
    "(Away) block net violation","(Away) failed block","(Away) failed setting")


full_zdts_only_skills_c_1_summary<-summary(full_zdts_only_skills_c_1,
                                           pars=c("beta_home","beta_away"))
full_zdts_only_skills_c_2_summary<-summary(full_zdts_only_skills_c_2,
                                           pars=c("beta_home","beta_away"))
full_zdts_only_skills_c_5_summary<-summary(full_zdts_only_skills_c_5,
                                           pars=c("beta_home","beta_away"))
full_zdts_only_skills_c_10_summary<-summary(full_zdts_only_skills_c_10,
                                            pars=c("beta_home","beta_away"))
posterior_summary_betas_zdts_only_skills<-round(
  cbind(full_zdts_only_skills_c_1_summary$summary[,c(1,3)],full_zdts_only_skills_c_2_summary$summary[,c(1,3)],
        full_zdts_only_skills_c_5_summary$summary[,c(1,3)],full_zdts_only_skills_c_10_summary$summary[,c(1,3)]),2)
xtable(posterior_summary_betas_zdts_only_skills)

##---- Convertion to array (necessary for the summaries)
array_posterior_full_zdts_only_skills_c_2<-as.array(full_zdts_only_skills_c_2)
array_posterior_full_zdts_only_skills_c_5<-as.array(full_zdts_only_skills_c_5)
array_posterior_full_zdts_only_skills_c_10<-as.array(full_zdts_only_skills_c_10)





### 95% Posterior Density Areas for skills parameters across all candidate models (c=4,25,100 in variances)
plot_beta_home_posterior_full_zdts_only_skills_c_2<-mcmc_areas(beta_home_full_zdts_only_skills_c_2,
                                                               prob = 0.95,prob_outer=0.95,
                                                               point_est = c( "mean"))+ggtitle("betas_home_c_2")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))


plot_beta_away_posterior_full_zdts_only_skills_c_2<-mcmc_areas(beta_away_full_zdts_only_skills_c_2,
                                                               prob = 0.95,prob_outer=0.95,
                                                               point_est = c( "mean"))+ggtitle("betas_away_c_2")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))

plot_beta_home_posterior_full_zdts_only_skills_c_5<-mcmc_areas(beta_home_full_zdts_only_skills_c_5,
                                                               prob = 0.95,prob_outer=0.95,
                                                               point_est = c( "mean"))+ggtitle("betas_home_c_5")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))


plot_beta_away_posterior_full_zdts_only_skills_c_5<-mcmc_areas(beta_away_full_zdts_only_skills_c_5,
                                                               prob = 0.95,prob_outer=0.95,
                                                               point_est = c( "mean"))+ggtitle("betas_away_c_5")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))



plot_beta_home_posterior_full_zdts_only_skills_c_10<-mcmc_areas(beta_home_full_zdts_only_skills_c_10,
                                                                prob = 0.95,prob_outer=0.95,
                                                                point_est = c( "mean"))+ggtitle("betas_home_c_10")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))


plot_beta_away_posterior_full_zdts_only_skills_c_10<-mcmc_areas(beta_away_full_zdts_only_skills_c_10,
                                                                prob = 0.95,prob_outer=0.95,
                                                                point_est = c( "mean"))+ggtitle("betas_away_c_10")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))

#----Plots Combinations
full_zdts_only_skills_c_2_5_areas<-ggarrange(plot_beta_home_posterior_full_zdts_only_skills_c_2, plot_beta_away_posterior_full_zdts_only_skills_c_2,
                                             plot_beta_home_posterior_full_zdts_only_skills_c_5, plot_beta_away_posterior_full_zdts_only_skills_c_5,
                                             ncol = 2, nrow = 2)
full_zdts_only_skills_c_5_10_areas<-ggarrange(plot_beta_home_posterior_full_zdts_only_skills_c_5, plot_beta_away_posterior_full_zdts_only_skills_c_5,
                                              plot_beta_home_posterior_full_zdts_only_skills_c_10, plot_beta_away_posterior_full_zdts_only_skills_c_10,
                                              ncol = 2, nrow = 2)












###-------------------------Leonardo Plots




#-------------------------------------------------------
#---Autocorrelation, Trace, Cumsum plots (Not ready yet)
##-------------------------------------------------------

# MCMC Convergence diagnostics
# a) Firstly, for gammas and betas indicators
# 
# convert them to a mcmc object in terms of our convenience


mcmc_beta_home_full_zdts_only_skills_c_2<-as.mcmc(beta_home_full_zdts_only_skills_c_2)
mcmc_beta_away_full_zdts_only_skills_c_2<-as.mcmc(beta_away_full_zdts_only_skills_c_2)


##--- MCMC Convergence Diagnostics
##--- Autocorrelation, Trace and Cumsum Plots
autocorr.plot(mcmc_beta_home_full_zdts_only_skills_c_2)
autocorr.plot(mcmc_beta_away_full_zdts_only_skills_c_2)


traceplot(mcmc_final_posterior_values_betas_home)
traceplot(mcmc_beta_away_full_zdts_only_skills_c_2)


cumsumplot(mcmc_beta_home_full_zdts_only_skills_c_2)
cumsumplot(mcmc_beta_away_full_zdts_only_skills_c_2)


