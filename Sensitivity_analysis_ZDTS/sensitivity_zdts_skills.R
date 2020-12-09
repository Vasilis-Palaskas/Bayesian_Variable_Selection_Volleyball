# Load the proper libraries.
library(rstan)
library(coda)

# Choose the working directory of this file (...\\Submitted_Appendix\\ZDTS\\)
 setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/ZDTS_Skills")

# Load the properly prepared data for both home and away skill events as well as
# both home and away teams in each match
load("X_home")
load("X_away")
load("data_zdts_skills")
# 
# load("X_home")
# load("X_away")
# load("data_zdts_skills")
#Model matrices for home and away sets scored, respectively
X_home<-data_zdts_skills[,c(1:17)]
X_away<-data_zdts_skills[,c(18:34)]

#Rename the columns
colnames(X_home)<-c("home_perfect_serve","home_very_good_serve","home_failed_serve",
                    "home_perfect_pass","home_very_good_pass","home_poor_pass","home_failed_pass",
                    "home_perfect_att1","home_blocked_att1","home_failed_att1",
                    "home_perfect_att2","home_blocked_att2","home_failed_att2",
                    "home_perfect_block","home_net_violation_block","home_failed_block","home_failed_setting")

colnames(X_away)<-c("away_perfect_serve","away_very_good_serve","away_failed_serve",
                    "away_perfect_pass","away_very_good_pass","away_poor_pass","away_failed_pass",
                    "away_perfect_att1","away_blocked_att1","away_failed_att1",
                    "away_perfect_att2","away_blocked_att2","away_failed_att2",
                    "away_perfect_block","away_net_violation_block","away_failed_block","away_failed_setting")

#### Standardization of the Model Matrices for numerical convenience
X_away_std_new<-X_home_std_new<-X_home_std<-X_away_std<-matrix(NA,nrow=132,ncol=17)
# for (i in 1:dim(X_home)[2]){
#   X_home_std[,i]<-(X_home[,i]-mean(X_home[,i]))/sd(X_home[,i])
#   X_away_std[,i]<-(X_away[,i]-mean(X_away[,i]))/sd(X_away[,i])
# }
for (i in 1:dim(X_home)[2]){
  X_home_std[,i]<-scale(X_home[,i])
  X_away_std[,i]<-scale(X_away[,i])
}
colnames(X_home_std)<-c("home_perfect_serve","home_very_good_serve","home_failed_serve",
                        "home_perfect_pass","home_very_good_pass","home_poor_pass","home_failed_pass",
                        "home_perfect_att1","home_blocked_att1","home_failed_att1",
                        "home_perfect_att2","home_blocked_att2","home_failed_att2",
                        "home_perfect_block","home_net_violation_block","home_failed_block","home_failed_setting")

colnames(X_away_std)<-c("away_perfect_serve","away_very_good_serve","away_failed_serve",
                        "away_perfect_pass","away_very_good_pass","away_poor_pass","away_failed_pass",
                        "away_perfect_att1","away_blocked_att1","away_failed_att1",
                        "away_perfect_att2","away_blocked_att2","away_failed_att2",
                        "away_perfect_block","away_net_violation_block","away_failed_block","away_failed_setting")

data_zdts_only_skills<-list(n_games=data_zdts_skills$N,
                            n_teams=data_zdts_skills$n_teams,
                            X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
                            home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)

#-------------------------------------------------------------------------------------------------------------------------------
########-------------- Sensitivity Analysis--------------#################


## Different Datalists for Bayesian Models based on the value of c
## multiplied by the variance of the parameter's (mu) prior distribution

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





# Define the model implemented for sensitivity analysis
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
  real<lower=0> c;//constant for the sensitivity analysis
}

parameters {
  
  vector[K] beta_home;
  vector[K] beta_away;
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
  for (g in 1:n_games) {
    if (lambda1_star[g]>150.0){
      lambda1[g]=150.0;
    } else {
      lambda1[g]=lambda1_star[g];
    }
    if (lambda2_star[g]>150.0){
      lambda2[g]=150.0;
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
  
  //Priors
  target+=normal_lpdf(beta_home|0,c*1);
  target+=normal_lpdf(beta_away|0,c*1);
  target+=normal_lpdf(mu|0,0.37);
  target+=normal_lpdf(home|0,0.37);
  
  
  
  //likelihood-systematic component
  for (g in 1:n_games) {
    target+=skellam_without_lpmf(sets_diff[g]|lambda1[g],lambda2[g]) ;
  }
  
}
generated quantities{
  vector[n_games] log_lik;
  vector[n_games] log_lik_star;
  real dev;
  real dev_star;

  dev_star=0;
  dev=0;
    for (g in 1:n_games) {
        log_lik[g] =skellam_without_lpmf(home_sets[g]-away_sets[g]|lambda1[g],lambda2[g]) ;
        dev=dev-2*log_lik[g];
        log_lik_star[g] =skellam_without_lpmf(home_sets[g]-away_sets[g]|lambda1_star[g],lambda2_star[g]) ;
        dev_star=dev_star-2*log_lik_star[g];

  }

  //overall=attack-defense;
  //DIC=mean(dev)+0.5*variance(dev);
}


"


# Extraction of the candidate models' deviances (Table 2)

# c=2
full_zdts_only_skills_c_2<-stan(model_code=sensit_betas_zdts_skills.stan,data=data_zdts_only_skills_c_2,thin=1,chains=2,
                             iter=10000,warmup=2000,seed="12345",init_r=1)




# c=5
full_zdts_only_skills_c_5<-stan(model_code=sensit_betas_zdts_skills.stan,data=data_zdts_only_skills_c_5,thin=1,chains=2,
                                iter=10000,warmup=2000,seed="12345",init_r=1)

# c=10
full_zdts_only_skills_c_10<-stan(model_code=sensit_betas_zdts_skills.stan,data=data_zdts_only_skills_c_10,thin=1,chains=2,
                                 iter=10000,warmup=2000,seed="12345",init_r=1)




save(full_zdts_only_skills_c_2,file="full_zdts_only_skills_c_2")
save(full_zdts_only_skills_c_2,file="full_zdts_only_skills_c_5")
save(full_zdts_only_skills_c_2,file="full_zdts_only_skills_c_10")

##----------------------------------------------------------------------------------------
##-------1) Compare the deviances (both star and without star) across several values of c=2,5,10

# Extraction of the candidate models' deviances
#c=2
dev_full_zdts_only_skills_c_2<-extract(full_zdts_only_skills_c_2,pars="dev")
dev_full_zdts_only_skills_c_2$dev

dev_full_zdts_only_skills_c_2<-extract(full_zdts_only_skills_c_2,pars="dev_star")
dev_full_zdts_only_skills_c_2$dev_star

mean(dev_full_zdts_only_skills_c_2$dev)#
sd(dev_full_zdts_only_skills_c_2$dev)#
mean(dev_full_zdts_only_skills_c_2$dev_star)#
sd(dev_full_zdts_only_skills_c_2$dev_star)#
#c=5
dev_full_zdts_only_skills_c_5<-extract(full_zdts_only_skills_c_5,pars="dev")
dev_full_zdts_only_skills_c_5$dev

dev_full_zdts_only_skills_c_5<-extract(full_zdts_only_skills_c_5,pars="dev_star")
dev_full_zdts_only_skills_c_5$dev_star


mean(dev_full_zdts_only_skills_c_5$dev)#
sd(dev_full_zdts_only_skills_c_5$dev)#
mean(dev_full_zdts_only_skills_c_5$dev_star)#
sd(dev_full_zdts_only_skills_c_5$dev_star)#
#Model c=10
dev_full_zdts_only_skills_c_10<-extract(full_zdts_only_skills_c_10,pars="dev")
dev_full_zdts_only_skills_c_10$dev

dev_full_zdts_only_skills_c_10<-extract(full_zdts_only_skills_c_10,pars="dev_star")
dev_full_zdts_only_skills_c_10$dev_star

mean(dev_full_zdts_only_skills_c_10$dev)#
sd(dev_full_zdts_only_skills_c_10$dev)#
mean(dev_full_zdts_only_skills_c_10$dev_star)#
sd(dev_full_zdts_only_skills_c_10$dev_star)#
#Model c=20
dev_sensit_att_def_c_20<-extract(ZDTS_paper_att_def_c_20,pars="dev")
dev_sensit_att_def_c_20$dev

mean(dev_sensit_att_def_c_20$dev)#358.6

sd(dev_sensit_att_def_c_20$dev)#7.3


##----------------------------------------------------------------------------------------
##-------2) Compare the Posterior distributions of betas across several values of c=2,5,10
beta_home_full_zdts_only_skills_c_2<-extract(full_zdts_only_skills_c_2,pars="beta_home")
beta_away_full_zdts_only_skills_c_2<-extract(full_zdts_only_skills_c_2,pars="beta_away")

beta_home_full_zdts_only_skills_c_5<-extract(full_zdts_only_skills_c_5,pars="beta_home")
beta_away_full_zdts_only_skills_c_5<-extract(full_zdts_only_skills_c_5,pars="beta_away")

beta_home_full_zdts_only_skills_c_10<-extract(full_zdts_only_skills_c_10,pars="beta_home")
beta_away_full_zdts_only_skills_c_10<-extract(full_zdts_only_skills_c_10,pars="beta_away")


beta_home_full_zdts_only_skills_c_2<-as.data.frame(beta_home_full_zdts_only_skills_c_2)
beta_away_full_zdts_only_skills_c_2<-as.data.frame(beta_away_full_zdts_only_skills_c_2)
beta_home_full_zdts_only_skills_c_5<-as.data.frame(beta_home_full_zdts_only_skills_c_5)
beta_away_full_zdts_only_skills_c_5<-as.data.frame(beta_away_full_zdts_only_skills_c_5)
beta_home_full_zdts_only_skills_c_10<-as.data.frame(beta_home_full_zdts_only_skills_c_10)
beta_away_full_zdts_only_skills_c_10<-as.data.frame(beta_away_full_zdts_only_skills_c_10)


## Convertion to array (necessary for the summaries)
array_posterior_full_zdts_only_skills_c_2<-as.array(full_zdts_only_skills_c_2)
array_posterior_full_zdts_only_skills_c_5<-as.array(full_zdts_only_skills_c_5)
array_posterior_full_zdts_only_skills_c_10<-as.array(full_zdts_only_skills_c_10)

### 95% Posterior Intervals for ability parameters

# Figure 2
plot_beta_home_posterior_full_zdts_only_skills_c_2<-mcmc_intervals(
  array_posterior_full_zdts_only_skills_c_2[,,c(16+team_abil_order_final_ordered_logistic)],
                               prob = 0.95,prob_outer=0.95,
                              point_est = c( "mean"))+ggtitle("betas_home_c_2")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 23))

plot_beta_away_posterior_full_zdts_only_skills_c_2<-mcmc_intervals(
  array_posterior_full_zdts_only_skills_c_2[,,c(16+team_abil_order_final_ordered_logistic)],
  prob = 0.95,prob_outer=0.95,
  point_est = c( "mean"))+ggtitle("betas_away_c_2")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 23))

plot_beta_home_posterior_full_zdts_only_skills_c_2<-mcmc_intervals(beta_home_full_zdts_only_skills_c_2,
                      prob = 0.95,prob_outer=0.95,
                      point_est = c( "mean"))+ggtitle("betas_home_c_2")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))


plot_beta_away_posterior_full_zdts_only_skills_c_2<-mcmc_intervals(beta_away_full_zdts_only_skills_c_2,
                       prob = 0.95,prob_outer=0.95,
                       point_est = c( "mean"))+ggtitle("betas_away_c_2")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))


##----------------------------------------------------------------------------------------
##-------3) Count the iterations that lambda1_star,lambda2_star exceedes  the maximum therehold
##----------   across several prior variances (c=2,5,10)
sum(lambda1_star>150)/(dim(lambda1_star)[1]*dim(lambda1_star)[2])
sum(lambda2_star>150)/(dim(lambda2_star)[1]*dim(lambda2_star)[2])

###-------------------------Leonardo Plots

#----------------------------
###---------------------------------Autocorrelation, Trace, Cumsum plots (Not ready yet)
# MCMC Convergence diagnostics
# a) Firstly, for gammas and betas indicators
# 
# convert them to a mcmc pobject in terms of our convenience
mcmc_final_posterior_values_gammas_home<-as.mcmc(final_posterior_values_gammas_home)
mcmc_final_posterior_values_gammas_away<-as.mcmc(final_posterior_values_gammas_away)

mcmc_final_posterior_values_betas_home<-as.mcmc(final_posterior_values_betas_home)
mcmc_final_posterior_values_betas_away<-as.mcmc(final_posterior_values_betas_away)

autocorr.plot(mcmc_final_posterior_values_gammas_home)
autocorr.plot(mcmc_final_posterior_values_gammas_away)
autocorr.plot(mcmc_final_posterior_values_betas_home)
autocorr.plot(mcmc_final_posterior_values_betas_away)

traceplot(mcmc_final_posterior_values_gammas_home)
traceplot(mcmc_final_posterior_values_gammas_away)
traceplot(mcmc_final_posterior_values_betas_home)
traceplot(mcmc_final_posterior_values_betas_away)

cumsumplot(mcmc_final_posterior_values_gammas_home)
cumsumplot(mcmc_final_posterior_values_gammas_away)
cumsumplot(mcmc_final_posterior_values_betas_home)
cumsumplot(mcmc_final_posterior_values_betas_away)


