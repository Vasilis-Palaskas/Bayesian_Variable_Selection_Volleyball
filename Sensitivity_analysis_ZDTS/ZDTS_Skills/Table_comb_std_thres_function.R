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
c_thres<-c(1/20,1/10,1/2,1,2,5,8,9,10)
c_std<-c(1/20,1/10,1/2,1,2,5,8,9,10)


##---Initializing Matrix for minimum deviances-divergent transitions

min_dev_vector<-div_trans_vector<-NULL
for (i in c_thres){
  for (j in c_std){
    data_std_thres_zdts_skills<-list(c_thres=i,c_std=j,
                                     n_games=data_zdts_skills$N,
                                     n_teams=data_zdts_skills$n_teams,
                                     X_home=X_home_std,X_away=X_away_std,
                                     K=ncol(X_home_std),
                                     home_sets=data_zdts_skills$home_sets,
                                     away_sets=data_zdts_skills$away_sets)
    
    # Define the model implemented for sensitivity analysis
    sensit_betas_zdts_only_skills_std_thres.stan=
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

"

#---Bayesian ZDTS Models Running
warmup_iters<-2000
total_iters<-12000
full_zdts_only_skills<-stan(model_code=sensit_betas_zdts_only_skills_std_thres.stan,
                       data= data_std_thres_zdts_skills,thin=1,
                       chains=1,cores=1,
                       iter=total_iters,warmup=warmup_iters,
                       seed="1234",init_r=1)

# Extraction of the candidate models' deviances (Table 2)
dev_full_zdts_skills<-extract(full_zdts_only_skills,pars="dev")

##--Summary Statistics of Model Deviances-divergent transitions
min_dev_vector<-c(min_dev_vector,min(dev_full_zdts_skills$dev))
divergent <- get_sampler_params(full_zdts_only_skills, inc_warmup=FALSE)[[1]][,'divergent__']

div_trans_vector<-c(div_trans_vector,sum(divergent)/(total_iters-warmup_iters))

  }
}
#---Tables with minimum deviances-divergent transitions
table_min_dev<-matrix(min_dev_vector,ncol=9,nrow=9)
table_div_trans<-matrix(div_trans_vector,nrow=9,ncol=9)
table_min_dev<-round(table_min_dev,1)
table_div_trans<-round(table_div_trans,2)

write.csv(table_min_dev,file="table_deviances_zdts_skills.csv")
write.csv(table_div_trans,file="table_div_trans.csv")
library(xtable)
table_min_dev_new<-cbind(c(1/20,1/10,1/2,1,2,5,10),table_min_dev)
table_min_dev_new<-rbind(c(0,1/2,1,2,5,10),table_min_dev_new)
xtable(table_min_dev_new)
#-------------------------------------------------------
#---Autocorrelation, Trace, Cumsum plots (Not ready yet)
##-------------------------------------------------------


# MCMC Convergence diagnostics
# a) Firstly, for gammas and betas indicators
# 
# convert them to a mcmc object in terms of our convenience


mcmc_beta_home_full_zdts_skills<-as.mcmc(beta_home_full_zdts_skills)
mcmc_beta_away_full_zdts_skills<-as.mcmc(beta_away_full_zdts_skills)


##--- MCMC Convergence Diagnostics
##--- Autocorrelation, Trace and Cumsum Plots
autocorr.plot(mcmc_beta_home_full_zdts_skills)
autocorr.plot(mcmc_beta_away_full_zdts_skills)


traceplot(mcmc_final_posterior_values_betas_home)
traceplot(mcmc_beta_away_full_zdts_skills)


cumsumplot(mcmc_beta_home_full_zdts_skills)
cumsumplot(mcmc_beta_away_full_zdts_skills)


