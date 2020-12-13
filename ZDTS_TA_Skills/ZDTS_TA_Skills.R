# Load the proper libraries.
library(rstan)
library(coda)
library(shinystan)
# Choose the working directory of this file (...\\Submitted_Appendix\\ZDTS\\)
setwd("C:/Users/vasileios palaskas/Desktop/BVS_Paper/ZDTS_TA_Skills")

# Load the properly prepared data for both home and away skill events as well as
# both home and away teams in each match
load("X_home")
load("X_away")
load("data_zdts_skills")


#Model matrices for home and away sets scored, respectively


#Rename the columns
colnames(X_home)<-c("(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass","
                                 (Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
                    "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
                    "(Home) block net violation","(Home) failed block","(Home) failed setting")

colnames(X_away)<-c("(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass","
                                 (Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
                    "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
                    "(Away) block net violation","(Away) failed block","(Away) failed setting")

#### Standardization of the Model Matrices for numerical convenience
X_home_std<-X_away_std<-matrix(NA,nrow=132,ncol=17)
for (i in 1:dim(X_home)[2]){
  X_home_std[,i]<-(X_home[,i]-mean(X_home[,i]))/sd(X_home[,i])
  X_away_std[,i]<-(X_away[,i]-mean(X_away[,i]))/sd(X_away[,i])
}

colnames(X_home_std)<-c("(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass","
                                 (Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
                        "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
                        "(Home) block net violation","(Home) failed block","(Home) failed setting")

colnames(X_away_std)<-c("(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass","
                                 (Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
                        "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
                        "(Away) block net violation","(Away) failed block","(Away) failed setting")


data_zdts_skills<-list(n_games=132,
		away_team=as.numeric(data_zdts_skills$away_team),
			home_team=as.numeric(data_zdts_skills$home_team),
				n_teams=data_zdts_skills$n_teams,
				X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
		home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)


## Run full_zdts_only_skills.stan
full_zdts_skills<-stan("full_zdts_skills.stan",
                             data=data_zdts_skills,chains=1,init_r=0.5,
                             iter=12000,warmup=2000)### R

# Load the output from the full ZDTS model ("full_zdts_skills")
load("full_zdts_skills")

# Extract the posterior summary statistics of both candidate variables' parameters and rest of other parameters.

betas_summary<-summary(full_zdts_skills, pars = c("beta_home","beta_away"))$summary
attack_raw_summary<-summary(full_zdts_skills, pars = c("attack_raw"))$summary
defense_raw_summary<-summary(full_zdts_skills, pars = c("defense_raw"))$summary
mu_summary<-summary(full_zdts_skills, pars = c("mu"))$summary
home_summary<-summary(full_zdts_skills, pars = c("home"))$summary


# Use their posterior means and standard deviations for both initial values specification and prior specification.

post_mean_beta_home<-betas_summary[1:17,1]####posterior mean for beta home
post_mean_beta_away<-betas_summary[18:34,1]####posterior mean for beta

post_sd_beta_home<-betas_summary[1:17,3]### posterior sd for beta
post_sd_beta_away<-betas_summary[18:34,3]### posterior sd for beta

post_mean_beta<-c(post_mean_beta_home,post_mean_beta_away)
post_sd_beta<-c(post_sd_beta_home,post_sd_beta_away)

# Step 1: Initialization of the model parameters.

gammas_home<-rep(1,17)
gammas_away<-rep(1,17)
mu<-mu_summary[,1]
home<-home_summary[,1]
attack_raw<-attack_raw_summary[,1]
defense_raw<-defense_raw_summary[,1]
betas_home<-post_mean_beta_home
betas_away<-post_mean_beta_away


# Prepare the vectors with the posterior samples of all gammas and betas, respectively.
gammas_home_matrix<-gammas_away_matrix<-betas_home_matrix<-betas_away_matrix<-NULL


T<-60000 # Total MCMC iterations

# Step 2 
for (i in 1:T){
  print(i)
  
  # Step 3: Data input needed for running the model through RStan.
  data_varsel_zdts<-list(n_teams=12,n_games=132,
                         home_sets=data_zdts_skills$home_sets,
                         away_sets=data_zdts_skills$away_sets,
                         X_home=as.matrix(X_home_std),X_away=as.matrix(X_away_std),
                         K=ncol(X_home_std),
                         home_team=data_zdts_skills$home_team,
                         away_team=data_zdts_skills$away_team,
                         gammas_home=gammas_home,gammas_away=gammas_away,
                         post_mean_beta_home=post_mean_beta_home,post_mean_beta_away=post_mean_beta_away,
                         post_sd_beta_home=post_sd_beta_home,post_sd_beta_away=post_sd_beta_away)
  
   # Step 4:Run the model through RStan for one sampling iteration (20 warm up and 21 total iterations, 21-20=1 sampling iteration) in order to update the betas from the full conditional posterior distributions. 
  # Use the previous iteration's parameter values as initial parameter values so MCMC Algorithm can begin.

  zdts_volley_skills_all<-stan("ZDTS_BVS_TA_Skills.stan",
                               data=data_varsel_zdts,chains=1,
                                iter=21,warmup=20,init=list(list(beta_home=betas_home,beta_away=betas_away,
                                                     mu=mu,
                                                     home=home,
                                                     attack_raw=attack_raw,defense_raw=defense_raw)),
control=list(adapt_window=15,adapt_init_buffer=3,adapt_term_buffer=2))### R

  # Initialize the log-likelihood for both cases 0 and 1 for gammas indicators/coefficients.
  log_point_zero<-matrix(NA,nrow=data_varsel_zdts$n_games,ncol=(2*data_varsel_zdts$K)) # matrix with log likelihoods when gamma[j]=0
  log_point_one<-matrix(NA,nrow=data_varsel_zdts$n_games,ncol=(2*data_varsel_zdts$K))  #  matrix with log likelihoods when gamma[j]=1
  
  # Extract both model's parameters and log-likelihoods for both cases of gammas indicators.
  par<-extract(zdts_volley_skills_all)
  mu<-par$mu[1]
  home<-par$home[1]
  betas_home<-par$beta_home[1,]
  betas_away<-par$beta_away[1,]
  attack_raw<-par$attack_raw[1,]
  defense_raw<-par$defense_raw[1,]
  log_point_zero<-par$log_lik_zero[1,,]
  log_point_one<-par$log_lik_one[1,,]
  
  
  # Prepare the vector of gammas for all candidate variables.
  gammas_home<-gammas_away<-NULL
  
  
  # Step 5
  for (j in 1:(2*data_varsel_zdts$K)){# K candidate variables
    log_O_j<- O_j<-NULL
    # Step 6: Calculation of the logarithm (more convenient) of O_j quantity needed in order to update the gamma indicators.
    if (j<(data_varsel_zdts$K+1)){
      log_O_j<-sum(log_point_one[,j])-sum(log_point_zero[,j])+dnorm(betas_home[j],0,sqrt(data_varsel_zdts$n_games)*post_sd_beta_home[j],log=T)-dnorm(betas_home[j],post_mean_beta_home[j],post_sd_beta_home[j],log=T)
      
    # We specify an upper bound in order to avoid Nan potential problems due to overflow.
      log_O_j[log_O_j>700]<-700
      O_j<-exp(log_O_j)
      
      gammas_home<-c(gammas_home,rbinom(1,1,O_j/(1+O_j)))
      
      
    }else {
      # Step 6: Calculation of the logarithm (more convenient) of O_j quantity needed in order to update the gamma indicators.
      log_O_j<-sum(log_point_one[,j])-sum(log_point_zero[,j])+dnorm(betas_away[j-data_varsel_zdts$K],0,sqrt(data_varsel_zdts$n_games)*post_sd_beta_away[j-data_varsel_zdts$K],log=T)-dnorm(betas_away[j-data_varsel_zdts$K],post_mean_beta_away[j-data_varsel_zdts$K],post_sd_beta_away[j-data_varsel_zdts$K],log=T)
      
    # We specify an upper bound in order to avoid Nan potential problems due to overflow.
      log_O_j[log_O_j>700]<-700
      O_j<-exp(log_O_j)
      
      gammas_away<-c(gammas_away,rbinom(1,1,O_j/(1+O_j)))
      
      
    }
  }
  
  # Step 7: In each one of T iterations, store the values of both posterior gammas and betas values in the Txp matrices
  gammas_home_matrix<-c(gammas_home_matrix,gammas_home)
  betas_home_matrix<-c(betas_home_matrix,betas_home)
  gammas_away_matrix<-c(gammas_away_matrix,gammas_away)
  betas_away_matrix<-c(betas_away_matrix,betas_away)
  
}  

# Save these values in order to manipulate them in terms of convergence diagnostics,, posterior summary statistics, etc...
save(gammas_home_matrix,file="gammas_home_matrix")
save(gammas_away_matrix,file="gammas_away_matrix")

save(betas_home_matrix,file="betas_home_matrix")
save(betas_away_matrix,file="betas_away_matrix")


# Store both gammas and betas posterior values after discarding the warmup from T iterations (here, we have chosen 10% of total T iterations).
warmup<-6000
# Each column includes the gammas values of each candidate variable.
final_posterior_values_gammas_home<-matrix(gammas_home_matrix[(data_varsel_zdts$K*warmup+1):length(gammas_home_matrix)],
                                           nrow=T-warmup,ncol=data_varsel_zdts$K,byrow=TRUE)
# Each column includes the gammas values of each candidate variable.
final_posterior_values_gammas_away<-matrix(gammas_away_matrix[(data_varsel_zdts$K*warmup+1):length(gammas_away_matrix)],
                                           nrow=T-warmup,ncol=data_varsel_zdts$K,byrow=TRUE)

final_posterior_values_betas_home<-matrix(betas_home_matrix[(data_varsel_zdts$K*warmup+1):length(betas_home_matrix)],
                                          nrow=T-warmup,ncol=data_varsel_zdts$K,byrow=TRUE)
final_posterior_values_betas_away<-matrix(betas_away_matrix[(data_varsel_zdts$K*warmup+1):length(betas_away_matrix)],
                                          nrow=T-warmup,ncol=data_varsel_zdts$K,byrow=TRUE)
# Prepare a dataframe with column names the names of candidate variables.
df_final_posterior_values_gammas_home<-as.data.frame(final_posterior_values_gammas_home)
names(X_home)<-c("(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass","
                                 (Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
                 "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
                 "(Home) block net violation","(Home) failed block","(Home) failed setting")
colnames(df_final_posterior_values_gammas_home)<-names(X_home)

df_final_posterior_values_gammas_away<-as.data.frame(final_posterior_values_gammas_away)
names(X_away)<-c("(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass","
                                 (Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
                 "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
                 "(Away) block net violation","(Away) failed block","(Away) failed setting")
colnames(df_final_posterior_values_gammas_away)<-names(X_away)

# Step 8: Obtain the posterior inclusion probabilities for each one candidate variable
posterior_inclusion_probabilities_home<-round(apply(df_final_posterior_values_gammas_home,2,mean),3)
print(posterior_inclusion_probabilities_home)

# Step 8: Obtain the posterior inclusion probabilities for each one candidate variable
posterior_inclusion_probabilities_away<-round(apply(df_final_posterior_values_gammas_away,2,mean),3)
print(posterior_inclusion_probabilities_away)

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
