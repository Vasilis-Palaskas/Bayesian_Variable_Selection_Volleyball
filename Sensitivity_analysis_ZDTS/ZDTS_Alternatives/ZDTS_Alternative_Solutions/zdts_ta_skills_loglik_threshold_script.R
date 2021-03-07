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
c_std<-5

data_zdts_skills_std<-list(c_std=c_std,
                                   n_games=data_zdts_skills$N,
                                   away_team=as.numeric(data_zdts_skills$away_team),
                                   home_team=as.numeric(data_zdts_skills$home_team),
                                   n_teams=data_zdts_skills$n_teams,
                                   X_home=X_home_std,X_away=X_away_std,
                                   K=ncol(X_home_std),
                                   home_sets=data_zdts_skills$home_sets,
                                   away_sets=data_zdts_skills$away_sets)


# Extraction of the candidate models' deviances (Table 2)


# ###---- ZDTS Alternative Solutions 1
full_zdts_skills<-stan(file="zdts_ta_skills_loglik_threshold.stan",
                       data=data_zdts_skills_std,thin=1,chains=1,cores=1,
                       iter=12000,warmup=2000,seed="12345",init_r=1)#15948/16000=99.6% divergent transitions after warmup

save(full_zdts_skills,file="full_zdts_skills")
##---Deviances
dev_full_zdts_skills<-extract(full_zdts_skills,pars="dev")
dev_full_zdts_skills$dev
####----- Posterior distributions of betas parameters across several c values
####-----  c: prior standard deviation multiplicator for betas parameters


min(dev_full_zdts_skills$dev)#
sd(dev_full_zdts_skills$dev)
rstan::check_divergences(full_zdts_skills)#0%

# ###-----------Rename the coefficients of stan fit objects
beta_home_full_zdts_skills<-extract(full_zdts_skills,pars="beta_home")
beta_away_full_zdts_skills<-extract(full_zdts_skills,pars="beta_away")

   colnames(beta_home_full_zdts_skills)<-c(
    "(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass",
    "(Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
     "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
     "(Home) block net violation","(Home) failed block","(Home) failed setting")

   colnames(beta_away_full_zdts_skills)<-c(
     "(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass",
     "(Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
     "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
     "(Away) block net violation","(Away) failed block","(Away) failed setting")




##---------------------------------------------------------
##-------1) Deviances Comparison across several values of c

####------Minimum Deviances


#---- 2)Sensitivity Analysis:  



#---- 3) Posterior summary statistics of betas parameters
###---
names(full_zdts_skills)[c(1:17)]<-c(
    "(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass",
    "(Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
    "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
    "(Home) block net violation","(Home) failed block","(Home) failed setting")

<-names(full_zdts_skills)[c(18:34)]<-c(
    "(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass",
    "(Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
    "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
    "(Away) block net violation","(Away) failed block","(Away) failed setting")

   #---Summary statistics Tables

full_zdts_skills_summary<-summary(full_zdts_skills,
                                  pars=c("beta_home","beta_away"))

posterior_summary_betas_zdts_skills<-round(
  cbind(full_zdts_skills_summary$summary[,c(1,3)],full_zdts_skills_summary$summary[,c(1,3)],
        full_zdts_skills_summary$summary[,c(1,3)],full_zdts_skills_summary$summary[,c(1,3)]),2)



##4) Count how many times exceeded these values
lambda1_star_full_zdts_skills<-extract(full_zdts_skills,pars="lambda1_star")
lambda1_star_full_zdts_skills$lambda1_star

lambda2_star_full_zdts_skills<-extract(full_zdts_skills,pars="lambda2_star")
lambda2_star_full_zdts_skills$lambda2_star

array_posterior_full_zdts_skills<-as.array(full_zdts_skills)

### 95% Posterior Density Areas for skills parameters across all candidate models (c=4,25,100 in variances)
plot_beta_home_posterior_full_zdts_skills<-mcmc_areas(beta_home_full_zdts_skills,
                                                      prob = 0.95,prob_outer=0.95,
                                                      point_est = c( "mean"))+ggtitle("betas_home")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))


plot_beta_away_posterior_full_zdts_skills_c_2<-mcmc_areas(beta_away_full_zdts_skills,
                                                          prob = 0.95,prob_outer=0.95,
                                                          point_est = c( "mean"))+ggtitle("betas_away_c_2")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))

full_zdts_skills_c_5_10_areas<-ggarrange(plot_beta_home_posterior_full_zdts_skills_c_5, plot_beta_away_posterior_full_zdts_skills_c_5,
                                         plot_beta_home_posterior_full_zdts_skills_c_1_20_10, plot_beta_away_posterior_full_zdts_skills_c_1_20_10,
                                         ncol = 2, nrow = 2)



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


