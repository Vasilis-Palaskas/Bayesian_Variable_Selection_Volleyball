######---------------------------------------------------------------------------- 
######---------- SECTION 3.1.2.: ZDTS Model Formulations Comparisons
# Load the proper libraries.
library(rstan)
library(coda)
library(shinystan)
library(loo)
library(bayesplot)
library(coefplot)
library(reshape2)
library(gridExtra)
library(ggmcmc)
# Choose the working directory where multiple model formulations exist
setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/ZDTS_Mod_Comp")# both home and away teams in each match


# load("data_by_sets")
# load("X_home_skills")
# load("X_away_skills")


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
X_home_diff_std<-data.frame(scale(X_home-X_away,center=T,scale=T) )
X_away_diff_std<-data.frame(scale(X_away-X_home,center=T,scale=T) )


#---Datalist required for the ZDTS with TA+Skills
#
# For numerical convenience (identifiability purposes), we use standardzed skill actions
data_separate_skills<-list(c_thres=2,c_std=5,
                          n_games=dim(data_by_sets)[1],
                          away_team=as.numeric(data_by_sets$away_Team),
                          home_team=as.numeric(data_by_sets$home_Team),
                          n_teams=length(levels(data_by_sets$home_Team)),
                          X_home=X_home_std,X_away=X_away_std,
                          K_home=ncol(X_home_std),
                          K_away=ncol(X_away_std),
                          home_sets=data_by_sets$home_sets,
                          away_sets=data_by_sets$away_sets)

data_skills_differences<-list(c_thres=2,c_std=5,
                              n_games=dim(data_by_sets)[1],
                              away_team=as.numeric(data_by_sets$away_Team),
                              home_team=as.numeric(data_by_sets$home_Team),
                              n_teams=length(levels(data_by_sets$home_Team)),
                              X_home_diff=X_home_diff_std,
                              X_away_diff=X_away_diff_std,
                              K=ncol(X_home_std),
                              home_sets=data_by_sets$home_sets,
                              away_sets=data_by_sets$away_sets)

#------ZDTS (Only skills for convenience) Model Formulations Comparisons in order to decide the best model structure for ZDTS

##1) Different betas for different model matrices corresponding to
## home and away teams skill actions, respectively (Run Model1.stan)

Model1<-stan("Model1.stan",
        data=data_separate_skills,chains=2,
                               init_r=0.5,
                               iter=5000,warmup=1000,cores=2,seed="1234")### R

save(Model1,file="Model1")

##2) Different betas for Differences between home and away skill actions corresponding to
## , respectively (Run Model1.stan)

Model2<-stan("Model2.stan",
        data=data_skills_differences,chains=2,
        init_r=0.5,
        iter=10000,warmup=1000,cores=2,seed="1234")### Bulk Effective Samples Size (ESS) is too low,

save(Model2,file="Model2")
##3) Same betas for Differences between home and away skill actions corresponding to
## , respectively (Run Model1.stan)

Model3<-stan("Model3.stan",
    data=data_skills_differences,chains=2,
    init_r=0.5,
    iter=5000,warmup=1000,cores=2,seed="1234")### R

save(Model3,file="Model3")

##4) 
model4<-stan("Model4.stan",
                   data=data_separate_skills,chains=2,
             init_r=0.5,
             iter=10000,warmup=1000,cores=2,seed="1234")### Bulk Effective Samples Size (ESS) is too low,

save(model4,file="model4")

##4) 
model5<-stan("Model5.stan",
             data=data_skills_differences,chains=2,
             init_r=0.5,
             iter=10000,warmup=1000,cores=2,seed="1234")### Bulk Effective Samples Size (ESS) is too low,

save(model5,file="model5")

###--------------Predictive Model Performance Evaluation-------------------------########
# Calculation of the DIC (Gelman,2004)
DIC_Gelman<-function(dev){
  res<-mean(dev)+0.5*var(dev)
  return(res)
}
###----Extraction of several posterior quantities (deviances, log-likelihoods)

deviance_Model1<-rstan::extract(Model1,pars="dev")$dev
log_lik_Model1<- extract_log_lik(Model1)
r_eff_log_lik_Model1<- relative_eff(exp(log_lik_Model1),
                                                       chain_id=rep(1:2,each=4000))


deviance_Model2<-rstan::extract(Model2,pars="dev")$dev
log_lik_Model2<- extract_log_lik(Model2)
r_eff_log_lik_Model2<- relative_eff(exp(log_lik_Model2),
                       chain_id=rep(1:2,each=4000))

deviance_Model3<-rstan::extract(Model3,pars="dev")$dev
log_lik_Model3<- extract_log_lik(Model3)
r_eff_log_lik_Model3<- relative_eff(exp(log_lik_Model3),
                       chain_id=rep(1:2,each=4000))

deviance_model4<-rstan::extract(model4,pars="dev")$dev
log_lik_model4<- extract_log_lik(model4)
r_eff_log_lik_model4<- relative_eff(exp(log_lik_model4),
        chain_id=rep(1:2,each=4000))

deviance_model5<-rstan::extract(model5,pars="dev")$dev
log_lik_model5<- extract_log_lik(model5)
r_eff_log_lik_model5<- relative_eff(exp(log_lik_model5),
                                    chain_id=rep(1:2,each=4000))
# ---WAIC, DIC, LOOIC Bayesian Information Criteria for the ordered logistic with only skills as covariates

waic1<-waic(log_lik_Model1)#### 248.0 
waic2<-waic(log_lik_Model2)####  298.7
waic3<-waic(log_lik_Model3)####  257.9
waic4<-waic(log_lik_model4)####   333.3
waic5<-waic(log_lik_model5)####   299.6

compare(waic1,waic2,waic3,waic4,waic5)


DIC_Gelman(deviance_Model1)# 252.5
DIC_Gelman(deviance_Model2)# 301.4
DIC_Gelman(deviance_Model3)#261.3
DIC_Gelman(deviance_model4)# 395.67
DIC_Gelman(deviance_model5)#302.44

###--------------Posterior Summary Statistics-Analysis------------------------########

#####----------------------Posterior summary----------------------------######


names(Model1)[1:16]<-colnames(X_home_std)
names(Model1)[17:32]<-colnames(X_away_std)

  


print(Model1,
      pars=c("mu","home","beta_home",
             "beta_away"),probs = c(0.025,0.5,0.975), digits=2)

## Posterior summary diagnostics

## Stan interface for both summary and convergence diagnostics
launch_shinystan(Model1)

##----Parameters Names

# teams_index <- unique(dataList$home_team)
sims <- rstan::extract(Model1)

beta_home <- sims$beta_home
beta_away <- sims$beta_away
beta_home_away<-cbind(beta_home,beta_away)
mu<- sims$mu
home<- sims$home

## Order of ability parameters (based on the posterior means)
beta_home_hat <- apply(beta_home,2, median)
beta_away_hat <- apply(beta_away,2, median)
beta_home_away_hat <- apply(beta_home_away,2, median)

mu_hat <- apply(mu,1,median)
home_hat <- apply(home,1,median)


beta_home_hat_ord <- order(beta_home_hat, decreasing = TRUE)
beta_away_hat_ord <- order(beta_away_hat, decreasing = TRUE)
beta_home_away_hat_ord <- order(beta_home_away_hat, decreasing = TRUE)

mu_hat_ord <- order(mu_hat, decreasing = TRUE)
home_hat_ord <- order(home_hat, decreasing = TRUE)

#----------ZDTS TA +Skills
##---Parameters Names

## Data frame of parameters in terms of convenience in both tables and graphs

beta_home<-data.frame(beta_home)
beta_away<-data.frame(beta_away)
beta_home_away<-data.frame(beta_home_away)

mu<-data.frame(mu)
home<-data.frame(home)
##---Proper parameters renaming
colnames(beta_home)<-colnames(X_home_std)
colnames(beta_away)<-colnames(X_away_std)

colnames(beta_home_away)<-c(colnames(X_home_std),
                            colnames(X_away_std))

###-----MCMC Posterior 95% uncertainty intervals

##---Combine 2 plots
color_scheme_set("brightblue")

pdf(file="ZDTS_Model1_skills.pdf", width =12, height =7.5)



plot1<-mcmc_intervals(beta_home[,c(beta_home_hat_ord)],
                      prob = 0.95,prob_outer=0.95,
                      point_est = c( "mean"))+ggtitle("Home Skill Events")+xlim(-6,6)+
  scale_x_continuous(breaks = seq(from = -3, to = 3, by = 1))+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))

plot2<-mcmc_intervals(beta_away[,c(beta_away_hat_ord)],
                      prob = 0.95,prob_outer=0.95,
                      point_est = c( "mean"))+ggtitle("Away Skill Events")+xlim(-6,6)+
  scale_x_continuous(breaks = seq(from = -3, to = 3, by = 1))+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))
grid.arrange(plot1,plot2,ncol=2)

dev.off()
####------------------Model 3 summary
names(Model3)[1:16]<-c("Differences_perfect_serves",
                                                           "Differences_very_good_serves","Differences_failed_serves","Differences_poor_passes",
                                                           "Differences_failed_passes","Differences_precise_passes","Differences_modetate_passes",
                                                           "Differences_perfect_att1","Differences_blocked_att1","Differences_failed_att1",
                                                           "Differences_perfect_att2","Differences_blocked_att2","Differences_failed_att2",
                                                           "Differences_failed_blocks","Differences_perfect_blocks","Differences_failed_settings")


print(Model3,
      pars=c("mu","home","beta"),probs = c(0.025,0.5,0.975), digits=2)

## Posterior summary diagnostics

## Stan interface for both summary and convergence diagnostics
launch_shinystan(Model3)




##----Parameters Names

# teams_index <- unique(dataList$home_team)
sims <- rstan::extract(Model3)
beta<- sims$beta
mu<- sims$mu
home<- sims$home


## Order of ability parameters (based on the posterior means)
beta_hat <- apply(beta,2, median)
mu_hat <- apply(mu,1,median)
home_hat <- apply(home,1,median)
beta_posit<-beta[,beta_hat>0]
beta_negat<-beta[,beta_hat<0]
beta_posit_hat <- apply(beta_posit,2, median)
beta_negat_hat <- apply(beta_negat,2, median)

beta_hat_ord <- order(beta_hat, decreasing = TRUE)
beta_posit_hat_ord <- order(beta_posit_hat, decreasing = TRUE)
beta_negat_hat_ord <- order(beta_negat_hat, decreasing = TRUE)

mu_hat_ord <- order(mu_hat, decreasing = TRUE)
home_hat_ord <- order(home_hat, decreasing = TRUE)


#----------ZDTS TA +Skills
##---Parameters Names

## Data frame of parameters in terms of convenience in both tables and graphs

beta<-data.frame(beta)

mu<-data.frame(mu)
home<-data.frame(home)
##---Proper parameters renaming
colnames(beta)<-c("perfect_serves",
                  "very_good_serves","failed_serves","poor_passes",
                  "failed_passes","precise_passes","modetate_passes",
                  "perfect_att1","blocked_att1","failed_att1",
                  "perfect_att2","blocked_att2","failed_att2",
                  "failed_blocks","perfect_blocks","failed_settings")
beta_posit<-beta[,beta_hat>0]
beta_negat<-beta[,beta_hat<0]

###-----MCMC Posterior 95% uncertainty intervals



##---Combine 2 plots
color_scheme_set("brightblue")

pdf(file="Full_ZDTS_same_betas_positive.pdf", width =12, height =7.5)

plot1<-mcmc_intervals(beta_posit[,c(beta_posit_hat_ord)],
               prob = 0.95,prob_outer=0.95,
               point_est = c( "mean"))+ggtitle("Positive Skill Events")+xlim(-6,6)+
  scale_x_continuous(breaks = seq(from = -6, to = 6, by = 1))+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))
dev.off()

pdf(file="Full_ZDTS_same_betas_negative.pdf", width =12, height =7.5)

plot2<-mcmc_intervals(beta_negat[,c(beta_negat_hat_ord)],
               prob = 0.95,prob_outer=0.95,
               point_est = c( "mean"))+ggtitle("Negative Skill Events")+xlim(-6,6)+
  scale_x_continuous(breaks = seq(from = -6, to = 6, by = 1))+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))
dev.off()

pdf(file="Full_ZDTS_same_betas.pdf", width =12, height =7.5)
grid.arrange(plot1,plot2,ncol=2)
dev.off()
