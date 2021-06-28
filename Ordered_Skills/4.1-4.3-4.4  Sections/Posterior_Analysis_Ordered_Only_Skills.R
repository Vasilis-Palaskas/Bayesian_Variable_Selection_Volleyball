######---------------------------------------------------------------------------- 
######---------- SECTION 4.4: Ordered logistic model with only skills (Model estimation after  BVS)
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

#---Data Preparation
source(file.choose())#-Data_Preparation.R

# Choose the working directory of this file (.../BVS_Paper/Ordered_Skills/4.1-... Sections)

setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/Ordered_Skills/4.1-4.3-4.4  Sections")
#----Rename properly the skill variables

# names(dataList$X)<-c("perfect serve","very good serve","failed serve","perfect pass","very good pass",
#                      "poor pass","failed pass","perfect att1","blocked att1",
#                      "failed att1","perfect att2","blocked att2","failed att2","perfect block",
#                      "block net violation","failed block","failed setting")
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

# Load the properly prepared data ("Data_ordered_skills").
# load("datalist_ordered")

X_home_diff<-data.frame(X_home-X_away)
colnames(X_home_diff)<-c(
  "perfect_serves","very_good_serves",
  "failed_serves","perfect_passes","very_good_passes",
  "poor_passes","failed_passes","perfect_att1",
  "blocked_att1","failed_att1","perfect_att2",
  "blocked_att2","failed_att2","perfect_blocks",
  "net_violation_blocks","failed_blocks","failed_settings")
##  Model 5: After inclusion only variables not being collinear with each other (by excluding them based on VIFs)
##  we will run a model including only skill actions emerged by the Gibbs Variable Selection implemented to those not-collinear variables
skill_events<-X_home_diff
X_ordered_Skills<-skill_events[,colnames(skill_events)%in%
                                 c("failed_serves","failed_passes",
                                   "perfect_att1", "blocked_att1","failed_att1",
                                   "perfect_att2","blocked_att2", "failed_att2",
                                   "failed_settings") ]
model_data<-list(Y=data_by_sets$sets_difference_factor,X=X_ordered_Skills,n_teams=
                   length(levels(data_by_sets$home_Team)),
                 N=dim(data_by_sets)[1],K=ncol(X_ordered_Skills),ncat=6)


ordered_skills_after_BVS_model5<-stan("ordered_skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
                                      data=model_data,control=list(max_treedepth=15),cores=2)


save(ordered_skills_after_BVS_model5,file="ordered_skills_after_BVS")

###--------------Predictive Model Performance Evaluation-------------------------########
# Calculation of the DIC (Gelman,2004)
DIC_Gelman<-function(dev){
  res<-mean(dev)+0.5*var(dev)
  return(res)
}
###----Extraction of several posterior quantities (deviances, log-likelihoods)

deviance_ordered_skills_after_BVS<-extract(ordered_skills_after_BVS,pars="dev")$dev
log_lik_ordered_skills_after_BVS<- extract_log_lik(ordered_skills_after_BVS)
r_eff_log_lik_ordered_skills_after_BVS<- relative_eff(exp(log_lik_ordered_skills_after_BVS),
                                                      chain_id=rep(1:2,each=10000))

# ---WAIC, DIC, LOOIC Bayesian Information Criteria for the ordered logistic with only skills as covariates

waic(log_lik_ordered_skills_after_BVS)####246.8 (19.9)
loo(log_lik_ordered_skills_after_BVS,
    r_eff=r_eff_log_lik_ordered_skills_after_BVS)#247.0 (19.9)for model with proper thinning 379,9
DIC_Gelman(deviance_ordered_skills_after_BVS)#245.2







######-------------ordered Logistic model with only skills
##---Parameters Names

skill_events_differences <-  c("failed_serves","failed_passes",
                               "perfect_att1", "blocked_att1","failed_att1",
                               "perfect_att2","blocked_att2", "failed_att2",
                               "failed_settings")
cutpoints<-c("c_1","c_2","c_3","c_4","c_5","c_6")

###--------------Posterior Summary Statistics-Analysis------------------------########

#####----------------------Posterior summary----------------------------######

names(ordered_skills_after_BVS_model5)[1:9]<-skill_events_differences
names(ordered_skills_after_BVS_model5)[c(10,15:19)]<-cutpoints



print(ordered_skills_after_BVS_model5,
      pars=c(
             "beta","first_temp_Intercept",
             "temp_Intercept"),probs = c(0.025,0.5,0.975), digits=2)



## Extraction of model parameters
sims <- rstan::extract(ordered_skills_after_BVS_model5)

beta <- sims$beta
first_temp_Intercept<- sims$first_temp_Intercept
temp_Intercept<- sims$temp_Intercept
temp_intercepts<-cbind(first_temp_Intercept,temp_Intercept)

## Order of ability parameters (based on the posterior means)
beta_hat <- apply(beta,2, median)
first_temp_Intercept_hat <- apply(first_temp_Intercept,1,median)
temp_intercepts_hat <- apply(temp_intercepts,2,median)

beta_hat_ord <- order(beta_hat, decreasing = TRUE)
first_temp_Intercept_hat_ord <- order(first_temp_Intercept_hat, decreasing = TRUE)
temp_intercepts_hat_ord <- order(temp_intercepts_hat, decreasing = TRUE)


##---Proper parameters renaming
colnames(beta)<-skill_events_differences
colnames(temp_intercepts)<-cutpoints



## Data frame of parameters in terms of convenience in both tables and graphs

beta<-as.data.frame(beta)
temp_intercepts<-as.data.frame(temp_intercepts)

###-----MCMC Posterior 95% uncertainty intervals


color_scheme_set("brightblue")


pdf(file="Ordered_Only_Skills_cutpoints.pdf", width =12, height =7.5)

mcmc_intervals(temp_intercepts[,c(temp_intercepts_hat_ord)],
                      prob = 0.95,prob_outer=0.95,
                      point_est = c( "mean"))+ggtitle("Cutpoints")+xlim(-6,6)+
  scale_x_continuous(breaks = seq(from = -6, to = 6, by = 1))+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))
dev.off()

pdf(file="Ordered_Only_Skills_Skills_Differences.pdf", width =12, height =7.5)

mcmc_intervals(beta[,c(beta_hat_ord)],
                      prob = 0.95,prob_outer=0.95,
                      point_est = c( "mean"))+ggtitle("Skill Differences")+xlim(-1,1)+
  scale_x_continuous(breaks = seq(from = -1, to = 1, by = 0.1))+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))

dev.off()
