# Load the proper libraries.
library(rstan)
library(coda)
library(shinystan)
library(loo)
library(bayesplot)

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

#---Transform set difference values in terms of fitting for ordered multinomial model (requires positive integers or factors)
data_by_sets$sets_difference_factor<-data_by_sets$sets_difference
for (i in 1:dim(data_by_sets)[1]){
  if (data_by_sets$sets_difference[i]==(-3)){
    data_by_sets$sets_difference_factor[i]<-1
  } else if (data_by_sets$sets_difference[i]==(-2)){
    data_by_sets$sets_difference_factor[i]<-2
  } else if (data_by_sets$sets_difference[i]==(-1)){
    data_by_sets$sets_difference_factor[i]<-3
  } else if (data_by_sets$sets_difference[i]==(1)){
    data_by_sets$sets_difference_factor[i]<-4
  } else if (data_by_sets$sets_difference[i]==(2)){
    data_by_sets$sets_difference_factor[i]<-5
  } else if (data_by_sets$sets_difference[i]==(3)){
    data_by_sets$sets_difference_factor[i]<-6
  }
  
}
########----------------Table 3 Results


# Compare the different combinations of ordered logistic models with variables under consideration 
# the perfect serve and failed pass ( in total different models) includind as baseline all the other
# skill variables with PSI>0.5
skill_events<-X_home_diff

X_ordered_Skills<-skill_events[,colnames(skill_events)%in%
                                  c("perfect_serves","failed_serves","failed_passes",
                                    "perfect_att1", "failed_att1",
                                    "perfect_att2", "failed_att2",
                                    "perfect_blocks","failed_settings") ]

model_data<-list(Y=data_by_sets$sets_difference_factor,X=X_ordered_Skills,n_teams=
                   length(levels(data_by_sets$home_Team)),
               N=dim(data_by_sets)[1],K=ncol(X_ordered_Skills),ncat=6)


## Table 3: Model 1
skill_events<-X_home_diff
X_ordered_Skills<-skill_events[,colnames(skill_events)%in%
                                 c("perfect_serves","failed_serves",
                                   "perfect_att1", "failed_att1",
                                   "perfect_att2", "failed_att2",
                                   "perfect_blocks","failed_settings")]
model_data<-list(Y=data_by_sets$sets_difference_factor,X=X_ordered_Skills,n_teams=
                   length(levels(data_by_sets$home_Team)),
                 N=dim(data_by_sets)[1],K=ncol(X_ordered_Skills),ncat=6)


ordered_skills_after_BVS_model1<-stan("ordered_skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
                                  data=model_data,control=list(max_treedepth=15),cores=2)
## Table 3: Model 2
skill_events<-X_home_diff
X_ordered_Skills<-skill_events[,colnames(skill_events)%in%
                                 c("perfect_serves","failed_serves","failed_passes",
                                   "perfect_att1", "failed_att1",
                                   "perfect_att2", "failed_att2",
                                   "perfect_blocks","failed_settings") ]

model_data<-list(Y=data_by_sets$sets_difference_factor,X=X_ordered_Skills,n_teams=
                   length(levels(data_by_sets$home_Team)),
                 N=dim(data_by_sets)[1],K=ncol(X_ordered_Skills),ncat=6)


ordered_skills_after_BVS_model2<-stan("ordered_skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
                                      data=model_data,control=list(max_treedepth=15),cores=2)

## Table 3: Model 3
skill_events<-X_home_diff
X_ordered_Skills<-skill_events[,colnames(skill_events)%in%
                                 c("failed_serves","failed_passes",
                                   "perfect_att1", "failed_att1",
                                   "perfect_att2", "failed_att2",
                                   "perfect_blocks","failed_settings") ]
model_data<-list(Y=data_by_sets$sets_difference_factor,X=X_ordered_Skills,n_teams=
                   length(levels(data_by_sets$home_Team)),
                 N=dim(data_by_sets)[1],K=ncol(X_ordered_Skills),ncat=6)


ordered_skills_after_BVS_model3<-stan("ordered_skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
                                      data=model_data,control=list(max_treedepth=15),cores=2)

## Table 3: Model 4
skill_events<-X_home_diff
X_ordered_Skills<-skill_events[,colnames(skill_events)%in%
                                 c("failed_serves",
                                   "perfect_att1", "failed_att1",
                                   "perfect_att2", "failed_att2",
                                   "perfect_blocks","failed_settings") ]
model_data<-list(Y=data_by_sets$sets_difference_factor,X=X_ordered_Skills,n_teams=
                   length(levels(data_by_sets$home_Team)),
                 N=dim(data_by_sets)[1],K=ncol(X_ordered_Skills),ncat=6)


ordered_skills_after_BVS_model4<-stan("ordered_skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
                                      data=model_data,control=list(max_treedepth=15),cores=2)
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


 ##Table 3 :Information Criteria
# Calculation of the DIC (Gelman,2004)
DIC_Gelman<-function(dev){
  res<-mean(dev)+0.5*var(dev)
  return(res)
}


###----Extraction of several posterior quantities (deviances, log-likelihoods)
# Model 1

 deviance_ordered_skills_after_BVS_model1<-extract(ordered_skills_after_BVS_model1,pars="dev")$dev
log_lik_ordered_skills_after_BVS_model1<- extract_log_lik(ordered_skills_after_BVS_model1)
r_eff_log_lik_ordered_skills_after_BVS_model1<- relative_eff(exp(log_lik_ordered_skills_after_BVS_model1),chain_id=rep(1:2,each=6000))


# Model 2

deviance_ordered_skills_after_BVS_model2<-extract(ordered_skills_after_BVS_model2,pars="dev")$dev
log_lik_ordered_skills_after_BVS_model2<- extract_log_lik(ordered_skills_after_BVS_model2)
r_eff_log_lik_ordered_skills_after_BVS_model2<- relative_eff(exp(log_lik_ordered_skills_after_BVS_model2),chain_id=rep(1:2,each=6000))



# Model 3
deviance_ordered_skills_after_BVS_model3<-rstan::extract(ordered_skills_after_BVS_model3,pars="dev")$dev
log_lik_ordered_skills_after_BVS_model3<-extract_log_lik(ordered_skills_after_BVS_model3)
r_eff_log_lik_ordered_skills_after_BVS_model3<- relative_eff(exp(log_lik_ordered_skills_after_BVS_model3),chain_id=rep(1:2,each=6000))


# Model 4
deviance_ordered_skills_after_BVS_model4<-extract(ordered_skills_after_BVS_model4,pars="dev")$dev
log_lik_ordered_skills_after_BVS_model4<- extract_log_lik(ordered_skills_after_BVS_model4)
r_eff_log_lik_ordered_skills_after_BVS_model4<- relative_eff(exp(log_lik_ordered_skills_after_BVS_model4),chain_id=rep(1:2,each=6000))

# Model 5
deviance_ordered_skills_after_BVS_model5<-rstan::extract(ordered_skills_after_BVS_model5,pars="dev")$dev
log_lik_ordered_skills_after_BVS_model5<- extract_log_lik(ordered_skills_after_BVS_model5)
r_eff_log_lik_ordered_skills_after_BVS_model5<- relative_eff(exp(log_lik_ordered_skills_after_BVS_model5),chain_id=rep(1:2,each=6000))

# ---WAIC, DIC, LOOIC Criteria of Table 3 for the ordered logistic with only skills as covariates

waic(log_lik_ordered_skills_after_BVS_model1)
waic(log_lik_ordered_skills_after_BVS_model2)
waic(log_lik_ordered_skills_after_BVS_model3)#247.0 (20.0)
waic(log_lik_ordered_skills_after_BVS_model4)
waic(log_lik_ordered_skills_after_BVS_model5)#    249.0 20.1

DIC_Gelman(deviance_ordered_skills_after_BVS_model1)
DIC_Gelman(deviance_ordered_skills_after_BVS_model2)
DIC_Gelman(deviance_ordered_skills_after_BVS_model3)#245.4
DIC_Gelman(deviance_ordered_skills_after_BVS_model4)
DIC_Gelman(deviance_ordered_skills_after_BVS_model5)#    249.0 20.1

loo(log_lik_ordered_skills_after_BVS_model1,r_eff=r_eff_log_lik_ordered_skills_after_BVS_model1)#for model with proper thinning 379,9
loo(log_lik_ordered_skills_after_BVS_model2,r_eff=r_eff_log_lik_ordered_skills_after_BVS_model2)#for model with proper thinning 379,9
loo(log_lik_ordered_skills_after_BVS_model3,r_eff=r_eff_log_lik_ordered_skills_after_BVS_model3)# 247.2 20.0
loo(log_lik_ordered_skills_after_BVS_model4,r_eff=r_eff_log_lik_ordered_skills_after_BVS_model4)#for model with proper thinning 379,9
loo(log_lik_ordered_skills_after_BVS_model5,r_eff=r_eff_log_lik_ordered_skills_after_BVS_model5)#for model with proper thinning 249.2 20.1

