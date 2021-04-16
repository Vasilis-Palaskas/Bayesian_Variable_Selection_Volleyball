# Load the proper libraries.
library(rstan)
library(coda)
library(shinystan)
library(loo)
library(bayesplot)
# Choose the working directory of this file (.../BVS_Paper/Ordered_Skills)
setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/Ordered_Skills")# Load the properly full prepared data ("datalist_ordered") for the ordered logistic models.
load("datalist_ordered")

# Choose the working directory of this file (.../BVS_Paper/Ordered_Skills/4.1-... Sections)

setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/Ordered_Skills/4.1-4.3-4.4  Sections")
#----Rename properly the skill variables

names(dataList$X)<-c("perfect serve","very good serve","failed serve","perfect pass","very good pass",
                     "poor pass","failed pass","perfect att1","blocked att1",
                     "failed att1","perfect att2","blocked att2","failed att2","perfect block",
                     "block net violation","failed block","failed setting")


########----------------Table 3 Results

# Compare the different combinations of ordered logistic models with variables under consideration 
# the perfect serve and failed pass ( in total different models) includind as baseline all the other
# skill variables with PSI>0.5
skill_events<-dataList$X

X_ordered_Skills<-skill_events[,colnames(skill_events)%in%
                                  c("perfect serve","failed serve","failed pass",
                                    "perfect att1", "failed att1",
                                    "perfect att2", "failed att2",
                                    "perfect block","failed setting") ]
model_data<-list(Y=dataList$Y,X=X_ordered_Skills,n_teams=12,
               N=dataList$N,K=ncol(X_ordered_Skills),ncat=6)


## Table 3: Model 1
skill_events<-dataList$X

X_ordered_Skills<-skill_events[,colnames(skill_events)%in%
                               c("perfect serve","failed serve",
                                 "perfect att1", "failed att1",
                                 "perfect att2", "failed att2",
                                 "perfect block","failed setting") ]
model_data<-list(Y=dataList$Y,X=X_ordered_Skills,n_teams=12,
               N=dataList$N,K=ncol(X_ordered_Skills),ncat=6)

ordered_skills_after_BVS_model1<-stan("ordered_skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
                                  data=model_data,control=list(max_treedepth=15),cores=2)
## Table 3: Model 2
skill_events<-dataList$X

X_ordered_Skills<-skill_events[,colnames(skill_events)%in%
                               c("perfect serve","failed serve","failed pass",
                                 "perfect att1", "failed att1",
                                 "perfect att2", "failed att2",
                                 "perfect block","failed setting") ]

model_data<-list(Y=dataList$Y,X=X_ordered_Skills,n_teams=12,
               N=dataList$N,K=ncol(X_ordered_Skills),ncat=6)

ordered_skills_after_BVS_model2<-stan("ordered_skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
                                      data=model_data,control=list(max_treedepth=15),cores=2)

## Table 3: Model 3
skill_events<-dataList$X

X_ordered_Skills<-skill_events[,colnames(skill_events)%in%
                               c("failed serve","failed pass",
                                 "perfect att1", "failed att1",
                                 "perfect att2", "failed att2",
                                 "perfect block","failed setting") ]
model_data<-list(Y=dataList$Y,X=X_ordered_Skills,n_teams=12,
               N=dataList$N,K=ncol(X_ordered_Skills),ncat=6)

ordered_skills_after_BVS_model3<-stan("ordered_skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
                                      data=model_data,control=list(max_treedepth=15),cores=2)

## Table 3: Model 4
skill_events<-dataList$X

X_ordered_Skills<-skill_events[,colnames(skill_events)%in%
                               c("failed serve",
                                 "perfect att1", "failed att1",
                                 "perfect att2", "failed att2",
                                 "perfect block","failed setting") ]
model_data<-list(Y=dataList$Y,X=X_ordered_Skills,n_teams=12,
               N=dataList$N,K=ncol(X_ordered_Skills),ncat=6)

ordered_skills_after_BVS_model4<-stan("ordered_skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
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
deviance_ordered_skills_after_BVS_model3<-extract(ordered_skills_after_BVS_model3,pars="dev")$dev
log_lik_ordered_skills_after_BVS_model3<- extract_log_lik(ordered_skills_after_BVS_model3)
r_eff_log_lik_ordered_skills_after_BVS_model3<- relative_eff(exp(log_lik_ordered_skills_after_BVS_model3),chain_id=rep(1:2,each=6000))


# Model 4
deviance_ordered_skills_after_BVS_model4<-extract(ordered_skills_after_BVS_model4,pars="dev")$dev
log_lik_ordered_skills_after_BVS_model4<- extract_log_lik(ordered_skills_after_BVS_model4)
r_eff_log_lik_ordered_skills_after_BVS_model4<- relative_eff(exp(log_lik_ordered_skills_after_BVS_model4),chain_id=rep(1:2,each=6000))

# ---WAIC, DIC, LOOIC Criteria of Table 3 for the ordered logistic with only skills as covariates

waic(log_lik_ordered_skills_after_BVS_model1)
waic(log_lik_ordered_skills_after_BVS_model2)
waic(log_lik_ordered_skills_after_BVS_model3)
waic(log_lik_ordered_skills_after_BVS_model4)

DIC_Gelman(deviance_ordered_skills_after_BVS_model1)
DIC_Gelman(deviance_ordered_skills_after_BVS_model2)
DIC_Gelman(deviance_ordered_skills_after_BVS_model3)
DIC_Gelman(deviance_ordered_skills_after_BVS_model4)

loo(log_lik_ordered_skills_after_BVS_model1,r_eff=r_eff_log_lik_ordered_skills_after_BVS_model1)#for model with proper thinning 379,9
loo(log_lik_ordered_skills_after_BVS_model2,r_eff=r_eff_log_lik_ordered_skills_after_BVS_model2)#for model with proper thinning 379,9
loo(log_lik_ordered_skills_after_BVS_model3,r_eff=r_eff_log_lik_ordered_skills_after_BVS_model3)#for model with proper thinning 379,9
loo(log_lik_ordered_skills_after_BVS_model4,r_eff=r_eff_log_lik_ordered_skills_after_BVS_model4)#for model with proper thinning 379,9

