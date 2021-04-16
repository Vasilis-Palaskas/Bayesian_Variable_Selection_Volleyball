# Load the proper libraries.
library(rstan)
library(coda)
library(shinystan)
library(loo)
# Choose the working directory of this file (.../BVS_Paper/Ordered_TA_Skills)
 setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/Ordered_TA_Skills")
load("datalist_ordered")

# Choose the working directory of this file (.../BVS_Paper/Ordered_TA_Skills/4.1-... Sections)

setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/Ordered_TA_Skills/4.1 -4.3-4.4  Sections")#Rename the proper names

#Rename the proper names

names(dataList$X)<-c("perfect serve","very good serve","failed serve","perfect pass","very good pass",
                     "poor pass","failed pass","perfect att1","blocked att1",
                     "failed att1","perfect att2","blocked att2","failed att2","perfect block",
                     "block net violation","failed block","failed setting")

# Keep these skills variables with post.incl.prob.>0.5
# 

skill_events<-dataList$X

X_ordered_TA_Skills<-skill_events[,colnames(skill_events)%in%
                                    c("perfect serve","failed serve","failed pass",
                                      "perfect att1", "failed att1",
                                      "perfect att2", "failed att2",
                                      "perfect block","failed setting") ]
model_data<-list(Y=dataList$Y,X=X_ordered_TA_Skills,n_teams=12,
               home_team=as.numeric(dataList$home_team),
               away_team=as.numeric(dataList$away_team),
               N=dataList$N,K=ncol(X_ordered_TA_Skills),ncat=6)


## Table 3: Model 1
skill_events<-dataList$X

X_ordered_TA_Skills<-skill_events[,colnames(skill_events)%in%
                               c("perfect serve","failed serve",
                                 "perfect att1", "failed att1",
                                 "perfect att2", "failed att2",
                                 "perfect block","failed setting") ]
model_data<-list(Y=dataList$Y,X=X_ordered_TA_Skills,n_teams=12,
                 home_team=as.numeric(dataList$home_team),
                 away_team=as.numeric(dataList$away_team),
                 N=dataList$N,K=ncol(X_ordered_TA_Skills),ncat=6)

ordered_TA_skills_after_BVS_model1<-stan("Ordered_TA_Skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
                                      data=model_data,control=list(max_treedepth=15),cores=2)
## Table 3: Model 2

skill_events<-dataList$X

X_ordered_TA_Skills<-skill_events[,colnames(skill_events)%in%
                                      c("perfect serve","failed serve","failed pass",
                                        "perfect att1", "failed att1",
                                        "perfect att2", "failed att2",
                                        "perfect block","failed setting") ]
model_data<-list(Y=dataList$Y,X=X_ordered_TA_Skills,n_teams=12,
                 home_team=as.numeric(dataList$home_team),
                 away_team=as.numeric(dataList$away_team),
                 N=dataList$N,K=ncol(X_ordered_TA_Skills),ncat=6)

ordered_TA_skills_after_BVS_model2<-stan("Ordered_TA_Skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
                                      data=model_data,control=list(max_treedepth=15),cores=2)
## Table 3: Model 3
skill_events<-dataList$X

X_ordered_TA_Skills<-skill_events[,colnames(skill_events)%in%
                                      c("failed serve","failed pass",
                                        "perfect att1", "failed att1",
                                        "perfect att2", "failed att2",
                                        "perfect block","failed setting") ]
model_data<-list(Y=dataList$Y,X=X_ordered_TA_Skills,n_teams=12,
                 home_team=as.numeric(dataList$home_team),
                 away_team=as.numeric(dataList$away_team),
                 N=dataList$N,K=ncol(X_ordered_TA_Skills),ncat=6)

ordered_TA_skills_after_BVS_model3<-stan("Ordered_TA_Skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
                                      data=model_data,control=list(max_treedepth=15),cores=2)

## Table 3: Model 4
skill_events<-dataList$X

X_ordered_TA_Skills<-skill_events[,colnames(skill_events)%in%
                                      c("failed serve",
                                        "perfect att1", "failed att1", 
                                        "perfect att2", "failed att2",
                                        "perfect block","failed setting") ]
model_data<-list(Y=dataList$Y,X=X_ordered_TA_Skills,n_teams=12,
                 home_team=as.numeric(dataList$home_team),
                 away_team=as.numeric(dataList$away_team),
                 N=dataList$N,K=ncol(X_ordered_TA_Skills),ncat=6)

ordered_TA_skills_after_BVS_model4<-stan("Ordered_TA_Skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
                                      data=model_data,control=list(max_treedepth=15),cores=2)

##Table 3 :Information Criteria
# Calculation of the DIC (Gelman,2004)
DIC_Gelman<-function(dev){
  res<-mean(dev)+0.5*var(dev)
  return(res)
}



# Model 1
####Extraction of the log-likelihood, deviance quantities

deviance_ordered_TA_skills_after_BVS_model1<-extract(ordered_TA_skills_after_BVS_model1,pars="dev")$dev
log_lik_ordered_TA_skills_after_BVS_model1<- extract_log_lik(ordered_TA_skills_after_BVS_model1)
r_eff_log_lik_ordered_TA_skills_after_BVS_model1<- relative_eff(exp(log_lik_ordered_TA_skills_after_BVS_model1),chain_id=rep(1:2,each=6000))

##WAIC, LOO, DIC

waic(log_lik_ordered_TA_skills_after_BVS_model1)####283.4
loo(log_lik_ordered_TA_skills_after_BVS_model1)#for model with thin=1
loo(log_lik_ordered_TA_skills_after_BVS_model1,r_eff=r_eff_log_lik_final_ordered_logistic)#for model with proper thinning 379,9
DIC_Gelman(deviance_ordered_TA_skills_after_BVS_model1)#346.0

# Model 2
####Extraction of the log-likelihood, deviance quantities

deviance_ordered_TA_skills_after_BVS_model2<-extract(ordered_TA_skills_after_BVS_model2,pars="dev")$dev
log_lik_ordered_TA_skills_after_BVS_model2<- extract_log_lik(ordered_TA_skills_after_BVS_model2)
r_eff_log_lik_ordered_TA_skills_after_BVS_model2<- relative_eff(exp(log_lik_ordered_TA_skills_after_BVS_model2),chain_id=rep(1:2,each=6000))

##WAIC, LOO, DIC

waic(log_lik_ordered_TA_skills_after_BVS_model2)####283.5
loo(log_lik_ordered_TA_skills_after_BVS_model2)#for model with thin=1
loo(log_lik_ordered_TA_skills_after_BVS_model2,r_eff=r_eff_log_lik_final_ordered_logistic)#for model with proper thinning 379,9
DIC_Gelman(deviance_ordered_TA_skills_after_BVS_model2)#344.7



# Model 3
####Extraction of the log-likelihood, deviance quantities

deviance_ordered_TA_skills_after_BVS_model3<-extract(ordered_TA_skills_after_BVS_model3,pars="dev")$dev
log_lik_ordered_TA_skills_after_BVS_model3<- extract_log_lik(ordered_TA_skills_after_BVS_model3)
r_eff_log_lik_ordered_TA_skills_after_BVS_model3<- relative_eff(exp(log_lik_ordered_TA_skills_after_BVS_model3),chain_id=rep(1:2,each=6000))

##WAIC, LOO, DIC

waic(log_lik_ordered_TA_skills_after_BVS_model3)####283.2
loo(log_lik_ordered_TA_skills_after_BVS_model3)#for model with thin=1
loo(log_lik_ordered_TA_skills_after_BVS_model3,r_eff=r_eff_log_lik_final_ordered_logistic)#for model with proper thinning 379,9
DIC_Gelman(deviance_ordered_TA_skills_after_BVS_model3)#347.2


# Model 4
####Extraction of the log-likelihood, deviance quantities

deviance_ordered_TA_skills_after_BVS_model4<-extract(ordered_TA_skills_after_BVS_model4,pars="dev")$dev
log_lik_ordered_TA_skills_after_BVS_model4<- extract_log_lik(ordered_TA_skills_after_BVS_model4)
r_eff_log_lik_ordered_TA_skills_after_BVS_model4<- relative_eff(exp(log_lik_ordered_TA_skills_after_BVS_model4),chain_id=rep(1:2,each=6000))

##WAIC, LOO, DIC

waic(log_lik_ordered_TA_skills_after_BVS_model4)####296.5
loo(log_lik_ordered_TA_skills_after_BVS_model4)#for model with thin=1
loo(log_lik_ordered_TA_skills_after_BVS_model4,r_eff=r_eff_log_lik_final_ordered_logistic)#for model with proper thinning 379,9
DIC_Gelman(deviance_ordered_TA_skills_after_BVS_model4)#338.

###----DIC, WAIC Comparisons
DIC_Gelman(deviance_ordered_TA_skills_after_BVS_model1)#
DIC_Gelman(deviance_ordered_TA_skills_after_BVS_model2)#
DIC_Gelman(deviance_ordered_TA_skills_after_BVS_model3)#
DIC_Gelman(deviance_ordered_TA_skills_after_BVS_model4)#

waic(log_lik_ordered_TA_skills_after_BVS_model1)####
waic(log_lik_ordered_TA_skills_after_BVS_model2)####
waic(log_lik_ordered_TA_skills_after_BVS_model3)####
waic(log_lik_ordered_TA_skills_after_BVS_model4)####

###### SECTION 4.4 Final chosen ordered logistic model with team abilities
### MCMC Posterior Summary Plots
## Table 3: Model 3
skill_events<-dataList$X

X_ordered_TA_Skills<-skill_events[,colnames(skill_events)%in%
                                    c("failed serve","failed pass",
                                      "perfect att1", "failed att1",
                                      "perfect att2", "failed att2",
                                      "perfect block","failed setting") ]
model_data<-list(Y=dataList$Y,X=X_ordered_TA_Skills,n_teams=12,
                 home_team=as.numeric(dataList$home_team),
                 away_team=as.numeric(dataList$away_team),
                 N=dataList$N,K=ncol(X_ordered_TA_Skills),ncat=6)

ordered_TA_skills_after_BVS<-stan("Ordered_TA_Skills_after_BVS.stan",iter=14000, warmup=2000,chains=2,thin=2,
                                         data=model_data,control=list(max_treedepth=15),cores=2)
save(ordered_TA_skills_after_BVS,file="ordered_TA_skills_after_BVS")
### Bayesian Information Criteria 

####Extraction of the log-likelihood, deviance quantities

deviance_ordered_TA_skills_after_BVS<-extract(ordered_TA_skills_after_BVS,pars="dev")$dev
log_lik_ordered_TA_skills_after_BVS<- extract_log_lik(ordered_TA_skills_after_BVS)
r_eff_log_lik_ordered_TA_skills_after_BVS<- relative_eff(exp(ordered_TA_skills_after_BVS),chain_id=rep(1:2,each=6000))

##WAIC, LOO, DIC

waic(log_lik_ordered_TA_skills_after_BVS)####258.8
loo(log_lik_ordered_TA_skills_after_BVS)#260.0
loo(log_lik_ordered_skills_after_BVS,r_eff=r_eff_log_lik_ordered_skills_after_BVS)#for model with proper thinning 379,9
DIC_Gelman(deviance_ordered_TA_skills_after_BVS)#249.8

### MCMC Posterior Summary Plots

sims <- rstan::extract(ordered_TA_skills_after_BVS)

## coefplot for team abilities
teams <- c("Ethnikos Alexandroupolis", "Pamvochaikos",
           "Iraklis Petosfairishs",   "Kyzikos Peramou",   
           "Panachaiki",    "Foinikas Syrou",          
           "Kifisia",  "Orestiada",  "Olympiacos",              
           "Panathinaikos",  "Iraklis Chalkidas",  "Paok") 
teams_index <- unique(dataList$home_team)
gen_abil <- sims$gen_abil
gen_abil_hat <- apply(gen_abil,2, median)
gen_abil_sd <- apply(gen_abil,2, sd)
#teams_index <- match(squadre16_17, teams)
# att <- att[,3,teams_index]
# def <- def[,3,teams_index]
# att_hat <- apply(att,2,median)
# att_sd <- apply(att,2,sd)
# def_hat <- apply(def,2,median)
# def_sd <- apply(def,2,sd)
ord <- order(gen_abil_hat, decreasing = TRUE)
# ord_2 <- order(def_hat)

coefplot( rev(gen_abil_hat[ord]), 
          rev(gen_abil_sd[ord]), 
          CI=2, 
          varnames=rev(as.character(teams[teams_index])[ord]), 
          main="General abilities (estim. +/- 2 s.e.)\n", 
          cex.var=1.5, mar=c(1,6,4.5,1),
          cex.main=1.3,pch=16, cex=2, col="red")

## coefplot for skill events
skill_events_differences <-   c("failed serve","failed pass",
             "perfect att1", "failed att1",
             "perfect att2", "failed att2",
             "perfect block","failed setting")
# teams_index <- unique(dataList$home_team)
beta <- sims$beta
beta_hat <- apply(beta,2, median)
beta_sd <- apply(beta,2, sd)
#teams_index <- match(squadre16_17, teams)
# att <- att[,3,teams_index]
# def <- def[,3,teams_index]
# att_hat <- apply(att,2,median)
# att_sd <- apply(att,2,sd)
# def_hat <- apply(def,2,median)
# def_sd <- apply(def,2,sd)
ord <- order(beta_hat, decreasing = TRUE)
# ord_2 <- order(def_hat)

coefplot( rev(beta_hat[ord]), 
          rev(beta_sd[ord]), 
          CI=2, 
          varnames=rev(as.character(skill_events_differences)[ord]), 
          main="Skill Events Differences (estim. +/- 2 s.e.)\n", 
          cex.var=1.5, mar=c(1,6,4.5,1),
          cex.main=1.3,pch=16, cex=2, col="blue")