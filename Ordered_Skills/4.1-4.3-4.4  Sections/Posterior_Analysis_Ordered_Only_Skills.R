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
##----Skill events selected via the BVS process based on PSI Median Threshold
skill_events<-dataList$X

X_ordered_Skills<-skill_events[,colnames(skill_events)%in%
                                 c("failed serve","failed pass",
                                   "perfect att1", "failed att1",
                                   "perfect att2", "failed att2",
                                   "perfect block","failed setting") ]
model_data<-list(Y=dataList$Y,X=X_ordered_Skills,n_teams=12,
                 N=dataList$N,K=ncol(X_ordered_Skills),ncat=6)

# Model Run via Rstan
ordered_skills_after_BVS<-stan("ordered_skills_after_BVS.stan",iter=24000, 
                               warmup=4000,chains=2,thin=2,
                               data=model_data,control=list(max_treedepth=15),cores=2)

save(ordered_skills_after_BVS,file="ordered_skills_after_BVS")

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

skill_events_differences <-   c("failed serve","failed pass",
                                "perfect att1", "failed att1",
                                "perfect att2", "failed att2",
                                "perfect block","failed setting")
cutpoints<-c("c_1","c_2","c_3","c_4","c_5","c_6")

###--------------Posterior Summary Statistics-Analysis------------------------########

#####----------------------Posterior summary----------------------------######

names(ordered_skills_after_BVS)[1:8]<-skill_events_differences
names(ordered_skills_after_BVS)[c(9,14:18)]<-cutpoints



print(ordered_skills_after_BVS,
      pars=c(
             "beta","first_temp_Intercept",
             "temp_Intercept"),probs = c(0.025,0.5,0.975), digits=2)



## Extraction of model parameters
sims <- rstan::extract(ordered_skills_after_BVS)

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
