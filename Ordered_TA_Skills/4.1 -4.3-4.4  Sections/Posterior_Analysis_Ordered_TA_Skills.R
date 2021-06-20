######---------------------------------------------------------------------------- 
######---------- SECTION 4.4: Ordered logistic model with Team Abilities (TA) as well as skills (Model estimation after  BVS)
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
# Choose the working directory of this file (.../BVS_Paper/Ordered_TA_Skills)
 setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/Ordered_TA_Skills")
load("datalist_ordered")

# Choose the working directory of this file (.../BVS_Paper/Ordered_TA_Skills/4.1-... Sections)

setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/Ordered_TA_Skills/4.1 -4.3-4.4  Sections")#Rename the proper names

#----Rename properly the skill variables

names(dataList$X)<-c("perfect serve","very good serve","failed serve","perfect pass","very good pass",
                     "poor pass","failed pass","perfect att1","blocked att1",
                     "failed att1","perfect att2","blocked att2","failed att2","perfect block",
                     "block net violation","failed block","failed setting")

##----Skill events selected via the BVS process based on PSI Median Threshold

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

ordered_TA_skills_after_BVS<-stan("Ordered_TA_Skills_after_BVS.stan",iter=28000, 
                                  warmup=8000,chains=2,thin=2,
                                  data=model_data,control=list(max_treedepth=15),cores=2)
save(ordered_TA_skills_after_BVS,file="ordered_TA_skills_after_BVS")

###--------------Predictive Model Performance Evaluation-------------------------########
# Calculation of the DIC (Gelman,2004)

DIC_Gelman<-function(dev){
  res<-mean(dev)+0.5*var(dev)
  return(res)
}
###----Extraction of several posterior quantities (deviances, log-likelihoods)

deviance_ordered_TA_skills_after_BVS<-extract(ordered_TA_skills_after_BVS,pars="dev")$dev
log_lik_ordered_TA_skills_after_BVS<- extract_log_lik(ordered_TA_skills_after_BVS)
r_eff_log_lik_ordered_TA_skills_after_BVS<- relative_eff(exp(log_lik_ordered_TA_skills_after_BVS),
                                                         chain_id=rep(1:2,each=10000))

# ---WAIC, DIC, LOOIC Bayesian Information Criteria for the ordered logistic with only skills as covariates

waic(log_lik_ordered_TA_skills_after_BVS)####256.9
loo(log_lik_ordered_TA_skills_after_BVS,r_eff=r_eff_log_lik_ordered_TA_skills_after_BVS)#258.1for model with proper thinning 379,9
DIC_Gelman(deviance_ordered_TA_skills_after_BVS)#249.4


###--------------Posterior Summary Statistics-Analysis------------------------########
##----Parameters Names

## Vector of teams names along with
## their ranking positions, points, abilities
teams <- c("Ethnikos Alexandroupolis", "Foinikas Syrou", "Iraklis Chalkidas",       
           "Iraklis Petosfairishs" ,"Kifisia", "Kyzikos Peramou" ,        
           "Olympiacos" ,"Orestiada" ,"Pamvochaikos" ,           
           "Panachaiki",  "Panathinaikos","Paok")        
# teams  levels(datafr_teams_scores_set_skills$away_Team)
observed_positions<-c("(7)","(6)","(9)","(8)","(5)","(11)","(1)","(12)","(4)","(10)","(3)","(2)")
observed_points<-c("(36)","(37)","(16)","(28)","(38)","(14)","(62)","(7)","(39)","(16)","(50)","(53)")


teams_gen_abil<-paste0(teams," ","Ability")


teams_pos<-paste0(teams," ",observed_positions)
teams_points<-paste0(teams," ",observed_points)

skill_events_differences <- c("failed serve","failed pass",
                              "perfect att1", "failed att1",
                              "perfect att2", "failed att2",
                              "perfect block","failed setting")

cutpoints<-c("c_1","c_2","c_3","c_4","c_5","c_6")

#####----------------------Posterior summary----------------------------######

names(ordered_TA_skills_after_BVS)[1:8]<-skill_events_differences
names(ordered_TA_skills_after_BVS)[c(20,25:29)]<-cutpoints
names(ordered_TA_skills_after_BVS)[c(30:41)]<-teams
names(ordered_TA_skills_after_BVS)[c(9:19)]<-teams[c(1:11)]



print(ordered_TA_skills_after_BVS,
      pars=c(
        "first_temp_Intercept",
        "temp_Intercept","beta","gen_abil"),probs = c(0.025,0.5,0.975), digits=2)


## Extraction of model parameters
sims <- rstan::extract(ordered_TA_skills_after_BVS)

teams_index <- teams_pos


beta <- sims$beta
first_temp_Intercept<- sims$first_temp_Intercept
temp_Intercept<- sims$temp_Intercept
temp_intercepts<-cbind(first_temp_Intercept,temp_Intercept)
gen_abil <- sims$gen_abil


## Order of ability parameters (based on the posterior means)
beta_hat <- apply(beta,2, median)
first_temp_Intercept_hat <- apply(first_temp_Intercept,1,median)
temp_intercepts_hat <- apply(temp_intercepts,2,median)
gen_abil_hat <- apply(gen_abil,2, median)

beta_hat_ord <- order(beta_hat, decreasing = TRUE)
first_temp_Intercept_hat_ord <- order(first_temp_Intercept_hat, decreasing = TRUE)
temp_intercepts_hat_ord <- order(temp_intercepts_hat, decreasing = TRUE)
gen_abil_hat_ord <- order(gen_abil_hat, decreasing = TRUE)

##---Proper parameters renaming

colnames(beta)<-skill_events_differences
colnames(temp_intercepts)<-cutpoints
colnames(gen_abil)<-teams



## Data frame of parameters in terms of convenience in both tables and graphs

beta<-as.data.frame(beta)
temp_intercepts<-as.data.frame(temp_intercepts)
gen_abil<-as.data.frame(gen_abil)

###-----MCMC Posterior Intervals


color_scheme_set("brightblue")


pdf(file="Ordered_TA_Skills_cutpoints.pdf", width =12, height =7.5)

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

pdf(file="Ordered_TA_Skills_Skills_Differences.pdf", width =12, height =7.5)

mcmc_intervals(beta[,c(beta_hat_ord)],
                      prob = 0.95,prob_outer=0.95,
                      point_est = c( "mean"))+ggtitle("Skill Differences")+xlim(-1,1)+
  scale_x_continuous(breaks = seq(from = -1, to = 1, by = 0.2))+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))

dev.off()

pdf(file="Ordered_TA_Skills_General_Abilities.pdf", width =12, height =7.5)

mcmc_intervals(gen_abil[,c(gen_abil_hat_ord)],
                      prob = 0.95,prob_outer=0.95,
                      point_est = c( "mean"))+ggtitle("Net General Abilities")+xlim(-2.5,2.5)+
  scale_x_continuous(breaks = seq(from = -2.5, to = 2.5, by = 0.5))+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))

dev.off()

pdf(file="Ordered_TA_Skills_General_Abilities(Areas).pdf", width =12, height =7.5)

mcmc_areas(gen_abil[,c(gen_abil_hat_ord)],
                      prob = 0.95,prob_outer=0.95,
                      point_est = c( "mean"))+ggtitle("Net General Abilities")+xlim(-2.5,2.5)+
  scale_x_continuous(breaks = seq(from = -2.5, to = 2.5, by = 0.5))+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))
dev.off()



