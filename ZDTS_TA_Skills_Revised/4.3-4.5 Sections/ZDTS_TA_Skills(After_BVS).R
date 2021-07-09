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
setwd("C:/Users/Bill_/Desktop/Github Projects/Bayesian_Variable_Selection_Volleyball/ZDTS_TA_Skills")
# both home and away teams in each match
load("X_home")
load("X_away")
load("data_zdts_skills")


## Vector of teams names along with
## their ranking positions, points, abilities
teams <- teams <- c("Ethnikos Alexandroupolis", "Foinikas Syrou", "Iraklis Chalkidas",       
                    "Iraklis Petosfairishs" ,"Kifisia", "Kyzikos Peramou" ,        
                    "Olympiacos" ,"Orestiada" ,"Pamvochaikos" ,           
                    "Panachaiki",  "Panathinaikos","Paok") 
observed_positions<-c("(7)","(6)","(9)","(8)","(5)","(11)","(1)","(12)","(4)","(10)","(3)","(2)")
observed_points<-c("(36)","(37)","(16)","(28)","(38)","(14)","(62)","(7)","(39)","(16)","(50)","(53)")


  teams_attack<-paste0(teams," ","Attack")
  teams_defense<-paste0(teams," ","Defense")
  teams_over<-paste0(teams," ","Overall")
  
  teams_pos<-paste0(teams," ",observed_positions)
  teams_points<-paste0(teams," ",observed_points)



#----Rename properly the skill variables


##----Skill events selected via the BVS process based on PSI Median Threshold
#### Standardization of the Model Matrices for numerical convenience
X_home_std<-X_away_std<-matrix(NA,nrow=132,ncol=17)
for (i in 1:dim(X_home)[2]){
  X_home_std[,i]<-(X_home[,i]-mean(X_home[,i]))/sd(X_home[,i])
  X_away_std[,i]<-(X_away[,i]-mean(X_away[,i]))/sd(X_away[,i])
}

#Rename the columns
colnames(X_home_std)<-c("(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass",
                    "(Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
                    "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
                    "(Home) block net violation","(Home) failed block","(Home) failed setting") 

colnames(X_away_std)<-c("(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass",
                    "(Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
                    "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
                    "(Away) block net violation","(Away) failed block","(Away) failed setting")
#--Selected covariates via the Gibbs BVS Method

final_X_home_std<-X_home_std[,colnames(X_home_std)%in%c("(Home) poor pass",
                                                        "(Home) failed pass",
                                                        "(Home) failed block")]

final_X_away_std<-X_away_std[,colnames(X_away_std)%in%c("(Away) failed serve",
                                                  "(Away) poor pass",
                                              "(Away) failed pass",
                                              "(Away) blocked att1",
                                              "(Away) failed att1",
                                             "(Away) block net violation",
                                             "(Away) failed block")]
#---Datalist required for the ZDTS with TA+Skills
data_zdts_ta_skills<-list(c_thres=2,c_std=5,
                          n_games=data_zdts_skills$N,
                          away_team=as.numeric(data_zdts_skills$away_team),
                          home_team=as.numeric(data_zdts_skills$home_team),
                          n_teams=data_zdts_skills$n_teams,
                          X_home=final_X_home_std,X_away=final_X_away_std,
                          K_home=ncol(final_X_home_std),
                          K_away=ncol(final_X_away_std),
                          home_sets=data_zdts_skills$home_sets,
                          away_sets=data_zdts_skills$away_sets)

## Run ZDTS_TA_Skills_after_BVS.stan
ZDTS_TA_Skills_after_BVS<-stan("ZDTS_TA_Skills_after_BVS.stan",
                               data=data_zdts_ta_skills,chains=2,
                               init_r=0.5,
                               iter=12000,warmup=2000)### R

save(ZDTS_TA_Skills_after_BVS,file="ZDTS_TA_Skills_after_BVS")

###--------------Predictive Model Performance Evaluation-------------------------########
# Calculation of the DIC (Gelman,2004)
DIC_Gelman<-function(dev){
  res<-mean(dev)+0.5*var(dev)
  return(res)
}
###----Extraction of several posterior quantities (deviances, log-likelihoods)

deviance_ZDTS_TA_Skills_after_BVS<-rstan::extract(ZDTS_TA_Skills_after_BVS,pars="dev")$dev
log_lik_ZDTS_TA_Skills_after_BVS<- extract_log_lik(ZDTS_TA_Skills_after_BVS)
r_eff_log_lik_ZDTS_TA_Skills_after_BVS<- relative_eff(exp(log_lik_ZDTS_TA_Skills_after_BVS),chain_id=rep(1:2,each=10000))

# ---WAIC, DIC, LOOIC Bayesian Information Criteria for the ordered logistic with only skills as covariates

waic(log_lik_ZDTS_TA_Skills_after_BVS)####271.1 (20.7)
loo(log_lik_ZDTS_TA_Skills_after_BVS,
    r_eff=r_eff_log_lik_ZDTS_TA_Skills_after_BVS)# 274.8 (21.1)
DIC_Gelman(deviance_ZDTS_TA_Skills_after_BVS)#269.6

###--------------Posterior Summary Statistics-Analysis------------------------########

#####----------------------Posterior summary----------------------------######

names(ZDTS_TA_Skills_after_BVS)[35:46]<-c(teams_attack)[1:12]
names(ZDTS_TA_Skills_after_BVS)[47:58]<-c(teams_defense)[1:12]
names(ZDTS_TA_Skills_after_BVS)[719:730]<-c(teams_over)[1:12]
names(ZDTS_TA_Skills_after_BVS)[1:3]<-colnames(final_X_home_std)
names(ZDTS_TA_Skills_after_BVS)[4:10]<-colnames(final_X_away_std)


print(ZDTS_TA_Skills_after_BVS,
      pars=c("mu","home","attack","defense","overall",
             "beta_home","beta_away","dev"),probs = c(0.025,0.5,0.975), digits=2)
