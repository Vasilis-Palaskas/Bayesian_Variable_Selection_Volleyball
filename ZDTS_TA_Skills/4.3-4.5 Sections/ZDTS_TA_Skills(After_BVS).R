
# Load the proper libraries.
library(rstan)
library(coda)
library(shinystan)
library(loo)
# Choose the working directory of this file (.../Bayesian_Variable_Selection_Volleyball/ZDTS_Skills)
setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/ZDTS_TA_Skills")# Load the properly prepared data for both home and away skill events as well as
# both home and away teams in each match
load("X_home")
load("X_away")
load("data_zdts_skills")
# # Choose the working directory of this file (.../Bayesian_Variable_Selection_Volleyball/ZDTS_Skills/Sections 4.3-4.5")

setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/ZDTS_TA_Skills/4.3-4.5 Sections")


#Rename the columns
colnames(X_home)<-c("(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass",
                    " (Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
                    "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
                    "(Home) block net violation","(Home) failed block","(Home) failed setting") 
  
colnames(X_away)<-c("(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass",
                    "(Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
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

X_home_std<-X_home_std[,colnames(X_home_std)%in%c( "(Home) poor pass","(Home) failed pass",
                                                   "(Home) blocked att1"
)
]
X_away_std<-X_away_std[,colnames(X_away_std)%in%c("(Away) failed serve","(Away) poor pass","(Away) failed pass",
                                                  "(Away) blocked att1","(Away) failed att1",
                                                  "(Away) block net violation","(Away) failed block")
]
data_zdts_skills<-list(n_games=data_zdts_skills$N,
                       away_team=as.numeric(data_zdts_skills$away_team),
                       home_team=as.numeric(data_zdts_skills$home_team),
                       n_teams=data_zdts_skills$n_teams,
                       X_home=X_home_std,X_away=X_away_std,K_home=ncol(X_home_std),K_away=ncol(X_away_std),
                       home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)


## Run ZDTS_TA_Skills_after_BVS.stan
ZDTS_TA_Skills_after_BVS<-stan("ZDTS_TA_Skills_after_BVS.stan",
                       data=data_zdts_skills,chains=2,init_r=0.5,
                       iter=12000,warmup=2000)### R

save(ZDTS_TA_Skills_after_BVS,file="ZDTS_TA_Skills_after_BVS")


#### Model Predictive Performance
##Table 3 :Information Criteria
# Calculation of the DIC (Gelman,2004)
DIC_Gelman<-function(dev){
  res<-mean(dev)+0.5*var(dev)
  return(res)
}





####Extraction of the log-likelihood, deviance quantities

deviance_ZDTS_TA_Skills_after_BVS<-extract(ZDTS_TA_Skills_after_BVS,pars="dev")$dev
log_lik_ZDTS_TA_Skills_after_BVS<- extract_log_lik(ZDTS_TA_Skills_after_BVS)
r_eff_log_lik_ZDTS_TA_Skills_after_BVS<- relative_eff(exp(log_lik_ZDTS_TA_Skills_after_BVS),chain_id=rep(1:2,each=10000))

##WAIC, LOO, DIC

waic(log_lik_ZDTS_TA_Skills_after_BVS)####267.4
loo(log_lik_ZDTS_TA_Skills_after_BVS)#270.5
loo(log_lik_ZDTS_TA_Skills_after_BVS,r_eff=r_eff_log_lik_ZDTS_TA_Skills_after_BVS)#for model with proper thinning 379,9
DIC_Gelman(deviance_ZDTS_TA_Skills_after_BVS)#265.6






### MCMC Posterior Summary Plots


sims <- rstan::extract(ordered_TA_skills_after_BVS)


# mu, home
mu <- sims$mu
home <- sims$home

mu_hat <- apply(home,2, median)
mu_sd <- apply(home,2, sd)

home_hat <- apply(home,2, median)
home_sd <- apply(home,2, sd)
## coefplot for team abilities
teams <- c("Ethnikos Alexandroupolis", "Pamvochaikos",
           "Iraklis Petosfairishs",   "Kyzikos Peramou",   
           "Panachaiki",    "Foinikas Syrou",          
           "Kifisia",  "Orestiada",  "Olympiacos",              
           "Panathinaikos",  "Iraklis Chalkidas",  "Paok") 
teams_index <- unique(dataList$home_team)

gen_abil <- sims$overall
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
Home_skill_events<-   c( "(Home) perfect pass",
                         "(Home) very good pass","(Home) poor pass","(Home) failed pass",
                         "(Home) blocked att1"
)

#home team's skill events
beta_home <- sims$beta_home
beta_home_hat <- apply(beta_home,2, median)
beta_home_sd <- apply(beta_home,2, sd)
#teams_index <- match(squadre16_17, teams)
# att <- att[,3,teams_index]
# def <- def[,3,teams_index]
# att_hat <- apply(att,2,median)
# att_sd <- apply(att,2,sd)
# def_hat <- apply(def,2,median)
# def_sd <- apply(def,2,sd)
ord <- order(beta_home_hat, decreasing = TRUE)
# ord_2 <- order(def_hat)

coefplot( rev(beta_home_hat[ord]), 
          rev(beta_home_sd[ord]), 
          CI=2, 
          varnames=rev(as.character(skill_events_differences)[ord]), 
          main="Home_skill_events (estim. +/- 2 s.e.)\n", 
          cex.var=1.5, mar=c(1,6,4.5,1),
          cex.main=1.3,pch=16, cex=2, col="blue")
#away team's skill events
Away_skill_events<-c("(Away) failed serve","(Away) failed pass",
                     "(Away) blocked att1","(Away) failed att1",
                     "(Away) perfect att2","(Away) failed att2","(Away) block net violation",
                     "(Away) failed block")


#home team's skill events
beta_away <- sims$beta_away
beta_away_hat <- apply(beta_away,2, median)
beta_away_sd <- apply(beta_away,2, sd)
#teams_index <- match(squadre16_17, teams)
# att <- att[,3,teams_index]
# def <- def[,3,teams_index]
# att_hat <- apply(att,2,median)
# att_sd <- apply(att,2,sd)
# def_hat <- apply(def,2,median)
# def_sd <- apply(def,2,sd)
ord <- order(beta_away_hat, decreasing = TRUE)
# ord_2 <- order(def_hat)

coefplot( rev(beta_away_hat[ord]), 
          rev(beta_away_sd[ord]), 
          CI=2, 
          varnames=rev(as.character(Away_skill_events)[ord]), 
          main="Away skill events (estim. +/- 2 s.e.)\n", 
          cex.var=1.5, mar=c(1,6,4.5,1),
          cex.main=1.3,pch=16, cex=2, col="blue")






### MCMC Areas for both home and away teams' skill events

plot_home_skill_events<-	mcmc_areas(beta_home,
                                    prob = 0.95,point_est = c( "mean"))+ggtitle("Home Teams' Skill Events")+
  scale_x_continuous(lim=c(0,10))+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 23),axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),axis.line.y=element_blank())

plot_away_skill_events<-mcmc_areas(beta_away,
                                   prob = 0.95,point_est = c( "mean"))+ggtitle("Away Teams' Skill Events")+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 23))

grid.arrange(plot_home_skill_events,plot_away_skill_events,ncol=2)

ggarrange(plot_home_skill_events,plot_away_skill_events,nrow=2,ncol=2)