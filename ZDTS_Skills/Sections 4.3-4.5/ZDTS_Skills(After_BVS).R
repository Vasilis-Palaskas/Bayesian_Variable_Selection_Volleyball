
# Load the proper libraries.
library(rstan)
library(coda)
library(shinystan)
library(loo)
library(coefplot)
library(bayesplot)
# Choose the working directory of this file (.../Bayesian_Variable_Selection_Volleyball/ZDTS_Skills)
 setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/ZDTS_Skills")
# Load the properly prepared data for both home and away skill events as well as
# both home and away teams in each match
load("X_home")
load("X_away")
load("data_zdts_skills")
# # Choose the working directory of this file (.../Bayesian_Variable_Selection_Volleyball/ZDTS_Skills/Sections 4.3-4.5")

setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/ZDTS_Skills/Sections 4.3-4.5")
#### Standardization of the model Matrices for numerical convenience

#Rename the columns
colnames(X_home)<-c("(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass","
                                 (Home) very good pass","(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
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

colnames(X_home_std)<-c("(Home) perfect serve","(Home) very good serve","(Home) failed serve","(Home) perfect pass","(Home) very good pass",
                        "(Home) poor pass","(Home) failed pass","(Home) perfect att1","(Home) blocked att1",
                        "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2",
                            "(Home) perfect block",
                        "(Home) block net violation","(Home) failed block",
                        "(Home) failed setting") 

colnames(X_away_std)<-c("(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass",
        "(Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
                        "(Away) failed att1","(Away) perfect att2","(Away) blocked att2",
                      "(Away) failed att2","(Away) perfect block",
                        "(Away) block net violation","(Away) failed block","(Away) failed setting")

# Κeep these skill events with posterior inclusion probabilities >0.5
X_home_std<-X_home_std[,colnames(X_home_std)%in%c( "(Home) perfect pass",
                                                  "(Home) very good pass","(Home) poor pass",
                                                  "(Home) failed pass",
                                                  "(Home) blocked att1"
                                                 )
                       ]
X_away_std<-X_away_std[,colnames(X_away_std)%in%c("(Away) failed serve","(Away) failed pass",
                                                   "(Away) blocked att1","(Away) failed att1",
                                                  "(Away) perfect att2","(Away) failed att2",
                                                  "(Away) block net violation",
                                                  "(Away) failed block")
]

data_zdts_only_skills<-list(n_games=data_zdts_skills$N,
                            n_teams=data_zdts_skills$n_teams,
                            X_home=X_home_std,X_away=X_away_std,K_home=ncol(X_home_std),K_away=ncol(X_away_std),
                            home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)


## Run ZDTS_Skills_after_BVS.stan
ZDTS_Skills_after_BVS<-stan("ZDTS_Skills_after_BVS.stan",
                            data=data_zdts_only_skills,chains=2,init_r=0.5,
                            iter=12000,warmup=2000)### R

save(ZDTS_Skills_after_BVS,file="ZDTS_Skills_after_BVS")

#### Model Predictive Performance
##Table 3 :Information Criteria
# Calculation of the DIC (Gelman,2004)
DIC_Gelman<-function(dev){
  res<-mean(dev)+0.5*var(dev)
  return(res)
}




####Extraction of the log-likelihood, deviance quantities

deviance_ZDTS_Skills_after_BVS<-extract(ZDTS_Skills_after_BVS,pars="dev")$dev
log_lik_ZDTS_Skills_after_BVS<- extract_log_lik(ZDTS_Skills_after_BVS)
r_eff_log_lik_ZDTS_Skills_after_BVS<- relative_eff(exp(log_lik_ZDTS_Skills_after_BVS),chain_id=rep(1:2,each=12000))

##WAIC, LOO, DIC

waic(log_lik_ZDTS_Skills_after_BVS)####235.8
loo(log_lik_ZDTS_Skills_after_BVS)#236.3
loo(log_lik_ZDTS_Skills_after_BVS,r_eff=r_eff_log_lik_ordered_skills_after_BVS_model1)#for model with proper thinning 379,9
DIC_Gelman(deviance_ZDTS_Skills_after_BVS)#238.0


###------------ MCMC Posterior Summary Plots

sims <- rstan::extract(ZDTS_Skills_after_BVS)

# mu, home
 mu <- sims$mu
 home <- sims$home
 
 mu_hat <- apply(mu,1, median)
 mu_sd <- apply(mu,1, sd)
 
 home_hat <- apply(home,1, median)
 home_sd <- apply(home,1, sd)
## coefplot for skill events
Home_skill_events<-c("(Home) perfect pass",
                      "(Home) very good pass","(Home) poor pass","(Home) failed pass",
                                 "(Home) blocked att1"
)

#home team's skill events
beta_home <- sims$beta_home
beta_home_hat <- apply(beta_home,2, median)
beta_home_sd <- apply(beta_home,2, sd)
ord <- order(beta_home_hat, decreasing = TRUE)


coefplot( rev(beta_home_hat[ord]), 
          rev(beta_home_sd[ord]), 
          CI=2, 
          varnames=rev(Home_skill_events[ord]), 
          main="Home_skill_events (estim. +/- 2 s.e.)\n", 
          cex.var=1.5, mar=c(1,6,4.5,1),
          cex.main=1.3,pch=16, cex=2, col="blue")

#away team's skill events
Away_skill_events<-c("(Away) failed serve","(Away) failed pass",
                     "(Away) blocked att1","(Away) failed att1",
                     "(Away) perfect att2","(Away) failed att2","(Away) block net violation",
                     "(Away) failed block")


#away team's skill events
beta_away <- sims$beta_away
beta_away_hat <- apply(beta_away,2, median)
beta_away_sd <- apply(beta_away,2, sd)

ord <- order(beta_away_hat, decreasing = TRUE)


coefplot( rev(beta_away_hat[ord]), 
          rev(beta_away_sd[ord]), 
          CI=2, 
          varnames=rev(as.character(Away_skill_events)[ord]), 
          main="Away skill events (estim. +/- 2 s.e.)\n", 
          cex.var=1.5, mar=c(1,6,4.5,1),
          cex.main=1.3,pch=16, cex=2, col="blue")


### MCMC Areas for both home and away teams' skill events
#Posterior means for ability parameters
team_abil_final_ordered_logistic_means<-apply(team_abil_final_ordered_logistic,2,mean) ##

# Order of ability parameters (based on the posterior means)
team_abil_order_final_ordered_logistic<-order(team_abil_final_ordered_logistic_means,decreasing=T)


##Posterior 95% uncertainty intervals

## Renaming model parameters in terms of convenience in both tables and graphs
## Teams names

names(final_ordered_logistic)[1:5]<-c("c1","c2","c3","c4","c5")
names(final_ordered_logistic)[17:28]<-c(teams_names)[1:12]

colnames(team_abil_final_ordered_logistic)[1:12]<-c(teams_names)[1:12]
colnames(c_final_ordered_logistic)<-c("c1","c2","c3","c4","c5")

## Teams names+points
names(final_ordered_logistic)[17:28]<-c(teams_names_points)[1:12]

colnames(team_abil_final_ordered_logistic)<-c(teams_names_points)[1:12]
## Teams names+positions
names(final_ordered_logistic)[17:28]<-c(teams_names_pos)[1:12]
colnames(team_abil_final_ordered_logistic)<-c(teams_names_pos)[1:12]

## Data frame of parameters (merged chains) in terms of convenience in both tables and graphs

posterior_c<-as.data.frame(c_final_ordered_logistic)
posterior_team_abil<-as.data.frame(team_abil_final_ordered_logistic)

## Convertion to array (necessary for the summaries)
array_posterior_final_ordered_logistic<-as.array(final_ordered_logistic)


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

