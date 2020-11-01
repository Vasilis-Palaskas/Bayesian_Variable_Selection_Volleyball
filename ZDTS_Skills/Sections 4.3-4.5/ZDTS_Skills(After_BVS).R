
# Load the proper libraries.
library(rstan)
library(coda)
library(shinystan)
library(loo)
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

colnames(X_away)<-c("(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass","
                                 (Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
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
                        "(Home) failed att1","(Home) perfect att2","(Home) blocked att2","(Home) failed att2","(Home) perfect block",
                        "(Home) block net violation","(Home) failed block","(Home) failed setting") 

colnames(X_away_std)<-c("(Away) perfect serve","(Away) very good serve","(Away) failed serve","(Away) perfect pass","
                                 (Away) very good pass","(Away) poor pass","(Away) failed pass","(Away) perfect att1","(Away) blocked att1",
                        "(Away) failed att1","(Away) perfect att2","(Away) blocked att2","(Away) failed att2","(Away) perfect block",
                        "(Away) block net violation","(Away) failed block","(Away) failed setting")

# Κeep these skill events with posterior inclusion probabilities >0.5
X_home_std<-X_home_std[,colnames(X_home_std)%in%c( "(Home) perfect pass",
                                                  "(Home) very good pass","(Home) poor pass","(Home) failed pass",
                                                  "(Home) blocked att1"
                                                 )
                       ]
X_away_std<-X_away_std[,colnames(X_away_std)%in%c("(Away) failed serve","(Away) failed pass",
                                                   "(Away) blocked att1","(Away) failed att1",
                                                  "(Away) perfect att2","(Away) failed att2","(Away) block net violation",
                                                  "(Away) failed block")
]
data_zdts_only_skills<-list(n_games=data_zdts_skills$N,
                            n_teams=data_zdts_skills$n_teams,
                            X_home=X_home_std,X_away=X_away_std,K_home=ncol(X_home_std),K_away=ncol(X_away_std),
                            home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)

# data_zdts_skills<-list(n_games=data_zdts_skills$N,
#                        n_teams=data_zdts_skills$n_teams,
#                        X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
#                        home_sets=data_zdts_skills$home_sets,away_sets=data_zdts_skills$away_sets)

## Run ZDTS_Skills_after_BVS.stan
ZDTS_Skills_after_BVS<-stan("ZDTS_Skills_after_BVS.stan",
                            data=data_zdts_only_skills,chains=2,init_r=0.5,
                            iter=12000,warmup=2000)### R

save(ZDTS_Skills_after_BVS,file="ZDTS_Skills_after_BVS")

### MCMC Posterior Summary Plots

sims <- rstan::extract(ZDTS_Skills_after_BVS)

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