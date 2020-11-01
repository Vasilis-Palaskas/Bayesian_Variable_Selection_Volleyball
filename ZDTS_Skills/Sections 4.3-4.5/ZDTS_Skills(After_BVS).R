
# Load the proper libraries.
library(rstan)
library(coda)
library(shinystan)
library(loo)
# Choose the working directory of this file (...\\Submitted_Appendix\\ZDTS\\)
setwd("C:/Users/vasileios palaskas/Desktop/BVS_Paper/ZDTS_Skills")
# Load the properly prepared data for both home and away skill events as well as
# both home and away teams in each match
load("X_home")
load("X_away")
load("data_zdts_skills")
# 
# load("X_home")
# load("X_away")
# load("data_zdts_skills")
#Model matrices for home and away sets scored, respectively


#### Standardization of the model Matrices for numerical convenience

#Rename the columns
colnames(X_home)<-c("home_perfect_serve","home_very_good_serve","home_failed_serve",
                    "home_perfect_pass","home_very_good_pass","home_poor_pass","home_failed_pass",
                    "home_perfect_att1","home_blocked_att1","home_failed_att1",
                    "home_perfect_att2","home_blocked_att2","home_failed_att2",
                    "home_perfect_block","home_net_violation_block","home_failed_block","home_failed_setting")

colnames(X_away)<-c("away_perfect_serve","away_very_good_serve","away_failed_serve",
                    "away_perfect_pass","away_very_good_pass","away_poor_pass","away_failed_pass",
                    "away_perfect_att1","away_blocked_att1","away_failed_att1",
                    "away_perfect_att2","away_blocked_att2","away_failed_att2",
                    "away_perfect_block","away_net_violation_block","away_failed_block","away_failed_setting")

#### Standardization of the Model Matrices for numerical convenience
X_home_std<-X_away_std<-matrix(NA,nrow=132,ncol=17)
for (i in 1:dim(X_home)[2]){
  X_home_std[,i]<-(X_home[,i]-mean(X_home[,i]))/sd(X_home[,i])
  X_away_std[,i]<-(X_away[,i]-mean(X_away[,i]))/sd(X_away[,i])
}

colnames(X_home_std)<-c("home_perfect_serve","home_very_good_serve","home_failed_serve",
                        "home_perfect_pass","home_very_good_pass","home_poor_pass","home_failed_pass",
                        "home_perfect_att1","home_blocked_att1","home_failed_att1",
                        "home_perfect_att2","home_blocked_att2","home_failed_att2",
                        "home_perfect_block","home_net_violation_block","home_failed_block","home_failed_setting")

colnames(X_away_std)<-c("away_perfect_serve","away_very_good_serve","away_failed_serve",
                        "away_perfect_pass","away_very_good_pass","away_poor_pass","away_failed_pass",
                        "away_perfect_att1","away_blocked_att1","away_failed_att1",
                        "away_perfect_att2","away_blocked_att2","away_failed_att2",
                        "away_perfect_block","away_net_violation_block","away_failed_block","away_failed_setting")

# Κeep these skill events with posterior inclusion probabilities >0.5
X_home_std<-X_home_std[,colnames(X_home_std)%in%c( "home_perfect_pass",
                                                  "home_very_good_pass","home_poor_pass","home_failed_pass",
                                                  "home_blocked_att1"
                                                 )
                       ]
X_away_std<-X_away_std[,colnames(X_away_std)%in%c("away_failed_serve","away_failed_pass",
                                                   "away_blocked_att1","away_failed_att1",
                                                  "away_perfect_att2","away_failed_att2","away_net_violation_block","away_failed_block")
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
ZDTS_Skills_after_BVS<-stan(file.choose(),
                            data=data_zdts_only_skills,chains=1,init_r=0.5,
                            iter=12000,warmup=2000)### R

save(ZDTS_Skills_after_BVS,file="ZDTS_Skills_after_BVS")
