# Load the proper libraries.
library(rstan)
library(coda)

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
X_home<-data_zdts_skills[,c(1:17)]
X_away<-data_zdts_skills[,c(18:34)]

#### Standardization of the model Matrices for numerical convenience
X_home_std<-X_away_std<-matrix(NA,nrow=132,ncol=17)
for (i in 1:dim(X_home)[2]){
  X_home_std[,i]<-(X_home[,i]-mean(X_home[,i]))/sd(X_home[,i])
  X_away_std[,i]<-(X_away[,i]-mean(X_away[,i]))/sd(X_away[,i])
}

# Κeep these skill events with posterior inclusion probabilities >0.5


data_zdts_only_skills<-list(n_games=data_zdts_skills$N,
                            n_teams=data_zdts_skills$n_teams,
                            X_home=X_home_std,X_away=X_away_std,K=ncol(X_home_std),
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