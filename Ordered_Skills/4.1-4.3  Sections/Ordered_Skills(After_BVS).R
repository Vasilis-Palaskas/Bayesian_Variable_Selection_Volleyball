# Load the proper libraries.
library(rstan)
library(coda)
library(shinystan)
# Choose the working directory of this file (.../BVS_Paper/Ordered_TA_Skills)
setwd()
# Load the properly full prepared data ("datalist_ordered") for the ordered logistic models.
load("datalist_ordered")


#Rename the proper names

names(dataList$X)<-c("perfect serve","very good serve","failed serve"," perfect pass","
                                  very good pass","poor pass","failed pass","perfect att1","blocked att1",
                     "failed att1","perfect att2","blocked att2","failed att2","perfect block",
                     "block net violation","failed block","failed setting")

# Keep these skills variables with post.incl.prob.>0.5
# 
X_ordered_TA_Skills<-dataList$X[,colnames(dataList$X)%in%
                                  c("failed serve","failed pass",
                                    "perfect att1", "failed att1",
                                    "perfect att2", "failed att2",
                                    "perfect block","failed setting") ]
dataList<-list(Y=dataList$Y,X=X_ordered_TA_Skills,n_teams=12,
               N=dataList$N,K=ncol(X_ordered_TA_Skills),ncat=6)


## Run 4.1 BVS for the Ordered multinomial/Ordered_TA_Skills(After BVS).stan
ordered_skills_after_BVS<-stan(file.choose(),iter=12000, warmup=2000,chains=1,thin=2,
                                  data=dataList,control=list(max_treedepth=15),cores=1)
save(ordered_skills_after_BVS,file="ordered_skills_after_BVS")
