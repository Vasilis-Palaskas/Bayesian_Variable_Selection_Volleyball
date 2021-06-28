# Load the proper libraries.
library(rstan)
library(coda)
library(bayesplot)
library(ggmcmc)
# Choose the working directory of this file (.../BVS_Paper/Ordered_TA_Skills)
setwd("C:\\Users\\vasileios palaskas\\Desktop\\Github folder\\Bayesian_Variable_Selection_Volleyball\\Ordered_TA_Skills")
# Load the properly full prepared data ("datalist_ordered") for the ordered logistic models.
#------Skills for both Home and Away Teams
X_home<-data_by_sets[c(
  "Home_perfect_serves","Home_very_good_serves",
  "Home_failed_serves","Home_perfect_passes","Home_very_good_passes",
  "Home_poor_passes","Home_failed_passes","Home_perfect_att1",
  "Home_blocked_att1","Home_failed_att1","Home_perfect_att2",
  "Home_blocked_att2","Home_failed_att2","Home_perfect_blocks",
  "Home_net_violation_blocks","Home_failed_blocks","Home_failed_settings")
]

X_away<-data_by_sets[c(
  "Away_perfect_serves","Away_very_good_serves",
  "Away_failed_serves","Away_perfect_passes","Away_very_good_passes",
  "Away_poor_passes","Away_failed_passes","Away_perfect_att1",
  "Away_blocked_att1","Away_failed_att1","Away_perfect_att2",
  "Away_blocked_att2","Away_failed_att2","Away_perfect_blocks",
  "Away_net_violation_blocks","Away_failed_blocks","Away_failed_settings")
]



#----Rename properly the skill variables
##----Skill events selected via the BVS process based on PSI Median Threshold
#### Standardization of the Model Matrices for numerical convenience

# Load the properly prepared data ("Data_ordered_skills").
# load("datalist_ordered")

X_home_diff<-data.frame(X_home-X_away)
colnames(X_home_diff)<-c(
  "perfect_serves","very_good_serves",
  "failed_serves","perfect_passes","very_good_passes",
  "poor_passes","failed_passes","perfect_att1",
  "blocked_att1","failed_att1","perfect_att2",
  "blocked_att2","failed_att2","perfect_blocks",
  "net_violation_blocks","failed_blocks","failed_settings")

#---Transform set difference values in terms of fitting for ordered multinomial model (requires positive integers or factors)
data_by_sets$sets_difference_factor<-data_by_sets$sets_difference
for (i in 1:dim(data_by_sets)[1]){
  if (data_by_sets$sets_difference[i]==(-3)){
    data_by_sets$sets_difference_factor[i]<-1
  } else if (data_by_sets$sets_difference[i]==(-2)){
    data_by_sets$sets_difference_factor[i]<-2
  } else if (data_by_sets$sets_difference[i]==(-1)){
    data_by_sets$sets_difference_factor[i]<-3
  } else if (data_by_sets$sets_difference[i]==(1)){
    data_by_sets$sets_difference_factor[i]<-4
  } else if (data_by_sets$sets_difference[i]==(2)){
    data_by_sets$sets_difference_factor[i]<-5
  } else if (data_by_sets$sets_difference[i]==(3)){
    data_by_sets$sets_difference_factor[i]<-6
  }
  
}
## Vector of teams names along with
## their ranking positions, points, abilities
teams <- levels(data_by_sets$home_Team)
observed_positions<-c("(7)","(6)","(9)","(8)","(5)","(11)","(1)","(12)","(4)","(10)","(3)","(2)")
observed_points<-c("(36)","(37)","(16)","(28)","(38)","(14)","(62)","(7)","(39)","(16)","(50)","(53)")


teams_attack<-paste0(teams," ","Attack")
teams_defense<-paste0(teams," ","Defense")
teams_over<-paste0(teams," ","Overall")

teams_pos<-paste0(teams," ",observed_positions)
teams_points<-paste0(teams," ",observed_points)
#Numerize the factors in terms of your convenience
dataList<-list(Y=data_by_sets$sets_difference_factor,X=X_home_diff,n_teams=length(levels(data_by_sets$home_Team)),
               N=dim(data_by_sets)[1],K=ncol(X_home_diff),ncat=6,
	        home_team=as.numeric(dataList$home_team),
      away_team=as.numeric(dataList$away_team))



## Run Full_ordered_skills.stan
Full_ordered_team_abilities_skills<-stan("Full_ordered_team_abilities_skills.stan",iter=10000, warmup=2000,
                                         chains=2,thin=2,
                                         data=dataList,control=list(max_treedepth=15),cores=2)

#save(Full_ordered_team_abilities_skills,file="Full_ordered_team_abilities_skills")

# Load the results from the full ordered logistic model (including both team abilities and all candidate variables).
load(file="Full_ordered_team_abilities_skills")

# Extract the posterior summary statistics of both candidate variables' parameters and the rest of other parameters.

betas_summary<-summary(Full_ordered_team_abilities_skills,pars = c("beta"))$summary
gen_abil_raw_summary<-summary(Full_ordered_team_abilities_skills,pars = c("gen_abil_raw"))$summary
intercept_summary<-summary(Full_ordered_team_abilities_skills, pars = c("temp_Intercept"))$summary

# Use their posterior means and standard deviations for both initial values specification and prior specification.

post_mean_betas<-betas_summary[1:dataList$K,1]
post_sd_betas<-betas_summary[1:dataList$K,3]

# Step 1: Initialization of the model parameters.

gammas<-rep(1,dataList$K)# All the candidate variables included in the model
temp_Intercept<-intercept_summary[,1]
betas<-post_mean_betas
gen_abil_raw<-gen_abil_raw_summary[,1]

# Prepare the vectors with the posterior samples of dimension Txp (p=K during algorithm iterations) for all gammas and betas coefficients , respectively.
gammas_matrix<-betas_matrix<-NULL

setwd("C:/Users/vasileios palaskas/Desktop/Github folder/Bayesian_Variable_Selection_Volleyball/Ordered_TA_Skills")

T<-30000 # Total MCMC iterations
# Step 2  
for (i in 1:T){
  print(i)
  # Step 3: Data input needed for running the model through RStan.
  data_varsel<-list(Y=dataList$Y,X=dataList$X,
                    home_team=as.numeric(dataList$home_team),
                    away_team=as.numeric(dataList$away_team),n_teams=dataList$n_teams,
                    N=dataList$N,K=dataList$K,
                    ncat=6,gammas=gammas,post_mean_betas=post_mean_betas,
                    post_sd_betas=post_sd_betas)
  
  # Step 4:Run the model through RStan for one sampling iteration (20 warm up and 21 total iterations, 21-20=1 sampling iteration) in order to update the betas from the full conditional posterior distributions. 
  # Use the previous iteration's parameter values as initial parameter values so MCMC Algorithm can begin.
  ord_volley_skills_all<-stan("Ordered_BVS_TA_Skills.stan",
                              data=data_varsel,chains=1,
                              iter=21,warmup=20,init=list(list(betas=betas,temp_Intercept=temp_Intercept,gen_abil_raw=gen_abil_raw)),
                              control=list(adapt_window=15,adapt_init_buffer=3,adapt_term_buffer=2))
  
  # Initialize the log-likelihood for both cases 0 and 1 for gammas indicators/coefficients.
  log_point_zero<-matrix(NA,nrow=data_varsel$N,ncol=data_varsel$K) # matrix with log likelihoods when gamma[j]=0
  log_point_one<-matrix(NA,nrow=data_varsel$N,ncol=data_varsel$K)  #  matrix with log likelihoods when gamma[j]=1
  
  # Extract both model's parameters and log-likelihoods for both cases of gammas indicators.
  par<-extract(ord_volley_skills_all)
  temp_Intercept<-par$temp_Intercept[1,]
  betas<-par$betas[1,]
  gen_abil_raw<-par$gen_abil_raw[1,]
  log_point_zero<-par$log_lik_zero[1,,]
  log_point_one<-par$log_lik_one[1,,]
  
  # Prepare the vector of gammas for all candidate variables.
  gammas<-NULL
  
  # Step 5
  for (j in 1:data_varsel$K){# K candidate variables
    log_O_j<-O_j<-NULL
    
    # Step 6: Calculation of the logarithm (more convenient) of O_j quantity needed in order to update the gamma indicators.
    
    log_O_j<-sum(log_point_one[,j])-sum(log_point_zero[,j])+
      dnorm(betas[j],0,sqrt(data_varsel$N)*(post_sd_betas[j]),log=T)-
      dnorm(betas[j],post_mean_betas[j],post_sd_betas[j],log=T)
    # We specify an upper bound in order to avoid Nan potential problems due to overflow.
    log_O_j[log_O_j>700]<-700
    O_j<-exp(log_O_j)
    
    gammas<-c(gammas,rbinom(1,1,O_j/(1+O_j)))
    
  }
  
  # Step 7: In each one of T iterations, store the values of both posterior gammas and betas coefficients in the Txp matrices
  gammas_matrix<-c(gammas_matrix,gammas)
  betas_matrix<-c(betas_matrix,betas)
}
save(gammas_matrix,file="BVS_Ordered_TA_Skills_gammas")
save(betas_matrix,file="BVS_Ordered_TA_Skills_betas")
load("BVS_Ordered_TA_Skills_gammas")
load("BVS_Ordered_TA_Skills_betas")
# Save these values in order to manipulate them in terms of convergence diagnostics,, posterior summary statistics, etc...

# Store both gammas and betas posterior values after discarding the warmup from T iterations (here, we have chosen to discard the 20% of total T iterations).
warmup<-6000
# Each column includes the gammas values of each candidate variable.
final_posterior_values_gammas<-matrix(gammas_matrix[(dataList$K*warmup+1):length(gammas_matrix)],
                                      nrow=T-warmup,ncol=dataList$K,byrow=TRUE)
# Each column includes the gammas values of each candidate variable.
final_posterior_values_betas<-matrix(betas_matrix[(dataList$K*warmup+1):length(betas_matrix)],
                                     nrow=T-warmup,ncol=dataList$K,byrow=TRUE)
# Prepare a dataframe with column names the names of candidate variables.
df_final_posterior_values_gammas<-as.data.frame(final_posterior_values_gammas)
names(dataList$X)<-c("perfect serve","very good serve","failed serve"," perfect pass",
                     "very good pass","poor pass","failed pass","perfect att1","blocked att1",
                     "failed att1","perfect att2","blocked att2","failed att2","perfect block",
                     "block net violation","failed block","failed setting")
colnames(df_final_posterior_values_gammas)<-names(dataList$X)
# Step 8: Obtain the posterior inclusion probabilities for each one candidate variable
posterior_inclusion_probabilities<-round(apply(df_final_posterior_values_gammas,2,mean),3)
print(posterior_inclusion_probabilities)

### MCMC Convergence Checking
# 
# a) Firstly, for gammas and betas indicators
# 
## convert them to a mcmc pobject in terms of our convenience
mcmc_final_posterior_values_gammas<-as.mcmc(final_posterior_values_gammas)
colnames(mcmc_final_posterior_values_gammas)<-names(dataList$X)
mcmc_final_posterior_values_betas<-as.mcmc(final_posterior_values_betas)
colnames(mcmc_final_posterior_values_betas)<-names(dataList$X)




### MCMC Convergence Checking
# 
# a) Firstly, for gammas and betas indicators
# 
# convert them to a mcmc pobject in terms of our convenience
mcmc_final_posterior_values_gammas<-as.mcmc(final_posterior_values_gammas)
colnames(mcmc_final_posterior_values_gammas)<-names(dataList$X)
mcmc_final_posterior_values_betas<-as.mcmc(final_posterior_values_betas)
colnames(mcmc_final_posterior_values_betas)<-names(dataList$X)


####------GG-Plots for Convergence Diagnostics

#----Step 1: Convert the mcmc object to a ggmcmc object
gg_posterior_values_betas<- ggs(mcmc_final_posterior_values_betas)

gg_posterior_values_gammas<- ggs(mcmc_final_posterior_values_gammas)

#----Step2: Save in a single pdf all the necessary plots for the assessment of the convergence
ggmcmc(gg_posterior_values_betas, 
       file = "converg_betas_ordered_bvs_ta_skills.pdf", plot=c( "running",
                                                              "geweke","Rhat","autocorrelation"))



ggmcmc(gg_posterior_values_gammas, 
       file = "converg_gammas_ordered_bvs_ta_skills.pdf", plot=c( "running",
                                                               "geweke","Rhat","autocorrelation"))



# ###----Betas 
pdf(file="cumul_plots_betas_ordered_bvs_ta_skills.pdf", width =16, height =9)


par(mfrow=c(5,4))

cum_betas1<-cumsum(mcmc_final_posterior_values_betas[,1])/c(1:length(mcmc_final_posterior_values_betas[,1]))
plot(cum_betas1,cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",col="blue",xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[1])

cum_betas2<-cumsum(mcmc_final_posterior_values_betas[,2])/c(1:length(mcmc_final_posterior_values_betas[,2]))
plot(cum_betas2,type="l",cex=1.5,cex.lab=1.5,cex.axis=1.5,col="blue",xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[2])

cum_betas3<-cumsum(mcmc_final_posterior_values_betas[,3])/c(1:length(mcmc_final_posterior_values_betas[,3]))
plot(cum_betas3,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[3])

cum_betas4<-cumsum(mcmc_final_posterior_values_betas[,4])/c(1:length(mcmc_final_posterior_values_betas[,4]))
plot(cum_betas4,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[4])

cum_betas5<-cumsum(mcmc_final_posterior_values_betas[,5])/c(1:length(mcmc_final_posterior_values_betas[,5]))
plot(cum_betas5,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[5])

cum_betas6<-cumsum(mcmc_final_posterior_values_betas[,6])/c(1:length(mcmc_final_posterior_values_betas[,6]))
plot(cum_betas6,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[6])

cum_betas7<-cumsum(mcmc_final_posterior_values_betas[,7])/c(1:length(mcmc_final_posterior_values_betas[,7]))
plot(cum_betas7,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[7])

cum_betas8<-cumsum(mcmc_final_posterior_values_betas[,8])/c(1:length(mcmc_final_posterior_values_betas[,8]))
plot(cum_betas8,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[8])

cum_betas9<-cumsum(mcmc_final_posterior_values_betas[,9])/c(1:length(mcmc_final_posterior_values_betas[,9]))
plot(cum_betas9,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[9])

cum_betas10<-cumsum(mcmc_final_posterior_values_betas[,10])/c(1:length(mcmc_final_posterior_values_betas[,10]))
plot(cum_betas10,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[10])

cum_betas11<-cumsum(mcmc_final_posterior_values_betas[,11])/c(1:length(mcmc_final_posterior_values_betas[,11]))
plot(cum_betas11,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[11])

cum_betas12<-cumsum(mcmc_final_posterior_values_betas[,12])/c(1:length(mcmc_final_posterior_values_betas[,12]))
plot(cum_betas12,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[12])

cum_betas13<-cumsum(mcmc_final_posterior_values_betas[,13])/c(1:length(mcmc_final_posterior_values_betas[,13]))
plot(cum_betas13,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[13])

cum_betas14<-cumsum(mcmc_final_posterior_values_betas[,14])/c(1:length(mcmc_final_posterior_values_betas[,14]))
plot(cum_betas14,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[14])

cum_betas15<-cumsum(mcmc_final_posterior_values_betas[,15])/c(1:length(mcmc_final_posterior_values_betas[,15]))
plot(cum_betas15,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[15])

cum_betas16<-cumsum(mcmc_final_posterior_values_betas[,16])/c(1:length(mcmc_final_posterior_values_betas[,16]))
plot(cum_betas16,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[16])

cum_betas17<-cumsum(mcmc_final_posterior_values_betas[,17])/c(1:length(mcmc_final_posterior_values_betas[,17]))
plot(cum_betas17,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_betas)[17])

dev.off()




###----gammas home
pdf(file="cumul_plots_gammas_ordered_bvs_ta_skills.pdf", width =16, height =9)


par(mfrow=c(5,4))

cum_gammas1<-cumsum(mcmc_final_posterior_values_gammas[,1])/c(1:length(mcmc_final_posterior_values_gammas[,1]))
plot(cum_gammas1,cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",col="blue",xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[1])

cum_gammas2<-cumsum(mcmc_final_posterior_values_gammas[,2])/c(1:length(mcmc_final_posterior_values_gammas[,2]))
plot(cum_gammas2,type="l",cex=1.5,cex.lab=1.5,cex.axis=1.5,col="blue",xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[2])

cum_gammas3<-cumsum(mcmc_final_posterior_values_gammas[,3])/c(1:length(mcmc_final_posterior_values_gammas[,3]))
plot(cum_gammas3,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[3])

cum_gammas4<-cumsum(mcmc_final_posterior_values_gammas[,4])/c(1:length(mcmc_final_posterior_values_gammas[,4]))
plot(cum_gammas4,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[4])

cum_gammas5<-cumsum(mcmc_final_posterior_values_gammas[,5])/c(1:length(mcmc_final_posterior_values_gammas[,5]))
plot(cum_gammas5,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[5])

cum_gammas6<-cumsum(mcmc_final_posterior_values_gammas[,6])/c(1:length(mcmc_final_posterior_values_gammas[,6]))
plot(cum_gammas6,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[6])

cum_gammas7<-cumsum(mcmc_final_posterior_values_gammas[,7])/c(1:length(mcmc_final_posterior_values_gammas[,7]))
plot(cum_gammas7,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[7])

cum_gammas8<-cumsum(mcmc_final_posterior_values_gammas[,8])/c(1:length(mcmc_final_posterior_values_gammas[,8]))
plot(cum_gammas8,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[8])

cum_gammas9<-cumsum(mcmc_final_posterior_values_gammas[,9])/c(1:length(mcmc_final_posterior_values_gammas[,9]))
plot(cum_gammas9,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[9])

cum_gammas10<-cumsum(mcmc_final_posterior_values_gammas[,10])/c(1:length(mcmc_final_posterior_values_gammas[,10]))
plot(cum_gammas10,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[10])

cum_gammas11<-cumsum(mcmc_final_posterior_values_gammas[,11])/c(1:length(mcmc_final_posterior_values_gammas[,11]))
plot(cum_gammas11,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[11])

cum_gammas12<-cumsum(mcmc_final_posterior_values_gammas[,12])/c(1:length(mcmc_final_posterior_values_gammas[,12]))
plot(cum_gammas12,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[12])

cum_gammas13<-cumsum(mcmc_final_posterior_values_gammas[,13])/c(1:length(mcmc_final_posterior_values_gammas[,13]))
plot(cum_gammas13,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[13])

cum_gammas14<-cumsum(mcmc_final_posterior_values_gammas[,14])/c(1:length(mcmc_final_posterior_values_gammas[,14]))
plot(cum_gammas14,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[14])

cum_gammas15<-cumsum(mcmc_final_posterior_values_gammas[,15])/c(1:length(mcmc_final_posterior_values_gammas[,15]))
plot(cum_gammas15,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[15])

cum_gammas16<-cumsum(mcmc_final_posterior_values_gammas[,16])/c(1:length(mcmc_final_posterior_values_gammas[,16]))
plot(cum_gammas16,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[16])

cum_gammas17<-cumsum(mcmc_final_posterior_values_gammas[,17])/c(1:length(mcmc_final_posterior_values_gammas[,17]))
plot(cum_gammas17,type="l",col="blue",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlab="Iterations",ylab=colnames(mcmc_final_posterior_values_gammas)[17])

dev.off()




## NOTE (Leo_Egidi): I ran the model with more chains to check the Gelman-Rubin statistic
##       R_hat and the effective sample size and everything works well!




