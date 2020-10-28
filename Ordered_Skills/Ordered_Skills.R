# Load the proper libraries.
library(rstan)
library(coda)

# Choose the working directory of this file (...\\Submitted_Appendix\\Ordered\\)

setwd("C:/Users/vasileios palaskas/Desktop/BVS_Paper/Ordered_TA_Skills")

# Load the properly prepared data ("Data_ordered_skills").
load("datalist_ordered")

head(dataList)

#Numerize the factors in terms of your convenience
dataList<-list(Y=dataList$Y,X=dataList$X,n_teams=12,
      N=dataList$N,K=ncol(dataList$X),ncat=6)



## Run Full_ordered_skills.stan
Full_ordered_skills<-stan(file.choose(),iter=12000, warmup=2000,chains=1,thin=2,
                          data=dataList,control=list(max_treedepth=15),cores=1)

save(Full_ordered_skills,file="Full_ordered_skills")
# Load the results from the full ordered logistic model (with all candidate variables).
load(file="Full_ordered_skills")

# Extract the posterior summary statistics of both candidate variables' parameters and rest of other parameters.


betas_summary<-summary(Full_ordered_skills,pars = c("beta"))$summary
intercept_summary<-summary(Full_ordered_skills, pars = c("temp_Intercept"))$summary

# Use their posterior means and standard deviations for both initial values specification and prior specification.

post_mean_betas<-betas_summary[1:dataList$K,1]
post_sd_betas<-betas_summary[1:dataList$K,3]

# Step 1: Initialization of the model parameters.

gammas<-rep(1,dataList$K)# All the candidate variables included in the model
temp_Intercept<-intercept_summary[,1]
betas<-post_mean_betas

# Prepare the vectors with the posterior samples of dimension Txp (p=K during algorithm iterations) for all gammas and betas coefficients , respectively.
gammas_matrix<-betas_matrix<-NULL


T<-30000 # Total MCMC iterations
# Step 2  
for (i in 1:T){
  print(i)
  # Step 3: Data input needed for running the model through RStan.
  data_varsel<-list(Y=dataList$Y,X=dataList$X,
                    N=dataList$N,K=dataList$K,
                    ncat=6,gammas=gammas,post_mean_betas=post_mean_betas,
                    post_sd_betas=post_sd_betas)
  
  # Step 4:Run the model through RStan for one sampling iteration (20 warm up and 21 total iterations, 21-20=1 sampling iteration) in order to update the betas from the full conditional posterior distributions. 
  # Use the previous iteration's parameter values as initial parameter values so MCMC Algorithm can begin.
  ord_volley_skills_all<-stan("Ordered_BVS_Skills.stan",
                              data=data_varsel,chains=1,
                              iter=21,warmup=20,init=list(list(betas=betas,temp_Intercept=temp_Intercept)),
                              control=list(adapt_window=15,adapt_init_buffer=3,adapt_term_buffer=2))
  
  # Initialize the log-likelihood for both cases 0 and 1 for gammas indicators/coefficients.
  log_point_zero<-matrix(NA,nrow=data_varsel$N,ncol=data_varsel$K) # matrix with log likelihoods when gamma[j]=0
  log_point_one<-matrix(NA,nrow=data_varsel$N,ncol=data_varsel$K)  #  matrix with log likelihoods when gamma[j]=1
  
  # Extract both model's parameters and log-likelihoods for both cases of gammas indicators.
  par<-extract(ord_volley_skills_all)
  temp_Intercept<-par$temp_Intercept[1,]
  betas<-par$betas[1,]
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

# Save these values in order to manipulate them in terms of convergence diagnostics,, posterior summary statistics, etc...
save(gammas_matrix,file="BVS_Ordered_Skills_gammas")
save(betas_matrix,file="BVS_Ordered_Skills_betas")

# Store both gammas and betas posterior values after discarding the warmup from T iterations (here, we have chosen to discard the 20% of total T iterations).
warmup<-6000
# Each column includes the gammas values of each candidate variable.
final_posterior_values_gammas<-matrix(gammas_matrix[(dataList$K*warmup+1):length(gammas_matrix)],
                                      nrow=T-warmup,ncol=dataList$K,byrow=TRUE)
# Each column includes the gammas values of each candidate variable.
final_posterior_values_betas<-matrix(betas_matrix[(dataList$K*warmup+1):length(betas_matrix)],
                                     nrow=T-warmup,ncol=dataList$K,byrow=TRUE)
# Prepare a dataframe by assigning in the variables names the correspodnding column names.
df_final_posterior_values_gammas<-as.data.frame(final_posterior_values_gammas)
colnames(df_final_posterior_values_gammas)<-names(dataList_order_skills$X)
# Step 8: Obtain the posterior inclusion probabilities for each one candidate variable
posterior_inclusion_probabilities<-round(apply(df_final_posterior_values_gammas,2,mean),3)
print(posterior_inclusion_probabilities)


# MCMC Convergence diagnostics







