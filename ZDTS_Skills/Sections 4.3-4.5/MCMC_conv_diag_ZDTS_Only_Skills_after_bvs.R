####------GG-Plots for Convergence Diagnostics

#----Step 1: Convert the mcmc object to a ggmcmc object
ZDTS_only_Skills_after_BVS_parameters<-ggs(mcmc(cbind(mu,home,
                                                      beta_home_away)))
gg_ZDTS_only_Skills_after_BVS_beta_home<- ggs(ZDTS_only_Skills_after_BVS,
                                              family = "beta_home")
gg_ZDTS_only_Skills_after_BVS_beta_away<- ggs(ZDTS_only_Skills_after_BVS,
                                              family = "beta_away")

gg_ZDTS_only_Skills_after_BVS_mu<- ggs(ZDTS_only_Skills_after_BVS,
                                                      family = "mu")
gg_ZDTS_only_Skills_after_BVS_home<- ggs(ZDTS_only_Skills_after_BVS,
                                                            family = "home")

#----Step2: Save in a single pdf all the necessary plots for the assessment of the convergence


ggmcmc(ZDTS_only_Skills_after_BVS_parameters, 
       file = "converg_gg_ZDTS_only_Skills_after_BVS_parameters.pdf", 
       plot=c( "running","traceplot", "geweke","Rhat","autocorrelation"))

#----Separately across multiple chains
ggmcmc(gg_ZDTS_only_Skills_after_BVS_home,
       file = "converg_home_ZDTS_only_Skills_after_BVS.pdf", plot=c( "running","traceplot",
                                                                              "geweke","Rhat","autocorrelation"))


ggmcmc(gg_ZDTS_only_Skills_after_BVS_mu,
       file = "converg_mu_ZDTS_only_Skills_after_BVS.pdf", plot=c( "running","traceplot",
                                                               "geweke","Rhat","autocorrelation"))

ggmcmc(gg_ZDTS_only_Skills_after_BVS_beta_home,
       file = "converg_ZDTS_only_Skills_after_BVS_beta_home.pdf", plot=c( "running","traceplot",
                                                                    "geweke","Rhat","autocorrelation"))

ggmcmc(gg_ZDTS_only_Skills_after_BVS_beta_away,
       file = "converg_ZDTS_only_Skills_after_BVS_beta_away.pdf", plot=c( "running","traceplot",
                                                                          "geweke","Rhat","autocorrelation"))
