
## Posterior summary diagnostics

## Stan interface for both summary and convergence diagnostics
launch_shinystan(ZDTS_only_Skills_after_BVS)
##----Parameters Names

# teams_index <- unique(dataList$home_team)
sims <- rstan::extract(ZDTS_only_Skills_after_BVS)

beta_home <- sims$beta_home
beta_away <- sims$beta_away
beta_home_away<-cbind(beta_home,beta_away)
mu<- sims$mu
home<- sims$home


## Order of ability parameters (based on the posterior means)
beta_home_hat <- apply(beta_home,2, median)
beta_away_hat <- apply(beta_away,2, median)
beta_home_away_hat <- apply(beta_home_away,2, median)

mu_hat <- apply(mu,1,median)
home_hat <- apply(home,1,median)

beta_home_hat_ord <- order(beta_home_hat, decreasing = TRUE)
beta_away_hat_ord <- order(beta_away_hat, decreasing = TRUE)
beta_home_away_hat_ord <- order(beta_home_away_hat, decreasing = TRUE)

mu_hat_ord <- order(mu_hat, decreasing = TRUE)
home_hat_ord <- order(home_hat, decreasing = TRUE)

#----------ZDTS TA +Skills
##---Parameters Names


## Data frame of parameters in terms of convenience in both tables and graphs

beta_home<-data.frame(beta_home)
beta_away<-data.frame(beta_away)
beta_home_away<-data.frame(beta_home_away)
mu<-data.frame(mu)
home<-data.frame(home)


##---Proper parameters renaming
colnames(beta_home)<-colnames(final_X_home_std)
colnames(beta_away)<-colnames(final_X_away_std)
colnames(beta_home_away)<-c(colnames(final_X_home_std),
                            colnames(final_X_away_std))

###-----MCMC Posterior 95% uncertainty intervals

##---Combine 2 plots
color_scheme_set("brightblue")

pdf(file="ZDTS_only_Skills_betas.pdf", width =12, height =7.5)

mcmc_intervals(beta_home_away[,c(beta_home_away_hat_ord)],
                      prob = 0.95,prob_outer=0.95,
                      point_est = c( "mean"))+ggtitle("Skill Events")+xlim(-2,1.4)+
  scale_x_continuous(breaks = c(seq(from = -2, to = 1.4, by = 0.4)))+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))
dev.off()
##--Separate 2 plots
color_scheme_set("brightblue")


plot1<-mcmc_intervals(beta_home[,c(beta_home_hat_ord)],
               prob = 0.95,prob_outer=0.95,
               point_est = c( "mean"))+ggtitle("Home Skill Events")+xlim(-2,1.4)+
  scale_x_continuous(breaks = c(seq(from = -2, to = 1.4, by = 0.4)))+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))


plot2<-mcmc_intervals(beta_away[,c(beta_away_hat_ord)],
                      prob = 0.95,prob_outer=0.95,
                      point_est = c( "mean"))+ggtitle("Away Skill Events")+xlim(-2,1.4)+
  scale_x_continuous(breaks = c(seq(from = -2, to = 1.4, by = 0.4)))+
  theme(axis.text.x = element_text( size = 23, angle = 0, hjust = .5, vjust = .5),
        axis.text.y = element_text( size = 23, angle = 0, hjust = 1, vjust = 0),  
        axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust= 0),
        plot.title  =element_text( size = 20))

pdf(file="ZDTS_only_Skills_betas_separate.pdf", width =12, height =7.5)
grid.arrange(plot1,plot2,ncol=2)
dev.off()
