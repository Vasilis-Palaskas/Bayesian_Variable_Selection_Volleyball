# Load packages
library(latex2exp)
library(shape)
library(skellam)

####-------ZDTS Distr Plot
l1<-c(seq(1,50,by=7),seq(52,902,by=50))
l2<-c(seq(0,50,by=10),seq(51,901,by=50))

#l1-l2 calculation across values
l1_l2_diff<-NULL
zdts_dist<-NULL
zdts_support<-c(-3,-2,-1,1,2,3)

for (i in l1){
	for (j in l2){
		l1_l2_diff<-c(l1_l2_diff,i-j)
		zdts_dist<-c(zdts_dist,
			dskellam(zdts_support,i,j)/sum(dskellam(zdts_support,i,j))
			)
		

	}
}

##----Dataframe of the ZDTS distribution
	zdts_dist_matrix<-matrix(zdts_dist,ncol=length(zdts_support),byrow=T)# matrix of zdts acrtoss  all support values
zdts_dist_df<-data.frame(l1_l2_difference=l1_l2_diff,
		zdts_distribution=zdts_dist_matrix)
colnames(zdts_dist_df)<-c("l1_l2_Difference","zdts_minus_3","zdts_minus_2",
		"zdts_minus_1","zdts_1","zdts_2","zdts_3")

zdts_dist_df[is.na(zdts_dist_df)]<--3
#####-----------Plot of this relationship between l1-l2 and ZDTS Distribution
par(mfrow=c(2,3))
##----PLOT 1 F(1)
plot( zdts_dist_df$l1_l2_Difference,zdts_dist_df$zdts_1, col=2, lwd=2, lty=2,
		xlab=TeX('$\\lambda_{1}-\\lambda_{2}$'),ylim=c(-4,1),main=TeX('ZDTS Likelihood'),
	ylab=TeX('$f(1)_{ZDTS}$')	, xaxt='n',yaxt="n",cex.axis=2,cex.lab=2,cex.main=2)
# Now, define a custom axis
axis(side = 1, at=c(seq(min(zdts_dist_df$l1_l2_Difference),max(zdts_dist_df$l1_l2_Difference),by=50)),cex.axis=1.6 ,cex.lab=1.6)
axis(side = 2, at=c(seq(-3.1,-0.1,by=1),seq(0,1,by=0.1)),cex.axis=1.6 ,cex.lab=1.6)

##----PLOT 2 F(2)

plot( zdts_dist_df$l1_l2_Difference,zdts_dist_df$zdts_2, col=3, lwd=2, lty=2,
		xlab=TeX('$\\lambda_{1}-\\lambda_{2}$'),ylim=c(-4,1),main=TeX('ZDTS Likelihood'),
	ylab=TeX('$f(2)_{ZDTS}$')	, xaxt='n',yaxt="n",cex.axis=2,cex.lab=2,cex.main=2)
# Now, define a custom axis
axis(side = 1, at=c(seq(min(zdts_dist_df$l1_l2_Difference),max(zdts_dist_df$l1_l2_Difference),by=50)),cex.axis=1.6 ,cex.lab=1.6)
axis(side = 2, at=c(seq(-3.1,-0.1,by=1),seq(0,1,by=0.1)),cex.axis=1.6 ,cex.lab=1.6)
##----PLOT 3 F(3)

plot( zdts_dist_df$l1_l2_Difference,zdts_dist_df$zdts_3, col=4, lwd=2, lty=2,
		xlab=TeX('$\\lambda_{1}-\\lambda_{2}$'),ylim=c(-4,1),main=TeX('ZDTS Likelihood'),
	ylab=TeX('$f(3)_{ZDTS}$')	, xaxt='n',yaxt="n",cex.axis=2,cex.lab=2,cex.main=2)
# Now, define a custom axis
axis(side = 1, at=c(seq(min(zdts_dist_df$l1_l2_Difference),max(zdts_dist_df$l1_l2_Difference),by=50)),cex.axis=1.6 ,cex.lab=1.6)
axis(side = 2, at=c(seq(-3.1,-0.1,by=1),seq(0,1,by=0.1)),cex.axis=1.6 ,cex.lab=1.6)




##----PLOT 4 F(-1)
plot( zdts_dist_df$l1_l2_Difference,zdts_dist_df$zdts_minus_1, col=5, lwd=2, lty=2,
		xlab=TeX('$\\lambda_{1}-\\lambda_{2}$'),ylim=c(-4,1),main=TeX('ZDTS Likelihood'),
	ylab=TeX('$f(-1)_{ZDTS}$')	, xaxt='n',yaxt="n",cex.axis=2,cex.lab=2,cex.main=2)
# Now, define a custom axis
axis(side = 1, at=c(seq(min(zdts_dist_df$l1_l2_Difference),max(zdts_dist_df$l1_l2_Difference),by=50)),cex.axis=1.6 ,cex.lab=1.6)
axis(side = 2, at=c(seq(-3.1,-0.1,by=1),seq(0,1,by=0.1)),cex.axis=1.6 ,cex.lab=1.6)

##----PLOT 5 F(-2)

plot(zdts_dist_df$l1_l2_Difference,zdts_dist_df$zdts_minus_2, col=6, lwd=2, lty=2,
		xlab=TeX('$\\lambda_{1}-\\lambda_{2}$'),ylim=c(-4,1),main=TeX('ZDTS Likelihood'),
	ylab=TeX('$f(-2)_{ZDTS}$')	, xaxt='n',yaxt="n",cex.axis=2,cex.lab=2,cex.main=2)
# Now, define a custom axis
axis(side = 1, at=c(seq(min(zdts_dist_df$l1_l2_Difference),max(zdts_dist_df$l1_l2_Difference),by=50)),cex.axis=1.6 ,cex.lab=1.6)
axis(side = 2, at=c(seq(-3.1,-0.1,by=1),seq(0,1,by=0.1)),cex.axis=1.6 ,cex.lab=1.6)
##----PLOT 6 F(-3)

plot( zdts_dist_df$l1_l2_Difference,zdts_dist_df$zdts_minus_3, col=7, lwd=2, lty=2,
		xlab=TeX('$\\lambda_{1}-\\lambda_{2}$'),ylim=c(-4,1),main=TeX('ZDTS Likelihood'),
	ylab=TeX('$f(-3)_{ZDTS}$')	, xaxt='n',yaxt="n",cex.axis=2,cex.lab=2,cex.main=2)
# Now, define a custom axis
axis(side = 1, at=c(seq(min(zdts_dist_df$l1_l2_Difference),max(zdts_dist_df$l1_l2_Difference),by=50)),cex.axis=1.6 ,cex.lab=1.6)
axis(side = 2, at=c(seq(-3.1,-0.1,by=1),seq(0,1,by=0.1)),cex.axis=1.6 ,cex.lab=1.6)

