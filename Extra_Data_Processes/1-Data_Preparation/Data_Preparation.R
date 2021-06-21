# Import libraries for data manipulation
library(sqldf)
library(dplyr)
#Read and further processing of the data

#
volley<-read.csv(file.choose(),
                 header=T)#Load Data\Final_Regular_Season_data.csv
head(volley,30)
colnames(volley)[3]<-"home_team"
colnames(volley)[names(volley)=="Team.B"] <-"away_team"
names(volley)
str(volley)
volley$home_team<-factor(volley$home_team)
volley$away_team<-factor(volley$away_team)
dim(volley)###494 obs and 52 variab


#New dataframe with corresponding team per set
# and the corresponding score (for Binomial per set)

# Here we make binary variable results_set
# win set for  home team :1
# win set for away team :0
results_set<- c()
results_set[volley$Points>volley$Points.1]<- 1
results_set[volley$Points <volley$Points.1]<- 0
teams <- unique(volley$home_team)

team1 <- match(volley$home_team,teams)
team2 <- match(volley$away_team, teams)
nteams <- length(teams)###number of teams in the league : 12
cont <-1
match_index <-c()#matchindex=matchday which refers to each match event separately
match_index[1] <-1
for (n in 2:length(volley$DAY)){
  if(team1[n]==team1[n-1] ){
    match_index[n]<-cont
  }else{
    cont <-cont+1
    match_index[n] <- cont
  }
}
match_index
attach(volley)


# Proper processing of the data
new_volley<-data.frame(Week=DAY,match_day=match_index,set=Set,volley[,-c(1,2)],results_set=results_set)
detach(volley)
head(new_volley,30)


# Rename properly both home and away skill events
colnames(new_volley)[c(6:28)]<-c(
  "Home_total_serves","Home_perfect_serves","Home_very_good_serves",
  "Home_failed_serves","Home_total_passes","Home_perfect_passes",
  "Home_very_good_passes","Home_poor_passes","Home_failed_passes",
  "Home_total_att1","Home_perfect_att1",
  "Home_blocked_att1","Home_failed_att1","Home_total_att2",
  "Home_perfect_att2","Home_blocked_att2","Home_failed_att2",
  "Home_total_blocks", "Home_perfect_blocks", "Home_net_violation_blocks",
  "Home_failed_blocks","Home_total_settings","Home_failed_settings")

colnames(new_volley)[c(6:28)+25]<-c(
  "Away_total_serves","Away_perfect_serves","Away_very_good_serves",
  "Away_failed_serves","Away_total_passes","Away_perfect_passes",
  "Away_very_good_passes","Away_poor_passes","Away_failed_passes",
  "Away_total_att1","Away_perfect_att1",
  "Away_blocked_att1","Away_failed_att1","Away_total_att2",
  "Away_perfect_att2","Away_blocked_att2","Away_failed_att2",
  "Away_total_blocks", "Away_perfect_blocks", "Away_net_violation_blocks",
  "Away_failed_blocks","Away_total_settings","Away_failed_settings")


#-----Store only the skills, by each match
new_volley_skills<-new_volley[,c(2,c(6:28),c(6:28)+25)]

new_volley_skills<-aggregate(. ~ match_day,new_volley_skills, sum)
# Creating the following variables
# home score ,team ,skills and away score, team , skills 
# per match (132 in total)

home_score<-away_score<-home_Team<-away_Team<-NULL
home_points<-away_points<-NULL
levels(volley$home_team)
for (i in 1:dim(new_volley_skills)[1]) {
 home_score[i]<-length(new_volley$results_set[new_volley$match_day==i & new_volley$results_set==1])
 away_score[i]<-length(new_volley$results_set[new_volley$match_day==i & new_volley$results_set==0])
 home_Team[i]<-as.vector(new_volley$home_team[new_volley$match_day==i&new_volley$set==1])
 away_Team[i]<-as.vector(new_volley$away_team[new_volley$match_day==i&new_volley$set==1])
 home_points[i]<-sum(new_volley[,5][new_volley$match_day==i ])
 away_points[i]<-sum(new_volley[,30][new_volley$match_day==i ])
}

# A dataframe (without technical details like passes, serv,...) with only
# home and away teams as well as their corresponding sets scores ####
data_by_sets<-data.frame(new_volley_skills,home_Team,away_Team,
                         home_sets=home_score,away_sets=away_score,
                         home_points,away_points)
head(data_by_sets)



# colnames(data_by_sets)[c(2:24)]<-c(
#   "Home_total_serves","Home_perfect_serves","Home_very_good_serves",
#   "Home_failed_serves","Home_total_passes","Home_perfect_passes",
#   "Home_very_good_passes","Home_poor_passes","Home_failed_passes",
#   "Home_total_att1","Home_perfect_att1",
#   "Home_blocked_att1","Home_failed_att1","Home_total_att2",
#   "Home_perfect_att2","Home_blocked_att2","Home_failed_att2",
#   "Home_total_blocks", "Home_perfect_blocks", "Home_net_violation_blocks",
#   "Home_failed_blocks","Home_total_settings","Home_failed_settings")
# 
# colnames(data_by_sets)[c(25:47)]<-c(
#   "Away_total_serves","Away_perfect_serves","Away_very_good_serves",
#   "Away_failed_serves","Away_total_passes","Away_perfect_passes",
#   "Away_very_good_passes","Away_poor_passes","Away_failed_passes",
#   "Away_total_att1","Away_perfect_att1",
#   "Away_blocked_att1","Away_failed_att1","Away_total_att2",
#   "Away_perfect_att2","Away_blocked_att2","Away_failed_att2",
#   "Away_total_blocks", "Away_perfect_blocks", "Away_net_violation_blocks",
#   "Away_failed_blocks","Away_total_settings","Away_failed_settings")


# colnames(data_by_sets)[c(2:24)]<-c(
#   "(Home) total serves","(Home) perfect serves","(Home) very good serves",
#   "(Home) failed serves","(Home) total passes","(Home) perfect passes",
#   "(Home) very good passes","(Home) poor passes","(Home) failed passes",
#   "(Home) total att1","(Home) perfect att1",
#   "(Home) blocked att1","(Home) failed att1","(Home) total att2",
#   "(Home) perfect att2","(Home) blocked att2","(Home) failed att1",
#   "(Home) total blocks", "(Home) perfect blocks", "(Home) net violation blocks",
#   "(Home) failed blocks","(Home) total settings","(Home) failed settings")
# 
# colnames(data_by_sets)[c(25:47)]<-c(
#   "(Away) total serves","(Away) perfect serves","(Away) very good serves",
#   "(Away) failed serves","(Away) total passes","(Away) perfect passes",
#   "(Away) very good passes","(Away) poor passes","(Away) failed passes",
#   "(Away) total att1","(Away) perfect att1",
#   "(Away) blocked att1","(Away) failed att1","(Away) total att2",
#   "(Away) perfect att2","(Away) blocked att2","(Away) failed att1",
#   "(Away) total blocks", "(Away) perfect blocks", "(Away) net violation blocks",
#   "(Away) failed blocks","(Away) total settings","(Away) failed settings")

data_by_sets$points_difference<-data_by_sets$home_points-data_by_sets$away_points
data_by_sets$sets_difference<-data_by_sets$home_sets-data_by_sets$away_sets

#---Levels
data_by_sets$home_Team<-factor(data_by_sets$home_Team)
data_by_sets$away_Team<-factor(data_by_sets$away_Team)

levels(data_by_sets$home_Team)
levels(data_by_sets$away_Team)
head(data_by_sets)



