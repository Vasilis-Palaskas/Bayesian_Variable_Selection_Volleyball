

#Read and further processing of the data

#
volley<-read.csv(file.choose(),header=T)#Load Data\Final_Regular_Season_data.csv
head(volley,30)
is.data.frame(volley)
colnames(volley)[3]<-"home_team"
colnames(volley)[names(volley)=="Team.B"] <-"away_team"
names(volley)
str(volley)
volley$home_team<-as.factor(volley$home_team)
volley$away_team<-as.factor(volley$away_team)
View(volley)
str(volley)
dim(volley)###494 obs and 52 variab


#New dataframe with corresponding team per set
# and the corresponding score (for Binomial per set)

#Here we make binary variable results_set
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
match_index <-c()
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
names(new_volley)
dim(new_volley)


# Rename properly both home and away skill events
colnames(new_volley)[c(7,8,9,11,12,13,14,16,17,18,20,21,22,24,25,26,28)]<-c("serv_perfect","serv_vg","serv_error",
"pass_perfect","pass_vg","pass_poor","pass_error","att1_perfect",
"att1_blocked","att1_error","att2_perfect","att2_blocked","att2_error",
"block_perfect","block_net_violation","block_error","setting_error")

colnames(new_volley)[c(7,8,9,11,12,13,14,16,17,18,20,21,22,24,25,26,28)+25]<-c("serv_perfect","serv_vg","serv_error",
"pass_perfect","pass_vg","pass_poor","pass_error","att1_perfect",
"att1_blocked","att1_error","att2_perfect","att2_blocked","att2_error",
"block_perfect","block_net_violation","block_error","setting_error")



# Creating the following variables
# home score ,team ,skills and away score, team , skills 
# per match (132 in total)


home_score<-away_score<-home_Team<-away_Team<-NULL


home_serv_perfect<-home_serv_vg<-home_serv_error<-NULL
home_pass_perfect<-home_pass_vg<-home_pass_poor<-home_pass_error<-NULL
home_att1_perfect<-home_att1_blocked<-home_att1_error<-NULL
home_att2_perfect<-home_att2_blocked<-home_att2_error<-NULL
home_block_perfect<-home_block_net_violation<-home_block_error<-NULL
home_setting_error<-NULL

away_serv_perfect<-away_serv_vg<-away_serv_error<-NULL
away_pass_perfect<-away_pass_vg<-away_pass_poor<-away_pass_error<-NULL
away_att1_perfect<-away_att1_blocked<-away_att1_error<-NULL
away_att2_perfect<-away_att2_blocked<-away_att2_error<-NULL
away_block_perfect<-away_block_net_violation<-away_block_error<-NULL
away_setting_error<-NULL
home_points<-away_points<-NULL
levels(volley$home_team)
for (i in 1:132) {
 home_score[i]<-length(new_volley$results_set[new_volley$match_day==i & new_volley$results_set==1])
 away_score[i]<-length(new_volley$results_set[new_volley$match_day==i & new_volley$results_set==0])
 home_Team[i]<-as.vector(new_volley$home_team[new_volley$match_day==i&new_volley$set==1])
 away_Team[i]<-as.vector(new_volley$away_team[new_volley$match_day==i&new_volley$set==1])
 
 home_serv_perfect[i]<-sum(new_volley[,7][new_volley$match_day==i])
 home_serv_vg[i]<-sum(new_volley[,8][new_volley$match_day==i ])
 home_serv_error[i]<-sum(new_volley[,9][new_volley$match_day==i ])
 home_pass_perfect[i]<-sum(new_volley[,11][new_volley$match_day==i ])
 home_pass_vg[i]<-sum(new_volley[,12][new_volley$match_day==i ])
 home_pass_poor[i]<-sum(new_volley[,13][new_volley$match_day==i ])
 home_pass_error[i]<-sum(new_volley[,14][new_volley$match_day==i ])
 home_att1_perfect[i]<-sum(new_volley[,16][new_volley$match_day==i ])
 home_att1_blocked[i]<-sum(new_volley[,17][new_volley$match_day==i ])
 home_att1_error[i]<-sum(new_volley[,18][new_volley$match_day==i ])
 home_att2_perfect[i]<-sum(new_volley[,20][new_volley$match_day==i ])
 home_att2_blocked[i]<-sum(new_volley[,21][new_volley$match_day==i ])
 home_att2_error[i]<-sum(new_volley[,22][new_volley$match_day==i ])
 home_block_perfect[i]<-sum(new_volley[,24][new_volley$match_day==i ])
 home_block_net_violation[i]<-sum(new_volley[,25][new_volley$match_day==i ])
 home_block_error[i]<-sum(new_volley[,26][new_volley$match_day==i ])
 home_setting_error[i]<-sum(new_volley[,28][new_volley$match_day==i ])

 away_serv_perfect[i]<-sum(new_volley[,32][new_volley$match_day==i ])
 away_serv_vg[i]<-sum(new_volley[,33][new_volley$match_day==i ])
 away_serv_error[i]<-sum(new_volley[,34][new_volley$match_day==i ])
 away_pass_perfect[i]<-sum(new_volley[,36][new_volley$match_day==i ])
 away_pass_vg[i]<-sum(new_volley[,37][new_volley$match_day==i ])
 away_pass_poor[i]<-sum(new_volley[,38][new_volley$match_day==i ])
 away_pass_error[i]<-sum(new_volley[,39][new_volley$match_day==i ])
 away_att1_perfect[i]<-sum(new_volley[,41][new_volley$match_day==i ])
 away_att1_blocked[i]<-sum(new_volley[,42][new_volley$match_day==i ])
 away_att1_error[i]<-sum(new_volley[,43][new_volley$match_day==i ])
 away_att2_perfect[i]<-sum(new_volley[,45][new_volley$match_day==i ])
 away_att2_blocked[i]<-sum(new_volley[,46][new_volley$match_day==i ])
 away_att2_error[i]<-sum(new_volley[,47][new_volley$match_day==i ])
 away_block_perfect[i]<-sum(new_volley[,49][new_volley$match_day==i ])
 away_block_net_violation[i]<-sum(new_volley[,50][new_volley$match_day==i ])
 away_block_error[i]<-sum(new_volley[,51][new_volley$match_day==i ])
 away_setting_error[i]<-sum(new_volley[,53][new_volley$match_day==i ])
 home_points[i]<-sum(new_volley[,5][new_volley$match_day==i ])
 away_points[i]<-sum(new_volley[,30][new_volley$match_day==i ])
}

# A dataframe (without technical details like passes, serv,...) with only
# home and away teams as well as their corresponding sets scores ####


datafr_teams_scores_set<-data.frame(home_Team,away_Team,home_score,away_score)
names(datafr_teams_scores_set)
levels(datafr_teams_scores_set$away_Team)


#dim(datafr_teams_scores_set)
#head(datafr_teams_scores_set,132)

# A dataframe (without technical details like passes, serv,...) and with
# home , away teams and correspondings points scores. 

datafr_teams_scores_points<-data.frame(home_Team,away_Team,points_difference=home_points-away_points)
names(datafr_teams_scores_points)
levels(datafr_teams_scores_set$away_Team)
dim(datafr_teams_scores_points)
head(datafr_teams_scores_points,132)
points_132<-datafr_teams_scores_points
colnames(points_132)[3]<-"Point difference"



datafr_teams_scores_set_skills<-data.frame(datafr_teams_scores_set,
home_serv_perfect,home_serv_vg,home_serv_error,
home_pass_perfect,home_pass_vg,home_pass_poor,home_pass_error,
home_att1_perfect,home_att1_blocked,home_att1_error,
home_att2_perfect,home_att2_blocked,home_att2_error,
home_block_perfect,home_block_net_violation,home_block_error,home_setting_error,
away_serv_perfect,away_serv_vg,away_serv_error,
away_pass_perfect,away_pass_vg,away_pass_poor,away_pass_error,
away_att1_perfect,away_att1_blocked,away_att1_error,
away_att2_perfect,away_att2_blocked,away_att2_error,
away_block_perfect,away_block_net_violation,away_block_error,away_setting_error)


