#------------------------------------Revised Data Processing------------------------------------
#-----Home transformations
data_by_sets$Home_precise_passes<-data_by_sets[,"Home_very_good_passes"]+
  data_by_sets[,"Home_perfect_passes"]

data_by_sets$Home_modetate_passes<-data_by_sets[,"Home_total_passes"]-data_by_sets$Home_precise_passes-
  data_by_sets[,"Home_poor_passes"]-data_by_sets[,"Home_failed_passes"]

data_by_sets$Home_all_failed_blocks<-data_by_sets[,"Home_net_violation_blocks"]+data_by_sets[,"Home_failed_blocks"]


#-----away transformations
data_by_sets$Away_precise_passes<-data_by_sets[,"Away_very_good_passes"]+
  data_by_sets[,"Away_perfect_passes"]

data_by_sets$Away_modetate_passes<-data_by_sets[,"Away_total_passes"]-data_by_sets$Away_precise_passes-
  data_by_sets[,"Away_poor_passes"]-data_by_sets[,"Away_failed_passes"]

data_by_sets$Away_all_failed_blocks<-data_by_sets[,"Away_net_violation_blocks"]+data_by_sets[,"Away_failed_blocks"]

#---Renames in terms of convenience
# colnames(data_by_sets)[colnames(data_by_sets)%in%c("Home_precise_passes",
#                                                    "Home_modetate_passes",
#                                                   "Home_failed_blocks",
#                                                   "Away_precise_passes",
#                                                   "Away_modetate_passes",
#                                                   "Away_failed_blocks" )]<-c(
#                                                     "(Home) precise passes",
#                                                     "(Home) modetate passes",
#                                                     "(Home) failed blocks",
#                                                     "(Away) precise passes",
#                                                     "(Away) modetate passes",
#                                                     "(Away) failed blocks"
#                                                   )

save(data_by_sets,file="data_by_sets")

#-----New variables concerning the evaluation of skills
data_by_sets_new_variables<-data_by_sets[,c(2:47,56:59)]/(data_by_sets$home_points+data_by_sets$away_points)
colnames(data_by_sets_new_variables)<-paste0(colnames(data_by_sets[,c(2:47,56:59)]),"_ratio_total_points")

data_by_sets<-cbind(data_by_sets,data_by_sets_new_variables)


#------X_home skills
colnames(data_by_sets)
X_home<-data_by_sets[c("Home_perfect_serves",
                       "Home_very_good_serves","Home_failed_serves","Home_poor_passes",
                       "Home_failed_passes","Home_precise_passes","Home_modetate_passes",
                       "Home_perfect_att1","Home_blocked_att1","Home_failed_att1",
                       "Home_perfect_att2","Home_blocked_att2","Home_failed_att2",
                       "Home_failed_blocks","Home_perfect_blocks","Home_failed_settings")]

X_away<-data_by_sets[c("Away_perfect_serves",
                       "Away_very_good_serves","Away_failed_serves","Away_poor_passes",
                       "Away_failed_passes","Away_precise_passes","Away_modetate_passes",
                       "Away_perfect_att1","Away_blocked_att1","Away_failed_att1",
                       "Away_perfect_att2","Away_blocked_att2","Away_failed_att2",
                       "Away_failed_blocks","Away_perfect_blocks","Away_failed_settings")]
save(X_home,file="X_home_skills")
save(X_away,file="X_away_skills")