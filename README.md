# Bayesian_Variable_Selection_Volleyball
In this repository, we have both R and Stan codes available as well as the proper data in order to implement Bayesian variable Selection (Gibbs method) for Volleyball Set Determination.


In each file, there are available the proper data, both R and Stan codes in order to implement the Gibbs Bayesian Variable Selection. 
In folders with suffix "_Skills" we refer to models which include only skill events as covariates (and mu, home effect for ZDTS model as well as the intercepts of ordered logistic)
while these ones with suffix "TA_Skills" (TA=Team Abilities) we refer to models including also as additional the home and away team (sum to zero sontraint) predictors additionally to other ones.
