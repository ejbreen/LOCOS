# Written by: Evan Breen
#
# Last Updated: 4/18/18
#
#

library(optmatch)
library(survival)

Cox_test <- function(pop){
  # Tests to see if the implant of interest is significant in a cox proportional
  # hazard model
  #
  # Args:
  #   pop: A dataframe that contains an implant indicator, age, and sex variable
  #
  # Returns:
  #   the p-value associated with the implant of interest
  
  pop$followup_date <- 15
  pop$followup_rev <- pop$ttr < pop$followup_date
  
  cox <- coxph(Surv(followup_date,followup_rev) ~ implant + age + sex, data = pop)
  
  
  sum <- summary(cox1)
  return(sum$coefficients[1,5])
}

#complete
Basic_Cox <- function(pop){
  # Tests the significance of the implant of interest in a basic cox proportional
  # hazards model
  # 
  # Args:
  #   pop: A dataframe containing simulated MARCQI data
  #
  # Returns:
  #   the p-value associated with the implant of interest
  
  cox_p <- Cox_test(pop)
  
  return(cox_p)
}

#not started
DAG_Cox <- function(pop){
  # 
  #
  # Args:
  #   pop: A dataframe containing simulated MARCQI data
  #
  # Returns:
  #   the p-value associated with the implant of interest
  
  
  return(-1)
}

#not started
Greedy_Propensity_Cox <- function(pop){
  # 
  #
  # Args:
  #   pop: A dataframe containing simulated MARCQI data
  #
  # Returns:
  #   the p-value associated with the implant of interest
  
  return(-1)
}

#complete 
Optimal_Propensity_Cox <- function(pop){
  # Tests the significance of the implant of interest in a cox proportional 
  # hazards model after optimal subset matching on propensity score
  #
  # Args:
  #   pop: A dataframe containing simulated MARCQI data
  #
  # Returns:
  #   the p-value associated with the implant of interest
  
  match_dist <-  match_on(glm(implant ~ age + sex, data = Tot_pop, family = binomial),
                          data = pop)
  
  fm <- fullmatch(implant ~ age + sex, data = Tot_pop, min.controls = 5, max.controls = 5)
  
  matched_pop <- cbind(pop, group = fm, matched = matched(fm))
  matched_pop$group <-  factor(matched_pop$group, levels = c(levels(matched_pop$group), 0))
  
  matched_pop <- matched_pop[matched_pop$matched == TRUE,]
  
  cox_p <- Cox_test(matched_pop)
  
  return(cox_p)
}

#still needs work
Optimal_Covariate_Cox <- function(pop){
  # Tests the significance of the implant of interest in a cox proportional 
  # hazards model after optimal subset matching on propensity score with 
  # covariate restrictions
  #
  # Args:
  #   pop: A dataframe containing simulated MARCQI data
  #
  # Returns:
  #   the p-value associated with the implant of interest
  match_dist <-  match_on(glm(implant ~ age + sex, data = Tot_pop, family = binomial),
                          data = pop, 
                          within = exactMatch(implant ~ CoP, data = pop))
  
  fm <- fullmatch(implant ~ age + sex, data = Tot_pop, min.controls = 5, max.controls = 5)
  
  matched_pop <- cbind(pop, group = fm, matched = matched(fm))
  matched_pop$group <-  factor(matched_pop$group, levels = c(levels(matched_pop$group), 0))
  
  matched_pop <- matched_pop[matched_pop$matched == TRUE,]
  
  cox_p <- Cox_test(matched_pop)
  
  return(cox_p)
  
}

