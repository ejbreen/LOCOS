source('R/setup_home.R')

# combine_pops <- function(T_pop, C_pop, pct){
#   T_pop_s <- head(T_pop, nrow(T_pop)*pct)
#   C_pop_s <- head(C_pop, nrow(C_pop)*pct)
#   Tot_pop <- rbind(T_pop_s, C_pop_s)
#   return(Tot_pop)
# }
# 
# library('simstudy')
# 
# source("R/TKA population simulation parameters.R")
# # out <- pop_sim_paramaters()
# C_pop_def <- out$C_pop_def
# T_pop_def <- out$T_pop_def
# Levels <- out$Levels
# 
# C_pop <- genData(60000, C_pop_def)
# T_pop <- genData(2000, T_pop_def)
# 
# # C_pop <- subset.data.frame(C_pop, select = -c(B, POLY, HEAD, APP))
# # T_pop <- subset.data.frame(T_pop, select = -c(B, POLY, HEAD, APP))
# 
# Tot_pop <- combine_pops(T_pop, C_pop, .01)

Tot_pop <- read.csv('sample_set_10000.csv')

match_exact <- exactMatch(IMPLANT ~ POLY_XPLE, data = Tot_pop)

match_dist <- match_on(IMPLANT ~ AGE + BMI + strata(POLY_XPLE),
                       within = match_exact, data = Tot_pop)
match_dist2 <-  match_on(glm(IMPLANT ~ AGE + BMI, data = Tot_pop,family = binomial),
                         data = Tot_pop, 
                         within = exactMatch(IMPLANT ~ B + POLY + HEAD + APP, 
                                             data = Tot_pop)
                         )

fm1 <- fullmatch(match_dist, min.controls = 2, max.controls = 2)

matched_pop <- cbind(Tot_pop, group = fm1, matched = matched(fm1))

matched_pop[matched_pop$group==1.20]



