
source('R/setup_home.R')
# source('R/setup_chromebook.R')
# source('R/setup_flux.R')

source('R/match_analysis.R')

DF <- data.frame(population = character(),
                 basic = double(),
                 DAG = double(),
                 greedy = double(),
                 optimal_prop = double(),
                 optimal_cov = double(),
                 stringsAsFactors = FALSE)

set_pop <- function(complexity, effect, iteration){
  x = paste("complexity-",complexity,"_",
            "effect-",effect,"_",
            "iteration-",iteration,
            sep = "")
  return(x)
}

complexities  <- c(1)
effect_levels <- c(1)
iterations     <- c(0:10)

for (comp in complexities){
  for (effect in effect_levels){
    for (iter in iterations){
      pop_name <- set_pop(comp,effect,iter)
      pop <- read.csv(paste("Data/",pop_name,".csv",sep = ""))
      
      pop <- head(pop, 5000)
      
      p1 <- Basic_Cox(pop)
      p2 <- DAG_Cox(pop)
      p3 <- Greedy_Propensity_Cox(pop)
      p4 <- Optimal_Propensity_Cox(pop)
      p5 <- Optimal_Covariate_Cox(pop)
      
      DF <- rbind(DF, c(pop_name, p1, p2, p3, p4, p5))
    }
  }
}

write.csv(DF, file = paste("match_comparison",'.csv', sep = ""))

