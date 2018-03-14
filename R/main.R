
# source('R/setup_home.R')
source('R/setup_flux.R')

library('simstudy')

source("R/TKA population simulation parameters.R")
out <- pop_sim_paramaters()
C_pop_def <- out$C_pop_def
T_pop_def <- out$T_pop_def
Levels <- out$Levels

C_pop <- genData(60000, C_pop_def)
T_pop <- genData(2000, T_pop_def)

C_pop <- subset.data.frame(C_pop, select = -c(B, POLY, HEAD, APP))
T_pop <- subset.data.frame(T_pop, select = -c(B, POLY, HEAD, APP))

source('R/ScalingTest1.R')

scales <-  c(.5, .55, .6)
scales_small <- c(.01, .05, .1, .15)

TimingDF <- data.frame(scale_factor = 0,
                      section = 0,
                      user = 0.0,
                      system = 0.0,
                      total = 0.0)

for(i in scales){
  print(paste('Running at scaling factor- ', i, sep = ""))
  out <-  runScalingTest(T_pop, C_pop, i)
  TimingDF <- rbind(TimingDF, c(i, 1, unname(out$Timings$setup[1:3])))
  TimingDF <- rbind(TimingDF, c(i, 2, unname(out$Timings$run[1:3])))
  TimingDF <- rbind(TimingDF, c(i, 3, unname(out$Timings$all[1:3])))
}

write.csv(TimingDF, file = paste('Timings-',scales[1],'-',tail(scales, 1), '.csv', sep = ""))

