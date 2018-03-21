
# source('R/setup_home.R')
# source('R/setup_chromebook.R')
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

scales <-  c(.1, .15, .2, .25, .3)
scales_small <- c(.01, .05, .1, .15)

DF <- data.frame(scale_factor = 0,
                 section = 'text',
                 user = 0.0,
                 system = 0.0,
                 total = 0.0,
                 mem = 0.0,
                 stringsAsFactors = FALSE)

for(i in scales){
  print(paste('Running at scaling factor- ', i, sep = ""))
  out <-  runScalingTest(T_pop, C_pop, i)
  DF <- rbind(DF, c(i, 'T_setup', unname(out$Timings$setup[1:3]), 0))
  DF <- rbind(DF, c(i, 'T_run', unname(out$Timings$run[1:3]), 0))
  DF <- rbind(DF, c(i, 'T_all', unname(out$Timings$all[1:3]), 0))
  DF <- rbind(DF, c(i, 'M_start', 0, 0, 0, out$Mem$start))
  DF <- rbind(DF, c(i, 'M_combined', 0, 0, 0, out$Mem$combined))
  DF <- rbind(DF, c(i, 'M_distmat', 0, 0, 0, out$Mem$distmat))
  DF <- rbind(DF, c(i, 'M_setup', 0, 0, 0, out$Mem$setup))
  DF <- rbind(DF, c(i, 'M_final', 0, 0, 0, out$Mem$final))
}

write.csv(DF, file = paste('Scaling-',scales[1],'-',tail(scales, 1), '.csv', sep = ""))

