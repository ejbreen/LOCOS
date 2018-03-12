
source('R/setup_home.R')
# source('R/setup_flux.R')

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

Rprofmem(filename = 'Rprof.out', append = FALSE, threshold = 1000)
out <-  runScalingTest(T_pop, C_pop, .1)
Timings <- out$Timings
bmatch_out <- out$bmatch_out

