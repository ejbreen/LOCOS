# Required packages:    Rglpk
#                       designmatch
#                       gurobi      'C:/gurobi752/win64/R/gurobi_7.5-2.zip'

install.packages('C:/gurobi752/win64/R/gurobi_7.5-2.zip', repos = NULL)

library(gurobi)
library(Rglpk)
library(designmatch)


combine_pops <- function(T_pop, C_pop, pct){
  T_pop_s <- head(T_pop, nrow(T_pop)*pct)
  C_pop_s <- head(C_pop, nrow(C_pop)*pct)
  Tot_pop <- rbind(T_pop_s, C_pop_s)
  return(Tot_pop)
}

Tot_pop <- combine_pops(T_pop, C_pop, .05)

attach(Tot_pop)

t_ind <- IMPLANT

x_mat <- cbind(B_MOP, B_COP, B_COC, B_DM, POLY_UHWMPE, POLY_XPLE, POLY_A_XPLE,
               HEAD_22mm, HEAD_28mm, HEAD_32mm, HEAD_36mm, HEAD_40mm, HEAD_44mm,
               APP_anterior, APP_anterolateral, APP_posterior, APP_transtrochanteric,
               S_VOLLUME, AGE, FEMALE, BMI)

dist_mat <- distmat(t_ind, x_mat)

subset_weight <- median(dist_mat)

mom_tols <- round(absstddif(x_mat, t_ind, 10), 2)
mom <- list(covs = x_mat, tols = mom_tols)


t_max = 60*5
solver = "gurobi" 
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate,
              round_cplex = 0, trace_cplex = 0)

out = bmatch(t_ind = t_ind, dist_mat = dist_mat, solver = solver)



t_id <- out$t_id
c_id <- out$c_id
out$group_id
meantab(x_mat, t_ind, t_id, c_id)


