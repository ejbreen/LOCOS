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

# Balance the covariate moments
mom_covs <- cbind()
mom_tols <- round(absstddif(x_mat, t_ind, 10), 2)
mom <- list(covs = x_mat, tols = mom_tols)

# # make the distributions look similar
ks_covs <- NULL
ks_n_grid <-NULL
ks_tols <- NULL
ks <- list(covs = ks_covs, n_grid = ks_n_grid, tols = ks_tols)

# exact within group matching
exact_covs <- cbind()
exact <- list(covs = exact_covs)

# exact with an allowed difference
near_exact_covs <- NULL
near_exact_devs <- NULL
near_exact <- list(covs = near_exact_covs, devs = near_exact_devs)

# exact overall matching
fine_covs <- cbind()
find <- list(covs = fine_covs)

# fine but with an allowed differenec
near_fine_covs <- NULL
near_fine_devs <- NULL
near_fine <- list(covs = near_fine_covs, devs = near_fine_devs)

t_max = 60*5
solver = "gurobi" 
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate,
              round_cplex = 0, trace_cplex = 0)

out = bmatch(t_ind = t_ind, dist_mat = dist_mat, n_controls = 5, total_groups = 100, 
             solver = solver)



t_id <- out$t_id
c_id <- out$c_id
out$group_id
meantab(x_mat, t_ind, t_id, c_id)


