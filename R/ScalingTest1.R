library(slam)
library(gurobi)

source('R/designmatch_adapted.R')

runScalingTest <- function(Tot_pop, scaling_factor){
  
  T_all = proc.time()
  T_build = T_all
  
  attach(Tot_pop)
  
  t_ind = IMPLANT
  
  x_mat = cbind(B_MOP, B_COP, B_COC, B_DM, POLY_UHWMPE, POLY_XPLE, POLY_A_XPLE,
                 HEAD_22mm, HEAD_28mm, HEAD_32mm, HEAD_36mm, HEAD_40mm, HEAD_44mm,
                 APP_anterior, APP_anterolateral, APP_posterior, APP_transtrochanteric,
                 S_VOLLUME, AGE, FEMALE, BMI)
  dist_mat = distmat(t_ind, x_mat)
  
  
  subset_weight = median(dist_mat)
  
  # Balance the covariate moments
  mom_covs = NULL
  mom_tols = NULL
  mom = list(covs = x_mat, tols = mom_tols)
  
  # # make the distributions look similar
  ks_covs = NULL
  ks_n_grid= 10
  ks_tols = NULL
  ks = list(covs = ks_covs, n_grid = ks_n_grid, tols = ks_tols)
  
  # exact within group matching
  exact_covs = cbind(B_MOP, B_COP, POLY_XPLE)
  exact = list(covs = exact_covs)
  
  # exact with an allowed difference
  near_exact_covs = NULL
  near_exact_devs = NULL
  near_exact = list(covs = near_exact_covs, devs = near_exact_devs)
  
  # exact overall matching
  fine_covs = cbind()
  find = list(covs = fine_covs)
  
  # fine but with an allowed differenec
  near_fine_covs = NULL
  near_fine_devs = NULL
  near_fine = list(covs = near_fine_covs, devs = near_fine_devs)
  
  t_max = 60*5
  solver = "gurobi" 
  approximate = 1
  solver = list(name = solver, t_max = t_max, approximate = approximate,
                round_cplex = 0, trace_cplex = 0)
  
  T_build = proc.time()-T_build
  T_run = proc.time()
  
  out = bmatch(t_ind = t_ind, dist_mat = dist_mat, exact = exact, n_controls = 5, 
               total_groups = sum(t_ind), solver = solver)
  
  T_run = proc.time()-T_run
  T_all = proc.time()-T_all
  
  return()
}