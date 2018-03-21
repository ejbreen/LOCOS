library(slam)
library(gurobi)

source('R/designmatch_adapted.R')

combine_pops <- function(T_pop, C_pop, pct){
  T_pop_s <- head(T_pop, nrow(T_pop)*pct)
  C_pop_s <- head(C_pop, nrow(C_pop)*pct)
  Tot_pop <- rbind(T_pop_s, C_pop_s)
  return(Tot_pop)
}

runScalingTest <- function(T_pop, C_pop, scaling_factor){
  
  M_start = mem_used()
  
  p = combine_pops(T_pop, C_pop, scaling_factor)
  
  M_combined = mem_used()
  T_all = proc.time()
  T_build = T_all
  
  t_ind = p$IMPLANT
  
  x_mat = cbind(p$B_MOP, p$B_COP, p$B_COC, p$B_DM, p$POLY_UHWMPE, p$POLY_XPLE, p$POLY_A_XPLE,
                p$HEAD_22mm, p$HEAD_28mm, p$HEAD_32mm, p$HEAD_36mm, p$HEAD_40mm, p$HEAD_44mm,
                p$APP_anterior, p$APP_anterolateral, p$APP_posterior, p$APP_transtrochanteric,
                p$S_VOLLUME, p$AGE, p$FEMALE, p$BMI)
  dist_mat = distmat(t_ind, x_mat)
  
  M_distmat = mem_used()
  
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
  exact_covs = cbind(p$B_MOP, p$B_COP, p$POLY_XPLE)
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
  
  t_max = 60*60*2
  solver = "gurobi" 
  approximate = 1
  solver = list(name = solver, t_max = t_max, approximate = approximate,
                round_cplex = 0, trace_cplex = 0)
  
  M_setup = mem_used()
  T_build = proc.time()-T_build
  T_run = proc.time()
  
  out = bmatch(t_ind = t_ind, dist_mat = dist_mat, exact = exact, n_controls = 5, 
               total_groups = sum(t_ind), solver = solver)
  
  M_final = mem_used()
  T_run = proc.time()-T_run
  T_all = proc.time()-T_all
  Timings = list(all=T_all, setup=T_build, run=T_run)
  Mem = list(start = M_start, combined = M_combined, distmat = M_distmat, 
             setup = M_setup, final = M_final)
  Timings
  Mem
  return(list(Timings = Timings, Mem = Mem, bmatch_out = out))
}
