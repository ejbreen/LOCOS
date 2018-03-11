
####################
# Matching Functions
####################

#! bmatch
bmatch = function(t_ind, dist_mat = NULL, subset_weight = NULL, n_controls = 1, total_groups = NULL,
                  mom = NULL,
                  ks = NULL,
                  exact = NULL,
                  near_exact = NULL,
                  fine = NULL,
                  near_fine = NULL, 
                  near = NULL,
                  far = NULL,
                  #use_controls = NULL,
                  solver = NULL) {
  
  use_controls = NULL
  
  if (is.null(mom)) {
    mom_covs = NULL
    mom_tols = NULL
    mom_targets = NULL
  } else {
    mom_covs = mom$covs
    mom_tols = mom$tols
    mom_targets = mom$targets
  }
  
  if (is.null(ks)) {
    ks_covs = NULL
    ks_n_grid = 10
    ks_tols = NULL
  } else {
    ks_covs = ks$covs
    ks_n_grid = ks$n_grid
    ks_tols = ks$tols
  }
  
  if (is.null(exact)) {
    exact_covs = NULL
  } else {
    exact_covs = exact$covs
  }
  
  if (is.null(near_exact)) {
    near_exact_covs = NULL
    near_exact_devs = NULL
  } else {
    near_exact_covs = near_exact$covs
    near_exact_devs = near_exact$devs
  }
  
  if (is.null(fine)) {
    fine_covs = NULL
  } else {
    fine_covs = fine$covs
  }
  
  if (is.null(near_fine)) {
    near_fine_covs = NULL
    near_fine_devs = NULL
  } else {
    near_fine_covs = near_fine$covs
    near_fine_devs = near_fine$devs
  }
  
  if (is.null(near)) {
    near_covs = NULL
    near_pairs = NULL
    near_groups = NULL
  } else {
    near_covs = near$covs
    near_pairs = near$pairs
    near_groups = near$groups
  }
  
  if (is.null(far)) {
    far_covs = NULL
    far_pairs = NULL
    far_groups = NULL
  } else {
    far_covs = far$covs
    far_pairs = far$pairs
    far_groups = far$groups
  }
  
  if (is.null(solver)) {
    solver = 'gurobi'
    t_max = 60 * 15
    approximate = 1
  } else {
    t_max = solver$t_max
    approximate = solver$approximate
    trace = solver$trace
    round_cplex = solver$round_cplex
    solver = 'gurobi'
  }
  
  
  
  #! CALL ERROR HANDLING
  
  if (is.null(subset_weight)) {
    subset_weight = 0
  }
  
  #! Generate the parameters
  cat(format("  Building the matching problem..."), "\n")
  prmtrs = .problemparameters(t_ind, dist_mat, subset_weight, n_controls, total_groups,
                              mom_covs, mom_tols, mom_targets,
                              ks_covs, ks_n_grid, ks_tols,
                              exact_covs,
                              near_exact_covs, near_exact_devs,
                              fine_covs,
                              near_fine_covs, near_fine_devs,
                              near_covs, near_pairs, near_groups,
                              far_covs, far_pairs, far_groups,
                              use_controls,
                              approximate)
  n_t = prmtrs$n_t
  n_c = prmtrs$n_c
  
  cvec = prmtrs$cvec
  Amat = prmtrs$Amat
  bvec = prmtrs$bvec
  ub = prmtrs$ub 
  sense = prmtrs$sense
  vtype = prmtrs$vtype
  c_index = prmtrs$c_index
  
  
  #! Find matches and calculate the elapsed time
  #! Gurobi
  #library(gurobi)
  if (requireNamespace('gurobi', quietly=TRUE)) {
    cat(format("  Gurobi optimizer is open..."), "\n")
    model = list()
    model$obj = cvec
    model$A = Amat
    model$sense = rep(NA, length(sense))
    model$sense[sense=="E"] = '='
    model$sense[sense=="L"] = '<='
    model$sense[sense=="G"] = '>='
    model$rhs = bvec
    model$vtypes = vtype
    model$ub = ub
    
    t_lim = list(TimeLimit = t_max, OutputFlag = trace)
    
    cat(format("  Finding the optimal matches..."), "\n")
    ptm = proc.time()
    out = gurobi::gurobi(model, t_lim)
    time = (proc.time()-ptm)[3]
    
    if (out$status == "INFEASIBLE") {
      cat(format("  Error: problem infeasible!"), "\n")
      obj_total = NA
      obj_dist_mat = NA
      t_id = NA
      c_id = NA
      group_id = NA
      time = NA
    }
    
    if (out$status ==  "OPTIMAL" || out$status == "TIME_LIMIT") {
      
      if (out$status == "OPTIMAL") {
        cat(format("  Optimal matches found"), "\n")
      }
      
      else {
        cat(format("  Time limit reached, best suboptimal solution given"), "\n")
      }
      
      if (approximate == 1) {
        rel = .relaxation_b(n_t, n_c, out$x, dist_mat, subset_weight, "gurobi", round_cplex, trace)
        out$x = rel$sol
        out$objval = rel$obj
        time = time + rel$time
      }
      
      #! Matched units indexes
      t_id = unique(sort(rep(1:n_t, n_c))[out$x[1:(n_t*n_c)]==1])
      c_id = (c_index+n_t)[out$x[1:(n_t*n_c)]==1]
      
      #! Group (or pair) identifier
      group_id_t = 1:(length(t_id))
      group_id_c = sort(rep(1:(length(t_id)), n_controls))
      group_id = c(group_id_t, group_id_c)
      
      #! Optimal value of the objective function
      obj_total = out$objval
      
      if (!is.null(dist_mat)) {
        obj_dist_mat = sum(c(as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE)) * out$x[1:(n_t*n_c)]==1))
      } else {
        obj_dist_mat = NULL
      }
    }
  } else {
    stop('Required solver not installed')
  }
  
  #! Output
  return(list(obj_total = obj_total, obj_dist_mat = obj_dist_mat, 
              t_id = t_id, c_id = c_id, group_id = group_id, time = time, status=out$status))
}


#! cardmatch
cardmatch = function(t_ind, mom = NULL, fine = NULL, solver = NULL) {
  
  if (is.null(mom)) {
    mom_covs = NULL
    mom_tols = NULL
    mom_targets = NULL
  } else {
    mom_covs = mom$covs
    mom_tols = mom$tols
    mom_targets = mom$targets
  }    
  if (is.null(fine)) {
    fine_covs = NULL
  } else {
    fine_covs = fine$covs
  }
  
  if (is.null(solver)) {
    solver = 'glpk'
    t_max = 60 * 15
    approximate = 1
  } else {
    t_max = solver$t_max
    approximate = solver$approximate
    trace = solver$trace
    round_cplex = solver$round_cplex
    solver = solver$name
  }
  
  #! CALL ERROR HANDLING
  
  #! Generate the parameters
  cat(format("  Building the matching problem..."), "\n")
  prmtrs = .problemparameters_cardmatch(t_ind, mom_covs, mom_tols, mom_targets, fine_covs)
  n_t = prmtrs$n_t
  n_c = prmtrs$n_c
  n_dec_vars = prmtrs$n_dec_vars
  cvec = prmtrs$cvec
  Amat = prmtrs$Amat
  bvec = prmtrs$bvec
  sense = prmtrs$sense
  vtype = prmtrs$vtype
  
  #! Find matches and calculate the elapsed time
  #! Gurobi
  if (solver == "gurobi") {
    #library(gurobi)
    if (requireNamespace('gurobi', quietly=TRUE)) {
      cat(format("  Gurobi optimizer is open..."), "\n")
      model = list()
      model$modelsense = 'max'
      model$obj = cvec
      model$A = Amat
      model$sense = rep(NA, length(sense))
      model$sense[sense=="E"] = '='
      model$sense[sense=="L"] = '<='
      model$sense[sense=="G"] = '>='
      model$rhs = bvec
      model$vtypes = vtype
      
      t_lim = list(TimeLimit = t_max, OutputFlag = trace)
      
      cat(format("  Finding the optimal matches..."), "\n")
      ptm = proc.time()
      out = gurobi::gurobi(model, t_lim)
      time = (proc.time()-ptm)[3]
      
      if (out$status == "INFEASIBLE") {
        cat(format("  Error: problem infeasible!"), "\n")
        obj_total = NA
        obj_dist_mat = NA
        t_id = NA
        c_id = NA
        group_id = NA
        time = NA
      }
      
      if (out$status ==  "OPTIMAL" || out$status == "TIME_LIMIT") {
        if (out$status == "OPTIMAL") {
          cat(format("  Optimal matches found"), "\n")
        }
        
        else {
          cat(format("  Time limit reached, best suboptimal solution given"), "\n")
        }
        
        #! Matched units indexes
        t_id = (1:n_dec_vars)[t_ind==1 & out$x==1]
        c_id = (1:n_dec_vars)[t_ind==0 & out$x==1]
        
        #! Group (or pair) identifier
        group_id = c(1:(length(t_id)), 1:(length(c_id)))
        
        #! Optimal value of the objective function
        obj_total = out$objval   
      }
    } else {
      stop('Required solver not installed')
    }
    
  }
  #! Output
  return(list(obj_total = obj_total, t_id = t_id, c_id = c_id, group_id = group_id, time = time))
}



#! distmatch
distmatch = function(t_ind, dist_mat = NULL, solver = NULL) {
  
  # Subset matching weight
  subset_weight = 0
  
  # Total number of matched pairs
  total_groups = sum(t_ind)
  
  # Match
  out = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight, total_groups = total_groups, solver = solver)
  
  #! Output
  return(list(obj_total = out$obj_total, obj_dist_mat = out$obj_dist_mat, 
              t_id = out$t_id, c_id = out$c_id, group_id = out$group_id, time = out$time, status = out$status))
  
}


########################
# Input helper functions
########################

# Absolute Standard Difference
absstddif <-
  function(X_mat, t_ind, std_dif) {
    n_vrs = ncol(X_mat)
    n_obs = nrow(X_mat)
    mom_tols_out = NA
    for (j in 1:n_vrs) {
      yes_before = unlist(X_mat[t_ind == 1, j])
      no_before = unlist(X_mat[t_ind == 0, j])
      pooled_sd = sqrt((var(yes_before, na.rm = TRUE) + var(no_before, na.rm = TRUE))/2)
      mom_tols_out = c(mom_tols_out, pooled_sd*std_dif)
    }
    mom_tols_out = mom_tols_out[-1]
    mom_tols_out
  }

# distance matrix
distmat <-
  function (t_ind, X_mat, calip_cov = NULL, calip_size = NULL, calip_penalty = NULL, near_exact_covs = NULL, near_exact_penalties = NULL, digits = 1) {
    dist_mat = .smahal(t_ind, X_mat)
    dist_mat = dist_mat/mean(dist_mat)
    if (!is.null(calip_cov)) {
      dist_mat = .addcalip(dist_mat, t_ind, calip_cov, calip_size, calip_penalty)
    }
    if (!is.null(near_exact_covs)) {
      for (i in 1:ncol(near_exact_covs)) {
        penalty_mat = abs(outer(near_exact_covs[t_ind==1, i], near_exact_covs[t_ind==0, i], "-"))
        penalty_mat = (penalty_mat!=0)*(near_exact_penalties[i])
        dist_mat = dist_mat+penalty_mat
      }
    }
    dist_mat = round(dist_mat/mean(dist_mat), digits)
    dist_mat
  }


###########################
# Result analysis functions
###########################

ecdfplot <-
  function(x, t_id, c_id, main_title = "", legend_position = "right") {
    xlim_min = quantile(x, .01)
    xlim_max = quantile(x, .99)
    aux1 = ecdf(x[t_id])
    aux2 = ecdf(x[c_id])
    if (length(c_id) > 250) {
      aux2 = ecdf(sample(x, 250))
    }
    plot(aux1, xlim = c(xlim_min, xlim_max), ylim = c(0, 1), xlab = "x", ylab = "CDF(x)", main = main_title, col = "black", pch = 0)
    par(new=TRUE)
    plot(aux2, xlim = c(xlim_min, xlim_max), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "", main = "", col = "blue", pch = 8)
    legend_aux1 = c("Treated", "Controls")
    legend_aux2 = c("black", "blue")
    legend_aux3 = c(0, 8)
    legend(legend_position, legend_aux1, bty = "n", pch = legend_aux3, col = legend_aux2)
    par(new=FALSE)	
  }


finetab <-
  function(nom_cov, t_id, c_id) {
    Cat = nom_cov
    Units = rep(0, length(Cat))
    Units[t_id] = 1
    Units[c_id] = 2
    tab = table(Cat, Units)[, 2:3]
    colnames(tab) = c("T", "C")
    tab
  }


loveplot <-
  function (X_mat, t_id, c_id, v_line, legend_position = "topright") {
    X_mat_t = X_mat[t_id, ]
    X_mat_c_before = X_mat[-t_id, ]
    X_mat_c_before_mean = apply(X_mat_c_before, 2, mean)
    X_mat_t_mean = apply(X_mat_t, 2, mean)
    X_mat_t_var = apply(X_mat_t, 2, var)
    X_mat_c_before_var = apply(X_mat_c_before, 2, var)
    std_dif_before = (X_mat_t_mean - X_mat_c_before_mean)/sqrt((X_mat_t_var + X_mat_c_before_var)/2)
    X_mat_c_after = X_mat[c_id, ]
    X_mat_c_after_mean = apply(X_mat_c_after, 2, mean)
    std_dif_after = (X_mat_t_mean - X_mat_c_after_mean)/sqrt((X_mat_t_var + X_mat_c_before_var)/2)
    #library("lattice")
    abs_std_dif_before = abs(std_dif_before)
    n_aux = length(abs_std_dif_before)
    abs_std_dif_after = abs(std_dif_after)
    dotchart(abs_std_dif_before[n_aux:1], labels = colnames(X_mat)[n_aux:1], cex = 0.7, pch = "", color = , main = "", xlim = c(0, 1), xlab = "Absolute standardized differences in means")
    points(abs_std_dif_before[n_aux:1], y = 1:ncol(X_mat), cex = 0.9, pch = 0)
    points(abs_std_dif_after[n_aux:1], y = 1:ncol(X_mat), cex = 0.8, pch = 8, col = "blue")
    legend(legend_position, c("Before matching", "After matching"), cex = 0.8, bty = "n", pch = c(0, 8), col = c("black", "blue"))
    abline(v = v_line, lty = 2)
  }


meantab <-
  function(X_mat, t_ind, t_id, c_id, exact = NULL, digits = 2) {
    #t_ind = c(rep(1, length(t_id)), rep(0, length(c_id)))
    #X_mat = X_mat[c(t_id, c_id), ]
    n_vrs = ncol(X_mat)
    n_obs = nrow(X_mat)
    output = matrix(NA, n_vrs, 7)
    rownames(output) = colnames(X_mat)
    colnames(output) = c("Mis", "Min", "Max", "Mean T", "Mean C", "Std Dif", "P-val")
    for (j in 1:n_vrs) {
      n_miss = sum(is.na(X_mat[, j]))
      miss_pct = round(n_miss/n_obs, 2)
      if (n_miss < n_obs) {
        #yes = unlist(X_mat[t_ind==1, j])
        #no = unlist(X_mat[t_ind==0, j])
        #mean_yes = mean(yes, na.rm = TRUE)
        #mean_no = mean(no, na.rm = TRUE)
        #min_val = min(c(yes, no), na.rm = TRUE)
        #max_val = max(c(yes, no), na.rm = TRUE)
        #pooled_sd = sqrt((var(yes, na.rm = TRUE)+var(no, na.rm = TRUE))/2)
        yes = unlist(X_mat[t_id, j])
        no = unlist(X_mat[c_id, j])
        
        yes_before = unlist(X_mat[t_ind==1, j])
        no_before = unlist(X_mat[t_ind==0, j])
        
        mean_yes = mean(yes, na.rm = TRUE)
        mean_no = mean(no, na.rm = TRUE)
        min_val = min(c(yes, no), na.rm = TRUE)
        max_val = max(c(yes, no), na.rm = TRUE)
        pooled_sd = sqrt((var(yes_before, na.rm = TRUE)+var(no_before, na.rm = TRUE))/2)
        std_diff = NA
        p_val = NA
        if (max_val > min_val) {
          std_diff = (mean_yes-mean_no)/pooled_sd
          if (is.null(exact)) {
            p_val = t.test(yes, no)$p.value
          }
          if (!is.null(exact)) {
            if (exact[j]=="w") {
              p_val = wilcox.test(yes, no)$p.value
            }
            if (exact[j]=="f") {
              tab_aux = table(t_ind[c(t_id, c_id)], unlist(X_mat[c(t_id, c_id), j]) )
              if (ncol(tab_aux)==1) {
                p_val = 1
              }
              if (ncol(tab_aux)!=1) {
                p_val = fisher.test(tab_aux)$p.value
              }
            }	
          }
        }
        output[j, ] = c(miss_pct, min_val, max_val, mean_yes, mean_no, std_diff, p_val)
      }
    }
    round(output, digits)
  }


pairsplot <-
  function(cov1, cov2, t_id, c_id, xlab = "", ylab = "", main = "") {
    plot(cov1[t_id], cov2[t_id], xlim = c(min(cov1[c(t_id, c_id)]), max(cov1[c(t_id, c_id)])), ylim = c(min(cov2[c(t_id, c_id)]), max(cov2[c(t_id, c_id)])), xlab = xlab, ylab = ylab, main = main, pch = 1, col = "red")
    points(cov1[c_id], cov2[c_id], pch = 0, col = "blue")
    for (i in 1:length(t_id)) {
      segments(cov1[t_id[i]], cov2[t_id[i]], cov1[c_id[i]], 
               cov2[c_id[i]], col = "grey")
    }
    abline(v = mean(cov1[t_id]), b = 0, col = "red")
    abline(v = mean(cov1[c_id]), b = 0, col = "blue")
    abline(h = mean(cov2[t_id]), b = 0, col = "red")
    abline(h = mean(cov2[c_id]), b = 0, col = "blue")
    legend_aux1 = c("Treated", "Controls")
    legend_aux2 = c("red", "blue")
    legend_aux3 = c(1, 0)
    legend("topright", legend_aux1, bty = "n", pch = legend_aux3, col = legend_aux2)
  }


####################
# Internal functions
####################

#! From Paul's Design of Obs. Studies (p. 251): "...adds a penalty to a distance matrix dmat for violations of the calip_size." "The default calip_size width is 0.2*sd(p). The magnitude of the penalty is penalty multiplied by the magnitude of the violation, where penalty is set to 1000 by default."
.addcalip = function (dist_mat, t_ind, calip_cov, calip_size, calip_penalty) {
  sdp = sd(calip_cov)
  penalty_mat = abs(outer(calip_cov[t_ind == 1], calip_cov[t_ind == 0], "-"))
  penalty_mat = (penalty_mat>(calip_size*sdp))*calip_penalty
  dist_mat = dist_mat+penalty_mat
  dist_mat
}



#! Builds the constraint matrix
.constraintmatrix = function(t_ind, n_controls, total_groups,
                             mom_covs, mom_tols, mom_targets,
                             ks_covs, ks_covs_aux, ks_n_grid, ks_tols,
                             exact_covs,
                             near_exact_covs, near_exact_devs,
                             fine_covs,
                             near_fine_covs, near_fine_devs,
                             near_covs, near_pairs, near_groups,
                             far_covs, far_pairs, far_groups,
                             use_controls,
                             approximate) {
  
  #! Number of treated units, number of controls
  n_t = sum(t_ind)
  n_c = length(t_ind)-n_t	
  
  #! Total number of units
  n_tot = n_t*n_c
  
  #! Build parts of the constraint matrix
  #! Part 1
  if (approximate == 1 | n_controls == 1) {
    row_ind_1 = sort(rep(1:n_t, n_c))
    col_ind_1 = 1:n_tot
    ones_1 = rep(1, n_tot)
  }
  else {
    row_ind_1 = c(sort(rep(1:n_t, n_c)), 1:n_t)
    col_ind_1 = 1:(n_tot+n_t)
    ones_1 = c(rep(1, n_tot), rep(-1*n_controls, n_t))
  }
  
  #! Part 2
  row_ind_2 = sort(rep(1:n_c, n_t))+n_t
  col_ind_2 = rep(seq(1, n_t*n_c, n_c), n_c)+(sort(rep(1:n_c, n_t))-1)
  ones_2 = rep(1, n_tot)
  #! Current max row index
  row_ind_cur	= max(row_ind_2)
  
  #! Parts 3 and 4: moments and K-S
  mom_ks_covs = NULL
  if ((!is.null(mom_covs) & is.null(mom_targets)) | !is.null(ks_covs)) {
    row_ind_3.4 = 0
    #! Number of moment covariates
    n_mom_covs = 0
    if(!is.null(mom_covs) & is.null(mom_targets)) {
      n_mom_covs = ncol(mom_covs)
    }
    #! Number of K-S covariates
    n_ks_covs = 0
    if(!is.null(ks_covs)) {
      n_ks_covs = ncol(ks_covs)
    }
    # Bind moment and K-S covariates		
    if(!is.null(mom_covs) & is.null(mom_targets) & is.null(ks_covs_aux)) {
      mom_ks_covs = mom_covs
      mom_ks_tols = mom_tols
    }
    if((is.null(mom_covs) & is.null(mom_targets)) & !is.null(ks_covs_aux)) {
      mom_ks_covs = ks_covs_aux
      mom_ks_tols = NA
      for (i in 1:ncol(ks_covs)) {
        mom_ks_tols = c(mom_ks_tols, rep(ks_tols[i], ks_n_grid[i]), rep(0, max(ks_n_grid)-ks_n_grid[i]))
      }
      mom_ks_tols = mom_ks_tols[-1]
    }
    if((!is.null(mom_covs) & is.null(mom_targets)) & !is.null(ks_covs_aux)) {
      mom_ks_covs = cbind(mom_covs, ks_covs_aux)
      mom_ks_tols = mom_tols
      for (i in 1:ncol(ks_covs)) {
        mom_ks_tols = c(mom_ks_tols, rep(ks_tols[i], ks_n_grid[i]), rep(0, max(ks_n_grid)-ks_n_grid[i]))
      }
    }	
  }					 
  if (!is.null(mom_ks_covs)) {
    n_mom_ks_covs = ncol(mom_ks_covs)
    if ((!is.null(mom_tols) & is.null(mom_targets)) | !is.null(ks_tols)) {
      row_ind_3.4 = sort(rep(1:(2*n_mom_ks_covs)+n_t+n_c, n_tot))
    }	
    col_ind_3.4 = NA
    mom_ks_vals_3.4 = NA
    j = 1
    k = 0
    for (i in 1:n_mom_ks_covs) {
      if (n_mom_covs != 0 & i <= n_mom_covs) {
        if ((!is.null(mom_tols) & is.null(mom_targets)) | !is.null(ks_tols)) {
          col_ind_3.4 = c(col_ind_3.4, rep(1:n_tot, 2))
        }
      }
      if (n_ks_covs != 0 & i > n_mom_covs) {
        if ((!is.null(mom_tols) & is.null(mom_targets)) | !is.null(ks_tols)) {
          col_ind_3.4 = c(col_ind_3.4, rep(1:n_tot, 2))
          k = k+1
          if (k >= max(ks_n_grid)) {
            j = j+1
            k = 0	
          }
        }
      }
      temp_mean_1 = rep(mom_ks_covs[t_ind==0, i], n_t)-(mom_ks_covs[t_ind==1, i])[sort(rep(1:n_t, n_c))]
      if ((!is.null(mom_tols) & is.null(mom_targets)) | !is.null(ks_tols)) {
        temp_mean_2 = temp_mean_1-(mom_ks_tols[i]*rep(1, n_t*n_c))
        temp_mean_3 = -temp_mean_1-(mom_ks_tols[i]*rep(1, n_t*n_c))
      }	
      mom_ks_vals_3.4 = c(mom_ks_vals_3.4, temp_mean_2, temp_mean_3)
      if (i == 1) {
        col_ind_3.4 = col_ind_3.4[-1]
        mom_ks_vals_3.4 = mom_ks_vals_3.4[-1]
      }
    }
    #! Current max row index
    row_ind_cur	= max(row_ind_3.4)
  }
  
  #! Moment target part
  rows_target = NULL
  cols_target = NULL
  vals_target = NULL
  if (!is.null(mom_covs) & !is.null(mom_targets)) {
    n_mom_covs = ncol(mom_covs)
    rows_target = sort(rep(1:(4*n_mom_covs)+row_ind_cur, n_tot))
    
    for (i in 1:n_mom_covs) {
      cols_target = c(cols_target, rep(1:n_tot, 4))
      temp_treatment_1 = (mom_covs[t_ind==1, i])[sort(rep(1:n_t, n_c))] - (mom_targets[i] + mom_tols[i])
      temp_treatment_2 = -1*(mom_covs[t_ind==1, i])[sort(rep(1:n_t, n_c))] + (mom_targets[i] - mom_tols[i])
      temp_control_1 = rep(mom_covs[t_ind==0, i], n_t) - (mom_targets[i] + mom_tols[i])
      temp_control_2 = -1*rep(mom_covs[t_ind==0, i], n_t) + (mom_targets[i] - mom_tols[i])
      vals_target = c(vals_target, temp_treatment_1, temp_treatment_2, temp_control_1, temp_control_2)
    }
    row_ind_cur = max(rows_target)
  }
  
  
  #! Part 5: exact
  rows_exact = NULL
  cols_exact = NULL
  vals_exact = NULL
  if (!is.null(exact_covs)) {
    n_exact_cats = ncol(exact_covs)
    for (i in 1:n_exact_cats) {
      rows_exact = c(rows_exact, rep(row_ind_cur+i, n_t*n_c))
      cols_exact = c(cols_exact, 1:(n_t*n_c))
      dist_exact_cov = abs(outer(exact_covs[t_ind==1, i], exact_covs[t_ind==0, i], "-"))
      dist_exact_cov = t(dist_exact_cov)
      vals_exact = c(vals_exact, as.vector(dist_exact_cov))
    }	
    row_ind_5 = rows_exact
    col_ind_5 = cols_exact
    exact_vals_5 = vals_exact
    row_ind_cur	= max(row_ind_5)
  }
  
  #! Part 6: near-exact
  rows_near_exact = NULL
  cols_near_exact = NULL
  vals_near_exact = NULL
  if (!is.null(near_exact_covs)) {
    n_near_exact_cats = ncol(near_exact_covs)
    for (i in 1:n_near_exact_cats) {
      rows_near_exact = c(rows_near_exact, rep(row_ind_cur+i, n_t*n_c))
      cols_near_exact = c(cols_near_exact, 1:(n_t*n_c))
      dist_near_exact_cov = abs(outer(near_exact_covs[t_ind==1, i], near_exact_covs[t_ind==0, i], "-"))
      dist_near_exact_cov = t(dist_near_exact_cov)
      vals_near_exact = c(vals_near_exact, as.vector(dist_near_exact_cov))
    }	
    row_ind_6 = rows_near_exact
    col_ind_6 = cols_near_exact
    near_exact_vals_6 = vals_near_exact
    row_ind_cur	= max(row_ind_6)
  }
  
  #! Part 7: fine
  bvec_7 = NULL
  rows_fine = NULL
  cols_fine = NULL
  vals_fine = NULL
  if (!is.null(fine_covs)) {
    #! Transform fine_covs to a matrix of binary inds. for each cat. of each fine balancing covariate
    fine_covs_2 = rep(NA, nrow(fine_covs))
    n_fine_covs = ncol(fine_covs)
    j = 1
    for (i in 1:n_fine_covs) {	
      aux = factor(fine_covs[, i])
      fine_covs_2 = cbind(fine_covs_2, diag(nlevels(aux))[aux,])
      if (j == 1) {
        fine_covs_2 = fine_covs_2[, -1]
      }
      j = j+1
    }
    n_fine_cats = ncol(fine_covs_2)
    for (i in 1:n_fine_cats) {
      rows_fine = c(rows_fine, rep(row_ind_cur+i, n_t*n_c))
      cols_fine = c(cols_fine, 1:(n_t*n_c))
      dist_fine_cov = outer(fine_covs_2[t_ind==1, i], fine_covs_2[t_ind==0, i], "-")
      dist_fine_cov = t(dist_fine_cov)
      vals_fine = c(vals_fine, as.vector(dist_fine_cov))
    }	
    row_ind_7 = rows_fine
    col_ind_7 = cols_fine
    fine_vals_7 = vals_fine
    bvec_7 = rep(0, n_fine_cats)
    row_ind_cur	= max(row_ind_7)
  }
  
  #! Part 8: near-fine
  bvec_8 = NULL
  rows_near_fine = NULL
  cols_near_fine = NULL
  vals_near_fine = NULL
  if (!is.null(near_fine_covs)) {
    #! Transform fine_covs to a matrix of binary inds. for each cat. of each fine balancing covariate
    near_fine_covs_2 = rep(NA, nrow(near_fine_covs))
    n_near_fine_covs = ncol(near_fine_covs)
    j = 1
    for (i in 1:n_near_fine_covs) {  
      near_aux = factor(near_fine_covs[, i])
      near_fine_covs_2 = cbind(near_fine_covs_2, diag(nlevels(near_aux))[near_aux,])
      if (j == 1) {
        near_fine_covs_2 = near_fine_covs_2[, -1]
      }
      j = j+1
    }
    n_near_fine_cats = ncol(near_fine_covs_2)
    for (i in 1:n_near_fine_cats) {
      rows_near_fine = c(rows_near_fine, rep(row_ind_cur+i, n_t*n_c))
      cols_near_fine = c(cols_near_fine, 1:(n_t*n_c))
      dist_near_fine_cov = outer(near_fine_covs_2[t_ind==1, i], near_fine_covs_2[t_ind==0, i], "-")
      dist_near_fine_cov = t(dist_near_fine_cov)
      vals_near_fine = c(vals_near_fine, as.vector(dist_near_fine_cov))
    }
    row_ind_cur = max(rows_near_fine)
    for (i in 1:n_near_fine_cats) {
      rows_near_fine = c(rows_near_fine, rep(row_ind_cur+i, n_t*n_c))
    }
    row_ind_8 = rows_near_fine
    col_ind_8 = c(cols_near_fine, cols_near_fine)
    near_fine_vals_8 = c(vals_near_fine, vals_near_fine)
    bvec_8 = rep(0, n_near_fine_cats)
    row_ind_cur	= max(row_ind_8)
  }
  
  #! Part 9: Far
  rows_ind_far_pairs = list()
  if (!is.null(far_covs)) {
    row_ind_9 = NULL
    col_ind_9 = NULL
    far_cov_vals_9 = NULL
    n_far_covs = ncol(far_covs)
    for (j in 1:n_far_covs) {
      far_cov = far_covs[,j]
      #! Far on average constraints
      if (!is.null(far_groups)) {
        far_group = far_groups[j]
        row_ind_far_all = sort(c(rep(row_ind_cur+1, n_tot)))		
        col_ind_far_all = rep(1:n_tot, 1)			
        temp_mean_3 = (-rep(far_cov[t_ind==0], n_t)+((far_cov[t_ind==1])[sort(rep(1:n_t, n_c))]))-(far_group*rep(1, n_t*n_c))
        vals_far_all = c(temp_mean_3)		
        row_ind_cur	= max(row_ind_far_all)
      }
      #! Far on all pairs constraints
      if (!is.null(far_pairs)) {
        far_pair = far_pairs[j]
        aux = abs(outer(far_cov[t_ind==1], far_cov[t_ind==0], FUN = "-"))
        temp = as.vector(matrix(t(aux), nrow = 1, byrow = TRUE))
        cols_ind_far_pairs = which(temp<far_pair)
        if (length(cols_ind_far_pairs)>0) {
          rows_ind_far_pairs[[j]] = row_ind_cur+(1:length(cols_ind_far_pairs))
          vals_far_pairs = rep(1, length(cols_ind_far_pairs))
          row_ind_cur	= max(rows_ind_far_pairs[[j]])
        }
        if (length(cols_ind_far_pairs)==0) {
          cols_ind_far_pairs = NULL
          rows_ind_far_pairs[[j]] = -1
          vals_far_pairs = NULL
        }
      }		
      #! Put together
      if (!is.null(far_groups) && is.null(far_pairs)) {
        row_ind_9 = c(row_ind_9, row_ind_far_all)
        col_ind_9 = c(col_ind_9, col_ind_far_all)
        far_cov_vals_9 = c(far_cov_vals_9, vals_far_all)
      }
      if (is.null(far_groups) && !is.null(far_pairs) && rows_ind_far_pairs[[j]] != -1) {
        row_ind_9 = c(row_ind_9, rows_ind_far_pairs[[j]])
        col_ind_9 = c(col_ind_9, cols_ind_far_pairs)
        far_cov_vals_9 = c(far_cov_vals_9, vals_far_pairs)
      }
      if (!is.null(far_groups) && !is.null(far_pairs) && rows_ind_far_pairs[[j]] != -1) {
        row_ind_9 = c(row_ind_9, row_ind_far_all, rows_ind_far_pairs[[j]])
        col_ind_9 = c(col_ind_9, col_ind_far_all, cols_ind_far_pairs)
        far_cov_vals_9 = c(far_cov_vals_9, vals_far_all, vals_far_pairs)
      }
      if (!is.null(far_groups) && !is.null(far_pairs) && rows_ind_far_pairs[[j]] == -1) {
        row_ind_9 = c(row_ind_9, row_ind_far_all)
        col_ind_9 = c(col_ind_9, col_ind_far_all)
        far_cov_vals_9 = c(far_cov_vals_9, vals_far_all)
      }
    }
  }
  
  #! Part 10: Near
  rows_ind_near_pairs = list()
  if (!is.null(near_covs)) {
    row_ind_10 = NULL
    col_ind_10 = NULL
    near_cov_vals_10 = NULL
    n_near_covs = ncol(near_covs)
    for (j in 1:n_near_covs) {
      near_cov = near_covs[,j]
      #! Near on average constraints
      if (!is.null(near_groups)) {
        near_group = near_groups[j]
        row_ind_near_all = sort(c(rep(row_ind_cur+1, n_tot)))  	
        col_ind_near_all = rep(1:n_tot, 1)			
        temp_mean_4 = (-rep(near_cov[t_ind==0], n_t)+((near_cov[t_ind==1])[sort(rep(1:n_t, n_c))]))-(near_group*rep(1, n_t*n_c))
        vals_near_all = c(temp_mean_4)		
        row_ind_cur	= max(row_ind_near_all)
      }
      #! Near on all pairs constraints
      if (!is.null(near_pairs)) {
        near_pair = near_pairs[j]
        aux = abs(outer(near_cov[t_ind==1], near_cov[t_ind==0], FUN = "-"))
        temp = as.vector(matrix(t(aux), nrow = 1, byrow = TRUE))
        cols_ind_near_pairs = which(temp>near_pair)
        if (length(cols_ind_near_pairs)>0) {
          rows_ind_near_pairs[[j]] = row_ind_cur+(1:length(cols_ind_near_pairs))
          vals_near_pairs = rep(1, length(cols_ind_near_pairs))
          row_ind_cur	= max(rows_ind_near_pairs[[j]])
        }
        if (length(cols_ind_near_pairs)==0) {
          cols_ind_near_pairs = NULL
          rows_ind_near_pairs[[j]] = -1
          vals_near_pairs = NULL
        }
      }		
      #! Put together
      if (!is.null(near_groups) && is.null(near_pairs)) {
        row_ind_10 = c(row_ind_10, row_ind_near_all)
        col_ind_10 = c(col_ind_10, col_ind_near_all)
        near_cov_vals_10 = c(near_cov_vals_10, vals_near_all)
      }
      if (is.null(near_groups) && !is.null(near_pairs) && rows_ind_near_pairs[[j]] != -1) {
        row_ind_10 = c(row_ind_10, rows_ind_near_pairs[[j]])
        col_ind_10 = c(col_ind_10, cols_ind_near_pairs)
        near_cov_vals_10 = c(near_cov_vals_10, vals_near_pairs)
      }
      if (!is.null(near_groups) && !is.null(near_pairs) && rows_ind_near_pairs[[j]] != -1) {
        row_ind_10 = c(row_ind_10, row_ind_near_all, rows_ind_near_pairs[[j]])
        col_ind_10 = c(col_ind_10, col_ind_near_all, cols_ind_near_pairs)
        near_cov_vals_10 = c(near_cov_vals_10, vals_near_all, vals_near_pairs)
      }
      if (!is.null(near_groups) && !is.null(near_pairs) && rows_ind_near_pairs[[j]] == -1) {
        row_ind_10 = c(row_ind_10, row_ind_near_all)
        col_ind_10 = c(col_ind_10, col_ind_near_all)
        near_cov_vals_10 = c(near_cov_vals_10, vals_near_all)
      }
    }
  }
  
  # Part 11: use controls
  if (!is.null(use_controls)) {
    use_controls = use_controls[(n_t+1):(n_t+n_c)]
    use_controls_aux = rep(use_controls, n_t)
    col_ind_11 = (1:n_tot)[use_controls_aux==1]
    row_ind_11 = rep(row_ind_cur+1, length(col_ind_11))
    use_controls_vals_11 = rep(1, length(col_ind_11))
    
    row_ind_cur = max(row_ind_11)
  }
  
  # Part 12: total_groups
  if (!is.null(total_groups)) {
    row_ind_12 = rep(row_ind_cur+1, n_t*n_c)
    col_ind_12 = 1:(n_t*n_c)
    ones_12 = rep(1, n_t*n_c)
    
    row_ind_cur = max(row_ind_12)
  }
  
  #! Put all the parts of the constraint matrix together
  #! Parts 1 and 2
  row_ind = c(row_ind_1, row_ind_2)
  col_ind = c(col_ind_1, col_ind_2)
  vals = c(ones_1, ones_2)
  #! Parts 3 and 4
  if (!is.null(mom_ks_covs)) {
    row_ind = c(row_ind, row_ind_3.4)
    col_ind = c(col_ind, col_ind_3.4)
    vals = c(vals, mom_ks_vals_3.4)
  }
  #! Part 3b
  if (!is.null(mom_covs) & !is.null(mom_targets)) {
    row_ind = c(row_ind, rows_target)
    col_ind = c(col_ind, cols_target)
    vals = c(vals, vals_target)
  }
  #! Part 5
  if (!is.null(exact_covs)) {
    row_ind = c(row_ind, row_ind_5)
    col_ind = c(col_ind, col_ind_5)
    vals = c(vals, exact_vals_5)
  }
  #! Part 6
  if (!is.null(near_exact_covs)) {
    row_ind = c(row_ind, row_ind_6)
    col_ind = c(col_ind, col_ind_6)
    vals = c(vals, near_exact_vals_6)
  }
  #! Part 7
  if (!is.null(fine_covs)) {
    row_ind = c(row_ind, row_ind_7)
    col_ind = c(col_ind, col_ind_7)
    vals = c(vals, fine_vals_7)
  }
  #! Part 8
  if (!is.null(near_fine_covs)) {
    row_ind = c(row_ind, row_ind_8)
    col_ind = c(col_ind, col_ind_8)
    vals = c(vals, near_fine_vals_8)
  }
  #! Part 9
  if (!is.null(far_covs)) {
    row_ind = c(row_ind, row_ind_9)
    col_ind = c(col_ind, col_ind_9)
    vals = c(vals, far_cov_vals_9)	
  }	
  #! Part 10
  if (!is.null(near_covs)) {
    row_ind = c(row_ind, row_ind_10)
    col_ind = c(col_ind, col_ind_10)
    vals = c(vals, near_cov_vals_10)  
  }	
  #! Part 11				 
  if (!is.null(use_controls)) {
    row_ind = c(row_ind, row_ind_11)
    col_ind = c(col_ind, col_ind_11)
    vals = c(vals, use_controls_vals_11)	
  }
  #! Part 12
  if (!is.null(total_groups)) {
    row_ind = c(row_ind, row_ind_12)
    col_ind = c(col_ind, col_ind_12)
    vals = c(vals, ones_12)
  }
  
  aux = cbind(row_ind, col_ind, vals)[order(col_ind), ]
  cnstrn_mat = simple_triplet_matrix(i = aux[, 1], j = aux[, 2], v = aux[, 3])
  
  #! Output  
  return(list(cnstrn_mat = cnstrn_mat, bvec_7 = bvec_7, bvec_8 = bvec_8, rows_ind_far_pairs = rows_ind_far_pairs, rows_ind_near_pairs = rows_ind_near_pairs))
  
}



#! Generate the parameters for bpmatch
.problemparameters = function(t_ind, dist_mat, subset_weight, n_controls, total_groups,
                              mom_covs, mom_tols, mom_targets,
                              ks_covs, ks_n_grid, ks_tols,
                              exact_covs,
                              near_exact_covs, near_exact_devs,
                              fine_covs,
                              near_fine_covs, near_fine_devs,
                              near_covs, near_pairs, near_groups,
                              far_covs, far_pairs, far_groups,
                              use_controls,
                              approximate) {
  
  #! Number of treated units and controls
  n_t = sum(t_ind)
  n_c = length(t_ind)-n_t
  
  #! Number of dec. vars.
  n_dec_vars = n_t*n_c
  
  #! Number of moment covariates
  n_mom_covs = 0
  if(!is.null(mom_covs)) {
    n_mom_covs = ncol(mom_covs)
  }
  
  #! Number of K-S covariates
  n_ks_covs = 0
  if(!is.null(ks_covs)) {
    n_ks_covs = ncol(ks_covs)
    
    if ((length(ks_n_grid)==1) && (n_ks_covs > 1)) {
      ks_n_grid = rep(ks_n_grid, n_ks_covs)
    }
    
  }
  #! Parameters used to minimize the K-S statistic
  ks_covs_aux = NULL
  if (is.null(ks_covs)) {
    max_ks_n_grid = 0
  }
  if (!is.null(ks_covs)) {
    max_ks_n_grid = max(ks_n_grid)
    #! Grid of values	
    ks_grid = matrix(0, nrow = max_ks_n_grid, ncol = n_ks_covs)
    for (i in 1:n_ks_covs) {
      ks_covs_t_aux = ks_covs[, i][t_ind==1]
      ks_grid_aux = quantile(ks_covs_t_aux, probs = seq(1/ks_n_grid[i], 1, 1/ks_n_grid[i]))
      ks_grid_aux = c(ks_grid_aux, rep(0, max_ks_n_grid-ks_n_grid[i]))
      ks_grid[, i] = ks_grid_aux
    }	
    #! Auxiliary covariates
    ks_covs_aux = matrix(0, nrow = length(t_ind), ncol = max_ks_n_grid*n_ks_covs)
    for (i in 1:n_ks_covs) {
      k = (i-1)*max_ks_n_grid
      for (j in 1:max_ks_n_grid) {
        ks_covs_aux[, j+k][ks_covs[, i]<ks_grid[j, i]] = 1
      }
    }
  }
  
  #! Coeffs. of the obj. fun., cvec
  if (is.null(dist_mat)) {
    if (approximate == 1 | n_controls == 1) {
      cvec = -(1*rep(1, n_t*n_c))
    }
    else {
      cvec = c(-(1*rep(1, n_t*n_c)), rep(0, n_t)) 
    }
  }
  if (!is.null(dist_mat)) {
    if (approximate == 1 | n_controls == 1) {
      cvec = as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE))-(subset_weight*rep(1, n_t*n_c))
    }
    else {
      cvec = c(as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE))-(subset_weight*rep(1, n_t*n_c)), rep(0, n_t))
    }
    
  }
  
  #! Constraint matrix, Amat
  constraintmat_out = .constraintmatrix(t_ind, n_controls, total_groups,
                                        mom_covs, mom_tols, mom_targets,
                                        ks_covs, ks_covs_aux, ks_n_grid, ks_tols,
                                        exact_covs,
                                        near_exact_covs, near_exact_devs,
                                        fine_covs,
                                        near_fine_covs, near_fine_devs,
                                        near_covs, near_pairs, near_groups,
                                        far_covs, far_pairs, far_groups,
                                        use_controls,
                                        approximate)
  
  cnstrn_mat = constraintmat_out$cnstrn_mat
  bvec_7 = constraintmat_out$bvec_7
  bvec_8 = constraintmat_out$bvec_8
  rows_ind_far_pairs = constraintmat_out$rows_ind_far_pairs
  rows_ind_near_pairs = constraintmat_out$rows_ind_near_pairs
  
  # Constraint vector, bvec
  #! Parts 1 and 2
  if (approximate == 1 | n_controls == 1) {
    bvec = c(rep(n_controls, n_t), rep(1, n_c))
  }
  else {
    bvec = c(rep(0, n_t), rep(1, n_c))
  }
  
  #! Part 3: moments
  if (!is.null(mom_covs) & is.null(mom_targets)) {
    bvec = c(bvec, rep(0, 2*n_mom_covs))	
  }	
  
  #! Part 4: K-S
  if (!is.null(ks_covs)) {
    bvec = c(bvec, rep(0, 2*n_ks_covs*max_ks_n_grid))
  }
  
  #! Part 3b: Target
  if (!is.null(mom_covs) & !is.null(mom_targets)) {
    bvec = c(bvec, rep(0, 4*n_mom_covs))
  }
  
  #! Part 5: exact
  if (!is.null(exact_covs)) {
    bvec = c(bvec, rep(0, ncol(exact_covs))) 
  }	
  
  #! Part 6: near-exact
  if (!is.null(near_exact_covs)) {
    bvec = c(bvec, near_exact_devs) 
  }
  
  #! Part 7: fine
  if (!is.null(fine_covs)) {
    bvec = c(bvec, bvec_7) 
  }
  
  #! Part 8: near-fine
  if (!is.null(near_fine_covs)) {
    n_near_fine_covs = ncol(near_fine_covs)
    near_fine_devs_aux = NULL
    for (j in 1:n_near_fine_covs) {
      near_fine_cov = near_fine_covs[, j]
      near_fine_devs_aux = c(near_fine_devs_aux, rep(near_fine_devs[j], length(names(table(near_fine_cov))) ))
    }
    bvec_8_aux = rep(NA, length(bvec_8)*2)
    bvec_8_aux[1:length(bvec_8)] = -near_fine_devs_aux
    bvec_8_aux[(length(bvec_8)+1):(2*length(bvec_8))] = near_fine_devs_aux
    bvec = c(bvec, bvec_8_aux) 
  }
  
  #! Part 9: far
  if (!is.null(far_covs)) {
    n_far_covs = ncol(far_covs)
    for (j in 1:n_far_covs) {
      if (!is.null(far_groups)) {
        bvec = c(bvec, rep(0, 1))
      }
      if (!is.null(far_pairs) && rows_ind_far_pairs[[j]] != -1) {
        bvec = c(bvec, rep(0, length(table(rows_ind_far_pairs[[j]]))))
      }
    }
  }
  
  #! Part 10: near
  if (!is.null(near_covs)) {
    n_near_covs = ncol(near_covs)
    for (j in 1:n_near_covs) {
      if (!is.null(near_groups)) {
        bvec = c(bvec, rep(0, 1))
      }
      if (!is.null(near_pairs) && rows_ind_near_pairs[[j]] != -1) {
        bvec = c(bvec, rep(0, length(table(rows_ind_near_pairs[[j]]))))
      }
    }
  }
  
  #! Part 11: use controls
  if (!is.null(use_controls)) {
    bvec = c(bvec, sum(use_controls)) 
  }
  
  #! Part 12: total_groups
  if (!is.null(total_groups)) {
    if (!is.null(n_controls)) {
      bvec = c(bvec, total_groups*n_controls)
    }
    else {
      bvec = c(bvec, total_groups)
    }
    
  }
  
  # Upper bounds, ub
  #! Parts 1 and 2
  #! Part 3: moments
  #! Part 4: K-S
  if (approximate == 1 | n_controls == 1) {
    ub = rep(1, n_t*n_c)
  }	
  else {
    ub = c(rep(1, n_t*n_c), rep(1, n_t))
  }
  
  # Sense, sense
  #! Parts 1 and 2
  #! Part 3: moments
  #! Part 4: K-S
  if (approximate == 1 | n_controls == 1) {
    sense = c(rep("L", n_t), rep("L", n_c), rep("L", 2*n_mom_covs*(is.null(mom_targets))), rep("L", 2*n_ks_covs*max_ks_n_grid))
  }
  else {
    sense = c(rep("E", n_t), rep("L", n_c), rep("L", 2*n_mom_covs*(is.null(mom_targets))), rep("L", 2*n_ks_covs*max_ks_n_grid))
  }
  
  #! Part 3b: Target
  if (!is.null(mom_covs) & !is.null(mom_targets)) {
    sense = c(sense, rep("L", 4*n_mom_covs))
  }
  
  #! Part 5: exact
  if (!is.null(exact_covs)) {
    sense = c(sense, rep("E", ncol(exact_covs))) 
  }
  
  #! Part 6: near-exact
  if (!is.null(near_exact_covs)) {
    sense = c(sense, rep("L", ncol(near_exact_covs))) 
  }
  
  #! Part 7: fine
  if (!is.null(fine_covs)) {
    sense = c(sense, rep("E", length(bvec_7))) 
  }
  
  #! Part 8: near-fine
  if (!is.null(near_fine_covs)) {
    #sense = c(sense, rep(c("G", "L"), length(bvec_8))) 
    sense = c(sense, rep("G", length(bvec_8)), rep("L", length(bvec_8)))
  }
  
  #! Part 9: far
  if (!is.null(far_covs)) {
    n_far_covs = ncol(far_covs)
    for (j in 1:n_far_covs) {
      if (!is.null(far_groups)) { 
        sense = c(sense, rep("G", 1)) 
      }
      if (!is.null(far_pairs) && rows_ind_far_pairs[[j]] != -1) {
        sense = c(sense, rep("E", length(table(rows_ind_far_pairs[[j]]))))
      }
    }
  }
  
  #! Part 10: near
  if (!is.null(near_covs)) {
    n_near_covs = ncol(near_covs)
    for (j in 1:n_near_covs) {
      if (!is.null(near_groups)) { 
        sense = c(sense, rep("L", 1)) 
      }
      if (!is.null(near_pairs) && rows_ind_near_pairs[[j]] != -1) {
        sense = c(sense, rep("E", length(table(rows_ind_near_pairs[[j]]))))
      }
    }
  }
  
  #! Part 11: use controls
  if (!is.null(use_controls)) {
    sense = c(sense, "E") 
  }
  
  #! Part 12: total_groups
  if (!is.null(total_groups)) {
    sense = c(sense, "E")
  }
  
  # Variable types, vtype
  #! Parts 1 and 2
  #! Part 3: moments
  #! Part 4: K-S
  if (approximate == 1) {
    vtype = rep("C", n_t*n_c)
  }
  else if (n_controls == 1) {
    vtype = rep("B", n_t*n_c)
  }
  else {
    vtype = c(rep("B", n_t*n_c), rep("B", n_t))
  }
  
  #! Part 5: exact
  #! Part 6: near-exact
  #! Part 7: fine
  #! Part 8: near-fine
  #! Part 9: far
  #! Part 10: near
  #! Part 11: use controls
  #! Part 12: total_groups
  
  # c_index
  c_index = rep(1:n_c, n_t)	
  
  # Output
  return(list(n_t = n_t, n_c = n_c,  
              cvec = cvec, 
              Amat = cnstrn_mat, 
              bvec = bvec, 
              ub = ub, 
              sense = sense,
              vtype = vtype,
              c_index = c_index))
}

#! Generate the parameters for cardmatch
.problemparameters_cardmatch = function(t_ind, mom_covs, mom_tols, mom_targets, fine_covs) {
  
  #! Number of treated units and controls
  n_t = sum(t_ind)
  n_c = length(t_ind)-n_t
  
  #! Number of dec. vars.
  n_dec_vars = n_t+n_c
  
  #! Coeffs. of the obj. fun., cvec
  cvec = c(rep(1, n_t), rep(0, n_c))
  
  #! Constraint matrix, Amat
  row_ind_cur = 0
  #! Mom balance
  if (!is.null(mom_covs)) {
    rows_mom = NULL
    cols_mom = NULL
    vals_mom = NULL
    n_mom_covs = ncol(mom_covs)
    k = 1
    for (i in 1:n_mom_covs) {
      #! Treated
      rows_mom_plus = rep(row_ind_cur+k, n_t)
      rows_mom_minus = rep(row_ind_cur+k+1, n_t)  
      rows_mom = c(rows_mom, rows_mom_plus, rows_mom_minus)
      cols_mom = c(cols_mom, rep(1:n_t, 2))
      vals_plus = c(mom_covs[t_ind==1, i]-mom_targets[i]-mom_tols[i])
      vals_minus = c(mom_covs[t_ind==1, i]-mom_targets[i]+mom_tols[i])
      vals_mom = c(vals_mom, c(vals_plus, vals_minus)) 
      #! Controls
      rows_mom_plus = rep(row_ind_cur+k+2, n_c)
      rows_mom_minus = rep(row_ind_cur+k+3, n_c)  
      rows_mom = c(rows_mom, rows_mom_plus, rows_mom_minus)
      cols_mom = c(cols_mom, rep(n_t+(1:n_c), 2))
      vals_plus = c(mom_covs[t_ind==0, i]-mom_targets[i]-mom_tols[i])
      vals_minus = c(mom_covs[t_ind==0, i]-mom_targets[i]+mom_tols[i])
      vals_mom = c(vals_mom, c(vals_plus, vals_minus))       
      k = k+4
    }
    row_ind_cur = max(rows_mom)
  }
  #! Fine balance    
  if (!is.null(fine_covs)) {
    rows_fine = NULL
    cols_fine = NULL
    vals_fine = NULL
    n_fine_covs = ncol(fine_covs)
    k = 1
    for (i in 1:n_fine_covs) {
      fine_covs_cats = unique(fine_covs[, i])	
      for (j in fine_covs_cats) {
        cols_fine_aux = which(fine_covs[, i]==j)
        rows_fine = c(rows_fine, rep(row_ind_cur+k, length(cols_fine_aux)))
        cols_fine = c(cols_fine, cols_fine_aux)
        vals_fine = c(vals_fine, c(rep(1, sum(cols_fine_aux<=n_t)), rep(-1, sum(cols_fine_aux>n_t))))
        k = k+1
      }
    }
    row_ind_cur = max(rows_fine)
  }
  #! Equal number of matched controls and matched treated units
  if (!is.null(mom_covs) & is.null(fine_covs)) {
    rows_equal = rep(row_ind_cur+1, n_dec_vars)
    cols_equal = 1:n_dec_vars
    vals_equal = c(rep(1, n_t), rep(-1, n_c))
  }
  #! Put together
  rows_ind = NULL
  cols_ind = NULL
  vals = NULL
  #! Mom balance      
  if (!is.null(mom_covs)) {
    rows_ind = c(rows_ind, rows_mom)
    cols_ind = c(cols_ind, cols_mom)
    vals = c(vals, vals_mom)
  }
  #! Fine balance      
  if (!is.null(fine_covs)) {
    rows_ind = c(rows_ind, rows_fine)
    cols_ind = c(cols_ind, cols_fine)
    vals = c(vals, vals_fine)
  }
  #! Equal number of matched controls and matched treated units
  if (!is.null(mom_covs) & is.null(fine_covs)) {
    rows_ind = c(rows_ind, rows_equal)
    cols_ind = c(cols_ind, cols_equal)
    vals = c(vals, vals_equal)
  }    
  #! 
  aux = cbind(rows_ind, cols_ind, vals)[order(cols_ind), ]
  Amat = simple_triplet_matrix(i = aux[, 1], j = aux[, 2], v = aux[, 3])
  
  #! Constraint vector, bvec
  bvec = NULL
  #! Mom balance
  if (!is.null(mom_covs)) {
    bvec_mom = rep(0, length(unique(rows_mom)))
    bvec = c(bvec, bvec_mom)
  }
  #! Fine balance
  if (!is.null(fine_covs)) {
    bvec_fine = rep(0, length(unique(rows_fine)))
    bvec = c(bvec, bvec_fine)
  }
  #! Equal number of matched controls and matched treated units
  if (!is.null(mom_covs) & is.null(fine_covs)) {
    bvec = c(bvec, 0)
  }
  
  #! Sense, sense
  sense = NULL
  #! Mom balance
  if (!is.null(mom_covs)) {
    sense_covs = rep(c("L", "G", "L", "G"), length(unique(rows_mom))/4)
    sense = c(sense, sense_covs) 
  }
  #! Fine balance
  if (!is.null(fine_covs)) {
    sense_fine = rep("E", length(unique(rows_fine)))
    sense = c(sense, sense_fine) 
  }
  if (!is.null(mom_covs) & is.null(fine_covs)) {
    sense_equal = c("E")
    sense = c(sense, sense_equal) 
  }
  
  #! Variable types, vtype
  vtype = rep("B", n_dec_vars)
  
  # Output
  return(list(n_t = n_t, 
              n_c = n_c,  
              n_dec_vars = n_dec_vars,
              cvec = cvec, 
              Amat = Amat, 
              bvec = bvec, 
              sense = sense,
              vtype = vtype))
  
}

# Solves relaxation problem
.relaxation_b = function(n_t, n_c, coef, dist_mat, subset_weight, solver, round_cplex, trace) {
  
  n_tot = n_t*n_c
  #! Part 1
  row_ind_1 = sort(rep(1:n_t, n_c))
  col_ind_1 = 1:n_tot
  val_1 = rep(1, n_tot)
  sense_1 = rep("L", n_t)
  bvec_1 = rep(1, n_t)
  
  row_ind_cur  = max(row_ind_1)
  
  #! Part 2
  row_ind_2 = rep(1:n_c, n_t)+row_ind_cur
  col_ind_2 = 1:n_tot
  val_2 = rep(1, n_tot)
  sense_2 = rep("L", n_c)
  bvec_2 = rep(1, n_c)
  
  row_ind = c(row_ind_1, row_ind_2)
  col_ind = c(col_ind_1, col_ind_2)
  vals = c(val_1, val_2)
  sense = c(sense_1, sense_2)
  bvec = c(bvec_1, bvec_2)
  
  aux = cbind(row_ind, col_ind, vals)[order(col_ind), ]
  Amat = simple_triplet_matrix(i = aux[, 1], j = aux[, 2], v = aux[, 3])
  
  ub = Inf
  vtype = rep("B", n_tot)
  cvec = coef
  
  
  #### SOLVER PART #######
  
  if (solver == "cplex"){
    #library(Rcplex)
    if (requireNamespace('Rcplex', quietly = TRUE)) {
      ptm = proc.time()
      out = Rcplex::Rcplex(cvec, Amat, bvec, ub = ub, sense = sense, vtype = vtype, n = 1,
                           control = list(trace = trace, round = round_cplex), objsense = "max")
      time = (proc.time()-ptm)[3]
      
      if (out$status==108) {
        cat(format("  Error: time limit exceeded, no integer solution!"), "\n")
        obj = 0
        sol = NULL
      } else if (is.na(out$obj)) {
        cat(format("  Error: problem infeasible!"), "\n")
        obj = 0
        sol = NULL
      } else {
        sol = out$xopt
        if(is.null(dist_mat)) {
          obj = sum(-1*(rep(1, n_tot)) * out$xopt)
        } else {
          obj = sum((as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE))-(subset_weight*rep(1, n_t*n_c)))*out$xopt)
        }
      }
    } else {
      stop('suggested package not installed')
    }
  }
  
  if (solver == "gurobi") {
    #library(gurobi)
    if (requireNamespace('gurobi', quietly = TRUE)) {
      model = list()
      model$obj = cvec
      model$A = Amat
      model$sense = rep(NA, length(sense))
      model$sense[sense=="E"] = '='
      model$sense[sense=="L"] = '<='
      model$sense[sense=="G"] = '>='
      model$rhs = bvec
      model$vtypes = vtype
      model$ub = ub
      model$modelsense = "max"
      
      params = list(OutputFlag = trace)
      ptm = proc.time()
      out = gurobi::gurobi(model, params)
      time = (proc.time()-ptm)[3]
      
      if (out$status == "INFEASIBLE") {
        cat(format("  Error: problem infeasible!"), "\n")
        obj = 0
        sol = NULL
      }
      
      if (out$status == "OPTIMAL") {
        sol = out$x
        if(is.null(dist_mat)) {
          obj = sum(-1*(rep(1, n_tot)) * out$x)
        } else {
          obj = sum((as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE))-(subset_weight*rep(1, n_t*n_c)))*out$x)
        }
      }
    } else {
      stop('suggested package not installed')
    }
  }
  
  if (solver == "symphony") {
    #library(Rsymphony)
    if (requireNamespace('Rsymphony', quietly = TRUE)) {
      dir = rep(NA, length(sense))
      dir[sense=="E"] = '=='
      dir[sense=="L"] = '<='
      dir[sense=="G"] = '>='
      
      bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
                   upper = list(ind=c(1:length(ub)), val=ub))
      
      ptm = proc.time()
      out= Rsymphony::Rsymphony_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = vtype, max = TRUE)
      time = (proc.time()-ptm)[3]
      
      if (out$status!=0) {
        cat(format("  Error: problem infeasible!"), "\n")
        obj = 0
        sol = NULL
      } else {
        sol = out$solution
        if(is.null(dist_mat)) {
          obj = sum(-1*(rep(1, n_tot)) * out$solution)
        } else {
          obj = sum((as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE))-(subset_weight*rep(1, n_t*n_c)))*out$solution)
        }
      }
    } else {
      stop('suggested package not installed')
    }
  }
  
  # GLPK
  if (solver == "glpk") {
    #library(Rglpk)
    dir = rep(NA, length(sense))
    dir[sense=="E"] = '=='
    dir[sense=="L"] = '<='
    dir[sense=="G"] = '>='
    
    bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
                 upper = list(ind=c(1:length(ub)), val=ub))
    
    ptm = proc.time()
    out= Rglpk_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = vtype, max = TRUE)
    time = (proc.time()-ptm)[3]
    
    if (out$status!=0) {
      cat(format("  Error: problem infeasible!"), "\n")
      obj = 0
      sol = NULL
    } else {
      sol = out$solution
      if(is.null(dist_mat)) {
        obj = sum(-1*(rep(1, n_tot)) * out$solution)
      } else {
        obj = sum((as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE))-(subset_weight*rep(1, n_t*n_c)))*out$solution)
      }
    }
    
  }
  
  return(list(sol = sol, obj = obj, time = time))
  
}

# Solves relaxation problem
.relaxation_n = function(n_tot, coef, dist_mat, subset_weight, solver, round_cplex, trace) {
  
  n_dec = (n_tot*(n_tot-1))-sum(1:(n_tot-1))
  #! Nonbipartite matching constraints
  rows_nbm = sort(rep(1:n_tot, n_tot-1))
  temp = matrix(0, nrow = n_tot, ncol = n_tot)
  temp[lower.tri(temp)] = 1:n_dec
  temp = temp+t(temp)
  diag(temp) = NA
  cols_nbm = as.vector(t(temp))
  cols_nbm = cols_nbm[!is.na(cols_nbm)]
  vals_nbm = rep(1, (n_tot-1)*n_tot)
  
  bvec = rep(1,  length(table(rows_nbm)))
  sense = rep("L", length(table(rows_nbm)))
  
  aux = cbind(rows_nbm, cols_nbm, vals_nbm)[order(cols_nbm), ]
  Amat = simple_triplet_matrix(i = aux[, 1], j = aux[, 2], v = aux[, 3])
  
  ub = Inf
  vtype = rep("B", n_dec)
  cvec = coef
  
  #### SOLVER PART #######
  
  if (solver == "cplex"){
    #library(Rcplex)
    if (requireNamespace('Rcplex', quietly = TRUE)) {
      ptm = proc.time()
      out = Rcplex::Rcplex(cvec, Amat, bvec, ub = ub, sense = sense, vtype = vtype, n = 1,
                           control = list(trace = trace, round = round_cplex), objsense = "max")
      time = (proc.time()-ptm)[3]
      
      if (out$status==108) {
        cat(format("  Error: time limit exceeded, no integer solution!"), "\n")
        obj = 0
        sol = NULL
      } else if (is.na(out$obj)) {
        cat(format("  Error: problem infeasible!"), "\n")
        obj = 0
        sol = NULL
      } else {
        sol = out$xopt
        obj = sum((t(dist_mat)[lower.tri(dist_mat)]-(subset_weight*rep(1, n_dec))) * out$xopt)
      } 
    } else {
      stop('suggested package not installed')
    }
  }
  
  if (solver == "gurobi") {
    #library(gurobi)
    if (requireNamespace('gurobi', quietly=TRUE)) {
      model = list()
      model$obj = cvec
      model$A = Amat
      model$sense = rep(NA, length(sense))
      model$sense[sense=="E"] = '='
      model$sense[sense=="L"] = '<='
      model$sense[sense=="G"] = '>='
      model$rhs = bvec
      model$vtypes = vtype
      model$ub = ub
      model$modelsense = "max"
      
      params = list(OutputFlag = trace)
      ptm = proc.time()
      out = gurobi::gurobi(model, params)
      time = (proc.time()-ptm)[3]
      
      if (out$status == "INFEASIBLE") {
        cat(format("  Error: problem infeasible!"), "\n")
        obj = 0
        sol = NULL
      }
      
      if (out$status == "OPTIMAL") {
        sol = out$x
        obj = sum((t(dist_mat)[lower.tri(dist_mat)]-(subset_weight*rep(1, n_dec))) * out$x)
      }
    } else {
      stop('suggested package not installed')
    }
  }
  
  if (solver == "symphony") {
    #library(Rsymphony)
    if (requireNamespace('Rsymphony', quietly = TRUE)) {
      dir = rep(NA, length(sense))
      dir[sense=="E"] = '=='
      dir[sense=="L"] = '<='
      dir[sense=="G"] = '>='
      
      bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
                   upper = list(ind=c(1:length(ub)), val=ub))
      
      ptm = proc.time()
      out= Rsymphony::Rsymphony_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = vtype, max = TRUE)
      time = (proc.time()-ptm)[3]
      
      if (out$status!=0) {
        cat(format("  Error: problem infeasible!"), "\n")
        obj = 0
        sol = NULL
      } else {
        sol = out$solution
        obj = sum((t(dist_mat)[lower.tri(dist_mat)]-(subset_weight*rep(1, n_dec))) * out$solution)
      }
    } else {
      stop('suggested package not installed')
    }
  }
  
  # GLPK
  if (solver == "glpk") {
    #library(Rglpk)
    #dir = rep(NA, length(sense))
    #dir[sense=="E"] = '=='
    #dir[sense=="L"] = '<='
    #dir[sense=="G"] = '>='
    #
    #bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
    #             upper = list(ind=c(1:length(ub)), val=ub))
    #
    #ptm = proc.time()
    #out= Rglpk_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = vtype, max = TRUE)
    #time = (proc.time()-ptm)[3]
    #
    #if (out$status!=0) {
    #  cat(format("  Error: problem infeasible!"), "\n")
    #  obj = 0
    #  sol = NULL
    #} else {
    #  sol = out$solution
    #  obj = sum((t(dist_mat)[lower.tri(dist_mat)]-(subset_weight*rep(1, n_dec))) * out$solution)
    #}
    ptm = proc.time()
    Amat = matrix(0, nrow = n_tot, ncol = n_tot)
    res = matrix(0, nrow = n_tot, ncol = n_tot)
    
    Amat[upper.tri(Amat)] = coef
    
    max_edge = max(Amat)
    while (max_edge > 0) {
      row = which(Amat == max_edge, arr.ind=TRUE)[1,1]
      col = which(Amat == max_edge, arr.ind=TRUE)[1,2]
      res[row,col] = 1
      Amat[row,] = 0
      Amat[,row] = 0
      Amat[col,] = 0
      Amat[,col] = 0
      max_edge = max(Amat)
    }
    sol = res[upper.tri(res)]
    obj = sum((t(dist_mat)[lower.tri(dist_mat)]-(subset_weight*rep(1, n_dec))) * sol)
    time = (proc.time()-ptm)[3]
  }
  
  
  return(list(sol = sol, obj = obj, time = time))
  
}

#! Rank based Mahalanobis distance, from Paul Rosenbaum's Design of Observational Studies, p. 251
.smahal = function(z, X) { 
  X = as.matrix(X) 
  n = dim(X)[1] 
  rownames(X) = 1:n 
  k = dim(X)[2] 
  m = sum(z) 
  for (j in 1:k) X[,j] = rank(X[,j]) 
  cv = cov(X) 
  vuntied = var(1:n)
  #! ***PENDING: correct this	
  diag(cv)[diag(cv) == 0] = .01
  rat = sqrt(vuntied/diag(cv)) 
  cv = diag(rat)%*%cv%*%diag(rat) 
  out = matrix(NA,m,n-m) 
  Xc = X[z == 0, ] 
  Xt = X[z == 1, ] 
  rownames(out) = rownames(X)[z==1] 
  colnames(out) = rownames(X)[z==0] 
  #library(MASS) 
  icov = ginv(cv) 
  for (i in 1:m) {
    out[i,] = mahalanobis(Xc,Xt[i,],icov,inverted=TRUE)
  } 
  out
}
