# Required packages:    Rglpk
#                       designmatch
#                       gurobi      'C:/gurobi752/win64/R/gurobi_7.5-2.zip'

install.packages('C:/gurobi752/win64/R/gurobi_7.5-2.zip', repos = NULL)

library(gurobi)
library(Rglpk)
library(designmatch)


# Load and attach data
data(lalonde)
attach(lalonde)

# Example 1: cardinality matching
#################################
# Cardinality matching finds the largest matched sample of pairs that meets balance
# requirements. Here the balance requirements are mean balance, fine balance and
# exact matching for different covariates. The solver used is glpk with the
# approximate option.
# Treatment indicator
t_ind = treatment
# Distance matrix
dist_mat = NULL
# Subset matching weight
subset_weight = 1
# Moment balance: constrain differences in means to be at most .05 standard deviations apart
mom_covs = cbind(age, education, black, hispanic, married, nodegree, re74, re75)
mom_tols = round(absstddif(mom_covs, t_ind, .05), 2)
mom = list(covs = mom_covs, tols = mom_tols)
# Fine balance
fine_covs = cbind(black, hispanic, married, nodegree)
fine = list(covs = fine_covs)
# Exact matching
exact_covs = cbind(black)
exact = list(covs = exact_covs)
# Solver options
t_max = 60*5
solver = "gurobi"
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate,
              round_cplex = 0, trace = 0)
# Match
out = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight,
             mom = mom, fine = fine, exact = exact, solver = solver)
# Indices of the treated units and matched controls
t_id = out$t_id
c_id = out$c_id
# Time
out$time/60
# Matched group identifier (who is matched to whom)
out$group_id
# Assess mean balance
meantab(mom_covs, t_ind, t_id, c_id)
# Assess fine balance (note here we are getting an approximate solution)
for (i in 1:ncol(fine_covs)) {
  print(finetab(fine_covs[, i], t_id, c_id))
}
# Assess exact matching balance
table(exact_covs[t_id]==exact_covs[c_id])
