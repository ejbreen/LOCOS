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

#################################
# Example 3: optimal subset matching
#################################
# Optimal subset matching pursues two competing goals at
# the same time: to minimize the total sum of covariate distances
# while matching as many observations as possible. The trade-off
# between these two goals is regulated by the parameter subset_weight
# (see Rosenbaum 2012 and Zubizarreta et al. 2013 for a discussion).
# Here the balance requirements are mean balance, near-fine balance
# and near-exact matching for different covariates.
# Again, the solver used is glpk with the approximate option.
# Treatment indicator
t_ind = treatment
# Matrix of covariates
X_mat = cbind(age, education, black, hispanic, married, nodegree, re74, re75)
# Distance matrix
dist_mat = distmat(t_ind, X_mat)
# Subset matching weight
subset_weight = median(dist_mat)
# Moment balance: constrain differences in means to be at most .05 standard deviations apart
mom_covs = cbind(age, education, black, hispanic, married, nodegree, re74, re75)
mom_tols = round(absstddif(mom_covs, t_ind, .05), 2)
mom = list(covs = mom_covs, tols = mom_tols)
# Near-fine balance
near_fine_covs = cbind(married, nodegree)
near_fine_devs = rep(5, 2)
near_fine = list(covs = near_fine_covs, devs = near_fine_devs)
# Near-exact matching
near_exact_covs = cbind(black, hispanic)
near_exact_devs = rep(5, 2)
near_exact = list(covs = near_exact_covs, devs = near_exact_devs)
# Solver options
t_max = 60*5
solver = "gurobi" #switched to gurobi
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate,
round_cplex = 0, trace_cplex = 0)
# Match
out = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight,
mom = mom, near_fine = near_fine, near_exact = near_exact, solver = solver)
# Indices of the treated units and matched controls
t_id = out$t_id
c_id = out$c_id
# Time
out$time/60
# Matched group identifier (who is matched to whom)
out$group_id
# Assess mean balance (note here we are getting an approximate solution)
meantab(X_mat, t_ind, t_id, c_id)
# Assess fine balance
for (i in 1:ncol(near_fine_covs)) {
print(finetab(near_fine_covs[, i], t_id, c_id))
}
# Assess exact matching balance
for (i in 1:ncol(near_exact_covs)) {
print(table(near_exact_covs[t_id, i]==near_exact_covs[c_id, i]))
}