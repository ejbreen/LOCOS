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
# Example 2: minimum distance matching
#################################
# The goal here is to minimize the total of distances between matched pairs. In
# this example there are no covariate balance requirements. Again, the solver
# used is glpk with the approximate option
# Treatment indicator
t_ind = treatment
# Matrix of covariates
X_mat = cbind(age, education, black, hispanic, married, nodegree, re74, re75)
# Distance matrix
dist_mat = distmat(t_ind, X_mat)
# Subset matching weight
subset_weight = NULL
# Total pairs to be matched
total_groups = sum(t_ind)
# Solver options
t_max = 60*5
solver = "gurobi"
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate,
round_cplex = 0, trace_cplex = 0)
# Match
out = bmatch(t_ind = t_ind, dist_mat = dist_mat, total_groups = total_groups,
solver = solver)
# Indices of the treated units and matched controls
t_id = out$t_id
c_id = out$c_id
# Total of distances between matched pairs
out$obj_total
# Assess mean balance
meantab(X_mat, t_ind, t_id, c_id)