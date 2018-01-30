# Required packages:    Rglpk
#                       designmatch 
#                           install.packages("designmatch")
#                       gurobi      

# Windows
# install.packages('C:/gurobi752/win64/R/gurobi_7.5-2.zip', repos = NULL)
# Chromebook
install.packages('/opt/gurobi752/linux64/R/gurobi_7.5-2_R_x86_64-pc-linux-gnu.tar.gz', repos = NULL)

library(gurobi)
library(Rglpk)
library(designmatch)

data(lalonde)
attach(lalonde)
quantiles = function(covar, n_q) {
    p_q = seq(0, 1, 1 / n_q)
	val_q = quantile(covar, probs = p_q, na.rm = TRUE)
	covar_out = rep(NA, length(covar))
    for (i in 1:n_q) {
        if (i == 1) { covar_out[covar < val_q[i + 1]] = i }
        if (i > 1 & i < n_q) { covar_out[covar >= val_q[i] & covar < val_q[i + 1]] = i }
        if (i == n_q) { covar_out[covar >= val_q[i] & covar < val_q[i + 1]] = i}
	}
	covar_out
}

age_5 = quantiles(age, 5)
education_5 = quantiles(education, 5)
re74_5 = quantiles(re74, 5)
re75_5 = quantiles(re75, 5)

# Treatment indicator
t_ind = treatment

# Fine Balance 
fine_covs = cbind(black, hispanic, married, nodegree, age_5, education_5, re74_5, re75_5)
fine = list(covs = fine_covs)

# Solver options
t_max = 60 * 5
solver = 'gurobi'
solver = list(name = solver, t_max = t_max, approximate = 0, round_cplex = 0, trace = 0)

# Match it
out_1 = cardmatch(t_ind, fine = fine, solver = solver)



