# required for gurobi
install.packages('slam', repos = 'https://cloud.r-project.org/') # required for gurobi
library(slam)
install.packages('/sw/arc/centos7/gurobi/gurobi751/linux64/R/gurobi_7.5-1_R_x86_64-pc-linux-gnu.tar.gz', repos = NULL)
library(gurobi)

# still having trouble with this one.  It requires Rglpk so I need glpk installed but it isn't available on flux
install.packages('designmatch', repos = 'https://cloud.r-project.org/')
