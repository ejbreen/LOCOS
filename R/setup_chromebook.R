install.packages('simstudy', repos = 'https://cloud.r-project.org/')
library(simstudy)

install.packages('slam', repos = 'https://cloud.r-project.org/') # required for gurobi
library(slam)
install.packages('MASS',  repos = 'https://cloud.r-project.org/')
library(MASS)
install.packages('lattice',  repos = 'https://cloud.r-project.org/')
library(lattice)

# the location of gurobi on the pixelbook
install.packages('/opt/gurobi752/linux64/R/gurobi_7.5-2_R_x86_64-pc-linux-gnu.tar.gz', repos = NULL)
library(gurobi)
