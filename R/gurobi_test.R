install.packages('C:/gurobi752/win64/R/gurobi_7.5-2.zip', repos = NULL)
library(gurobi)

model <- list()

model$A <- matrix(c(1, 2, 3, 1, 1, 0), nrow = 2, ncol = 3, byrow = T)
model$obj <- c(1, 1, 2)
model$modelsense <- 'max'
model$rhs <- c(4, 1)
model$sense <- c('<=', '>=')
model$vtype <- 'B'

params <- list(OutputFlag = 0)

result <- gurobi(model, params)

result$x

