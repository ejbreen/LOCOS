

for (i in 1:3){
  print(paste('Generating dataset', i))
  
  C_pop <- genData(60000, C_pop_def)
  T_pop <- genData(2000, T_pop_def)
  
  C_pop <- subset.data.frame(C_pop, select = -c(B, POLY, HEAD, APP))
  T_pop <- subset.data.frame(T_pop, select = -c(B, POLY, HEAD, APP))
  
  
  write.csv(C_pop[,2:ncol(C_pop)], file = paste(
    "Data/C_pop",i,".csv"))
  write.csv(T_pop[,2:ncol(T_pop)], file = paste(
    "Data/T_pop",i,".csv"))
}

