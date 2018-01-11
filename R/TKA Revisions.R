# Defining new data set definitions with revision risk built in
# The 
T_rev_random_def = T_pop_def 
T_rev_causal_def = T_pop_def
factors = T_pop_def$varname

Causal_formula = 'round(1/(1+exp(-(5.3*B_COP - 5.295))))'

T_rev_random_def <- defData(T_rev_random_def,
                            varname = '5yr_Rev_Rate',
                            dist = 'binary',
                            formula = .005)

T_rev_causal_def <- defData(T_rev_causal_def,
                            varname = '5yr_Rev_Rate',
                            dist = 'nonrandom',
                            formula = Causal_formula)

for (i in 1:1){
  print(paste('Generating dataset', i))
  
  T_rev_random <- genData(2000, T_rev_random_def)
  T_rev_causal <- genData(2000, T_rev_causal_def)
  
  T_rev_random <- subset.data.frame(T_rev_random, 
                                    select = -c(B, POLY, HEAD, APP))
  T_rev_causal <- subset.data.frame(T_rev_causal, 
                                    select = -c(B, POLY, HEAD, APP))
  
  
  write.csv(T_rev_random[,2:ncol(C_pop)], 
            file = paste("Data/C_pop",i,".csv"))
  write.csv(T_rev_causal[,2:ncol(T_pop)], 
            file = paste("Data/T_pop",i,".csv"))
}


