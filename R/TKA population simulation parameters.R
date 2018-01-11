
print(paste('Defining population simulation paramaters'))

Bearing_levels <- c('Metal-on-plastic','Ceramic-on-plastic','Ceramic-on-ceramic','Dual_mobility')
Poly_Type_levels <- c('UHWMPE','XLPE','Antioxidant_XLPE')
Head_Size_levels <- c('22mm','28mm','32mm','36mm','40mm','44mm')
Approach_levels <- c('anterior','anterolateral','posterior','transtrochanteric')



C_pop_def <- defData(varname = 'IMPLANT', 
                     dist = 'nonrandom', formula = 0)
C_pop_def <- defData(C_pop_def, varname = 'B', 
                     dist = 'categorical', formula = '.449;.498;.007;.046')
C_pop_def <- defData(C_pop_def, varname = 'B_MOP', 
                     dist = 'nonrandom', formula = 'B==1')
C_pop_def <- defData(C_pop_def, varname = 'B_COP', 
                     dist = 'nonrandom', formula = 'B==2')
C_pop_def <- defData(C_pop_def, varname = 'B_COC', 
                     dist = 'nonrandom', formula = 'B==3')
C_pop_def <- defData(C_pop_def, varname = 'B_DM', 
                     dist = 'nonrandom', formula = 'B==4')
C_pop_def <- defData(C_pop_def, varname = 'POLY', 
                     dist = 'categorical', formula = '.003;.814;.183')
C_pop_def <- defData(C_pop_def, varname = 'POLY_UHWMPE', 
                     dist = 'nonrandom', formula = '(POLY==1) & (B!=3)')
C_pop_def <- defData(C_pop_def, varname = 'POLY_XPLE', 
                     dist = 'nonrandom', formula = '(POLY==2) & (B!=3)')
C_pop_def <- defData(C_pop_def, varname = 'POLY_A_XPLE', 
                     dist = 'nonrandom', formula = '(POLY==3) & (B!=3)')
C_pop_def <- defData(C_pop_def, varname = 'HEAD', 
                     dist = 'categorical', formula = '.000;.007;.245;.651;.090;.007')
C_pop_def <- defData(C_pop_def, varname = 'HEAD_22mm', 
                     dist = 'nonrandom', formula = 'HEAD==1')
C_pop_def <- defData(C_pop_def, varname = 'HEAD_28mm', 
                     dist = 'nonrandom', formula = 'HEAD==2')
C_pop_def <- defData(C_pop_def, varname = 'HEAD_32mm', 
                     dist = 'nonrandom', formula = 'HEAD==3')
C_pop_def <- defData(C_pop_def, varname = 'HEAD_36mm', 
                     dist = 'nonrandom', formula = 'HEAD==4')
C_pop_def <- defData(C_pop_def, varname = 'HEAD_40mm', 
                     dist = 'nonrandom', formula = 'HEAD==5')
C_pop_def <- defData(C_pop_def, varname = 'HEAD_44mm', 
                     dist = 'nonrandom', formula = 'HEAD==6')
C_pop_def <- defData(C_pop_def, varname = 'APP', 
                     dist = 'categorical', formula = '.437;.223;.332;.008')
C_pop_def <- defData(C_pop_def, varname = 'APP_anterior',
                     dist = 'nonrandom', formula = 'APP==1')
C_pop_def <- defData(C_pop_def, varname = 'APP_anterolateral',
                     dist = 'nonrandom', formula = 'APP==2')
C_pop_def <- defData(C_pop_def, varname = 'APP_posterior',
                     dist = 'nonrandom', formula = 'APP==3')
C_pop_def <- defData(C_pop_def, varname = 'APP_transtrochanteric',
                     dist = 'nonrandom', formula = 'APP==4')
C_pop_def <- defData(C_pop_def, varname = 'S_VOLLUME', 
                     dist = 'poisson', formula = 110.5)
C_pop_def <- defData(C_pop_def, varname = 'AGE', 
                     dist = 'normal', formula = 63.99, variance = 11.17*11.17)
C_pop_def <- defData(C_pop_def, varname = 'FEMALE', 
                     dist = 'binary', formula = .546)
C_pop_def <- defData(C_pop_def, varname = 'BMI', 
                     dist = 'normal', formula = 30.25, variance = 6.28*6.28)


T_pop_def <- defData(varname = 'IMPLANT',
                     dist = 'nonrandom', formula = 1)
T_pop_def <- defData(T_pop_def, varname = 'B',
                     dist = 'categorical', formula = '.264;.736;.000;.000')
T_pop_def <- defData(T_pop_def, varname = 'B_MOP', 
                     dist = 'nonrandom', formula = 'B==1')
T_pop_def <- defData(T_pop_def, varname = 'B_COP', 
                     dist = 'nonrandom', formula = 'B==2')
T_pop_def <- defData(T_pop_def, varname = 'B_COC', 
                     dist = 'nonrandom', formula = 'B==3')
T_pop_def <- defData(T_pop_def, varname = 'B_DM', 
                     dist = 'nonrandom', formula = 'B==4')
T_pop_def <- defData(T_pop_def, varname = 'POLY',
                     dist = 'categorical', formula = '.0007;.9993;.0000')
T_pop_def <- defData(T_pop_def, varname = 'POLY_UHWMPE', 
                     dist = 'nonrandom', formula = '(POLY==1) & (B!=3)')
T_pop_def <- defData(T_pop_def, varname = 'POLY_XPLE', 
                     dist = 'nonrandom', formula = '(POLY==2) & (B!=3)')
T_pop_def <- defData(T_pop_def, varname = 'POLY_A_XPLE', 
                     dist = 'nonrandom', formula = '(POLY==3) & (B!=3)')
T_pop_def <- defData(T_pop_def, varname = 'HEAD',
                     dist = 'categorical', formula = '.000;.007;.259;.600;.123;.011')
T_pop_def <- defData(T_pop_def, varname = 'HEAD_22mm', 
                     dist = 'nonrandom', formula = 'HEAD==1')
T_pop_def <- defData(T_pop_def, varname = 'HEAD_28mm', 
                     dist = 'nonrandom', formula = 'HEAD==2')
T_pop_def <- defData(T_pop_def, varname = 'HEAD_32mm', 
                     dist = 'nonrandom', formula = 'HEAD==3')
T_pop_def <- defData(T_pop_def, varname = 'HEAD_36mm', 
                     dist = 'nonrandom', formula = 'HEAD==4')
T_pop_def <- defData(T_pop_def, varname = 'HEAD_40mm', 
                     dist = 'nonrandom', formula = 'HEAD==5')
T_pop_def <- defData(T_pop_def, varname = 'HEAD_44mm', 
                     dist = 'nonrandom', formula = 'HEAD==6')
T_pop_def <- defData(T_pop_def, varname = 'APP',
                     dist = 'categorical', formula = '.227;.256;.505;.012')
T_pop_def <- defData(T_pop_def, varname = 'APP_anterior',
                     dist = 'nonrandom', formula = 'APP==1')
T_pop_def <- defData(T_pop_def, varname = 'APP_anterolateral',
                     dist = 'nonrandom', formula = 'APP==2')
T_pop_def <- defData(T_pop_def, varname = 'APP_posterior',
                     dist = 'nonrandom', formula = 'APP==3')
T_pop_def <- defData(T_pop_def, varname = 'APP_transtrochanteric',
                     dist = 'nonrandom', formula = 'APP==4')
T_pop_def <- defData(T_pop_def, varname = 'S_VOLLUME',
                     dist = 'poisson', formula = 38.4)
T_pop_def <- defData(T_pop_def,  varname = 'AGE',
                     dist = 'normal', formula = 64.02, variance = 11.44*11.44)
T_pop_def <- defData(T_pop_def, varname = 'FEMALE',
                     dist = 'binary', formula = .537)
T_pop_def <- defData(T_pop_def, varname = 'BMI',
                     dist = 'normal', formula = 31.40, variance = 7.05*7.05)


write.csv(C_pop_def, file = "Data/C_pop_def.csv")
write.csv(T_pop_def, file = "Data/T_pop_def.csv")



