ggplot(data = head(T_pop, 100), aes(x=))+geom_histogram(binwidth = 1)

plot(x=seq(1, 1000), y=out$group_id)

head(T_pop[["AGE"]], 10)

summary(T_pop)

colSums(T_pop)

out <- bmatch_fnc(Tot_pop)

Tot_pop[c_id]

out_x <- bmatch_fnc(Tot_pop_x)



###########################
source('R/Bmatch_fnc.R')
Tot_pop <- combine_pops(T_pop, C_pop,.1)

matched <- subset(Tot_pop, subset = IMPLANT==1)
out <- bmatch_fnc(Tot_pop)
matched <- rbind(matched, Tot_pop[out$c_id])
Tot_pop_x <- Tot_pop[-c(out$c_id),]
for (i in 2:5){
  out_x <- bmatch_fnc(Tot_pop_x)
  matched <- rbind(matched, Tot_pop_x[out_x$c_id])
  Tot_pop_x <- Tot_pop_x[-c(out_x$c_id)]
  print(i)
}


t_ind <- Tot_pop$IMPLANT
attach(Tot_pop)
x_mat <- cbind(B_MOP, B_COP, B_COC, B_DM, POLY_UHWMPE, POLY_XPLE, POLY_A_XPLE, HEAD_22mm, HEAD_28mm, HEAD_32mm, HEAD_36mm, HEAD_40mm, HEAD_44mm,APP_anterior, APP_anterolateral, APP_posterior, APP_transtrochanteric,S_VOLLUME, AGE, FEMALE, BMI)
c_all <- subset(Tot_pop, subset = IMPLANT == 0)$id + sum(t_ind)
t_id <- subset(matched, subset = IMPLANT == 1)$id
c_id <- subset(matched, subset = IMPLANT == 0)$id + sum(t_ind) #adjusting index for combined set
meantab(x_mat, t_ind, t_id, c_id)
meantab(x_mat, t_ind, t_id, c_all)


