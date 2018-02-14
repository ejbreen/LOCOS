ggplot(data = head(T_pop, 100), aes(x=))+geom_histogram(binwidth = 1)

plot(x=seq(1, 1000), y=out$group_id)

head(T_pop[["AGE"]], 10)

summary(T_pop)

colSums(T_pop)

out <- bmatch_fnc(Tot_pop)

Tot_pop[c_id]

Tot_pop_2 <- Tot_pop[-c(c_id),]

out_2 <- bmatch_fnc(Tot_pop_2)



matched <- subset(Tot_pop, subset = IMPLANT==1)
matched <- rbind(matched, Tot_pop[out$c_id])
