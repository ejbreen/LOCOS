ggplot(data = C_pop, aes(x=AGE))+geom_histogram(binwidth = 1)


head(T_pop[["AGE"]], 10)

summary(T_pop)
