source('R/setup_home.R')

library(survival)

Tot_pop <- read.csv('R/sample_set_1000.csv')

Tot_pop$followup_date <- 15
Tot_pop$followup_rev <- Tot_pop$ttr < Tot_pop$followup_date

cox1 <- coxph(Surv(followup_date,followup_rev) ~ implant + age + sex, data = Tot_pop)


sum1 <- summary(cox1)
sum1$coefficients
sum1$coefficients[1,5]
                  