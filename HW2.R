library("Sleuth3")
library("ggplot2")


housing <- read.csv("C:/Users/hammera/Downloads/Housing.csv")
head(housing)
View(housing)
summary(housing)

housing$index_sa_log <- log(housing$index_sa)
qplot(region, index_sa_log, data=housing, geom="boxplot")
with(housing, summary(index_sa_log[region=="South"]))
with(housing, summary(index_sa_log[region=="West"]))
qplot(index_sa_log, data=housing, geom="histogram") + facet_grid(region ~.)
new_vec <- subset(housing$index_sa_log, housing$region=="South")



log_t_test <- t.test(index_sa_log~region, data=housing, alternative="two.sided")
print(log_t_test)

back_transform_mean = ((exp(5.532352-5.833817)))
print(back_transform_mean)

exp(c(0.4462554, 0.1566740))








