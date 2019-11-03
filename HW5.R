library("Sleuth3")
library("ggplot2")

aov_data = ex0524
income_summary = summary(aov_data$Income2005)
qplot(x=IQquartile, y=Income2005, data=aov_data, geom=c("boxplot"))
aov_data$logincome <- log(aov_data$Income2005)
qplot(x=IQquartile, y=logincome, data=aov_data, geom=c("boxplot"))
log_income_summary = summary(aov_data$logincome)

log_data_aov = aov(logincome~IQquartile, data=aov_data)
log_data_anova = anova(log_data_aov)

plot(log_data_aov,which=1)

quartile_means <- with(aov_data,unlist(lapply(split(logincome,IQquartile),mean)))
print(quartile_means)
t_df = qt(0.975, 2580)
print(t_df)
num_participants_in_quartile <- with(aov_data,unlist(lapply(split(logincome,IQquartile),length)))
print(num_participants_in_quartile)
conf_int_plus = (10.785 - 10.050) + (1.961)*(sqrt(0.937))*(sqrt((2)/1292))
print(conf_int_plus)
conf_int_minus = (10.785 - 10.050) - (1.961)*(sqrt(0.937))*(sqrt((2)/1292))
print(conf_int_minus)
pt_estimate = exp(10.785-10.050)
print(pt_estimate)
print(exp(c(0.66, 0.81)))


