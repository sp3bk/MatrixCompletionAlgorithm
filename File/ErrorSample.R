n.per.group<-10

alpha<-0.05 # for a (1.00-alpha)=95% confidence inte

# Simulate raw data for an experiment or observational study.
data.raw <- data.frame(
  treatment=rep(c('A','B'), each=n.per.group),
  value=c(rnorm(n.per.group, 2), rnorm(n.per.group, 3))   
)

value=c(rnorm(n.per.group, 2))
value
data.raw$treatment
# This data frame calculates statistics for each treatment.
data.summary <- data.frame(
  treatment=levels(data.raw$treatment),
  mean=tapply(data.raw$value, data.raw$treatment, mean),
  n=tapply(data.raw$value, data.raw$treatment, length),
  sd=tapply(data.raw$value, data.raw$treatment, sd)
)
data.summary
# Precalculate standard error of the mean (SEM)
data.summary$sem <- data.summary$sd/sqrt(data.summary$n)
# Precalculate margin of error for confidence interval
data.summary$me <- qt(1-alpha/2, df=data.summary$n)*data.summary$sem

g <- list(5,4,3,2,1)
g[[1]] <- append(g[[1]],9)
g[[1]] <- append(g[[1]],8)
g
