mydata     <- read.csv("data.csv")
reg.output <- lm(y ~ x, mydata)

# Summary:
summary(reg.output)
# I can see that R-squared is 77.1%, and the regression coefficient of x is 
#   significantly different from 0 at alpha = 5%.

# Plot y against x:
plot(mydata$x, mydata$y)

# Humm, let's try to use Bayesian regression, not sure....
library(rstanarm) # allows fitting Bayesian models, the easy way
bayesreg.output <- stan_glm(y ~ x, family = gaussian(), mydata)
summary(bayesreg.output)
plot(bayesreg.output)

# And from here on I wrote some really bad code, I just want to delete it in 
#  such a way that Git won't complain.

