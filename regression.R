mydata     <- read.csv("data.csv")
reg.output <- lm(y ~ x, mydata)

# Summary:
summary(reg.output)
# I can see that R-squared is 77.1%, and the regression coefficient of x is 
#   significantly different from 0 at alpha = 5%.
  
