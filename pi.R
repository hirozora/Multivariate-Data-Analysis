# For fun only. 
# See: https://academo.org/demos/estimating-pi-monte-carlo/

Nreps <- 1000000
x     <- runif(Nreps, -1, 1)
y     <- runif(Nreps, -1, 1)
mean(x^2 + y^2 < 1)
