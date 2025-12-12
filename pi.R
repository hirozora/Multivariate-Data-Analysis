# For fun only. 
# See: https://academo.org/demos/estimating-pi-monte-carlo/

Nreps <- 1000000
x     <- runif(Nreps, -1, 1)
y     <- runif(Nreps, -1, 1)
mean(x^2 + y^2 < 1) * 4

# We can also try using the formula of the volume of a sphere.
# For radius 1, the volume is (4/3)*pi.
# So the following should also work to approx. pi:
x     <- runif(Nreps, -1, 1)
y     <- runif(Nreps, -1, 1)
z     <- runif(Nreps, -1, 1)
mean(x^2 + y^2 + z^2 < 1) * 6
