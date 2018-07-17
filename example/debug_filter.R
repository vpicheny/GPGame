## Configuration exhibiting an issue with the filters.
## Should return several warning messages of the form:
##
## 1: In 1:n.s.target :
##    numerical expression has 2 elements: only the first used
##
## And the nb of integration / candidate pts should not match the targets

# Define objective function (R^2 -> R^2)
fun1 <- function(x, ...)
{
  if (is.null(dim(x)))    x <- matrix(x, nrow = 1)
  b1 <- 15 * x[, 1] - 5
  b2 <- 15 * x[, 2]
  return(cbind((b2 - 5.1*(b1/(2*pi))^2 + 5/pi*b1 - 6)^2 + 10*((1 - 1/(8*pi)) * cos(b1) + 1),
               -sqrt((10.5 - b1)*(b1 + 5.5)*(b2 + 0.5)) - 1/30*(b2 - 5.1*(b1/(2*pi))^2 - 6)^2-
                 1/3 * ((1 - 1/(8 * pi)) * cos(b1) + 1)))
}

# To use parallel computation (turn off on Windows)
library(parallel)
parallel <- TRUE #
if(parallel) ncores <- detectCores() else ncores <- 1

# Simple configuration: no filter, discretization is a 21x21 grid

# Grid definition
n.s <- rep(21, 2)
x.to.obj   <- c(1,2)
gridtype <- 'cartesian'
n.ite <- 2

filtercontrol <- list(nsimPoints=200, ncandPoints=100,
                      filter=c("window", "window"))

res <- solve_game(fun1, equilibrium = "KSE", crit = "sur", n.init=6, n.ite=n.ite,
                  d = 2, nobj=2, x.to.obj = x.to.obj,
                  integcontrol=list(n.s=n.s, gridtype=gridtype),
                  filtercontrol=filtercontrol,
                  ncores = ncores, trace=1, seed=1)

# Get estimated equilibrium and corresponding pay-off
NE <- res$Eq.design
Poff <- res$Eq.poff

# Draw results
plotGame(res, equilibrium = "KSE", Nadir=Nadir, calibcontrol=calibcontrol) #Nadir=c(Inf, -20))
