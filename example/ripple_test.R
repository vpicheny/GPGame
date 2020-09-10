## Tests on the switching ripple design problem

#' @title Switching ripple test function
#' Correspond a switching ripple suppressor design problem for voltage source inversion in powered system.
#' @param x vector specifying the location where the function is to be evaluated, of size k + 4, k > 1, see Details.
#' @param constrained optional boolean, with \code{TRUE} the last 5 columns correspond to the values of the constraints
#' @details 
#' Columns of \code{x} correspond to L1, L2, L3, C1, ..., Ck, Cf where k is an arbitrary integer > 1.
#' 
#' Parameters of the problems follow Table 2 in (Zhang et al, 2019).
#' 
#' @return Vector of objective values. The first k values are related to the suppression of harmonics while the k+1 one is the sum of inductors.
#' 
#' 
#' @references
#' Zhang, Z., He, C., Ye, J., Xu, J., & Pan, L. (2019). Switching ripple suppressor design of the grid-connected inverters: A perspective of many-objective optimization with constraints handling. Swarm and evolutionary computation, 44, 293-303. \cr
#' 
#' He, C., Tian, Y., Wang, H., & Jin, Y. (2019). A repository of real-world datasets for data-driven evolutionary multiobjective optimization. Complex & Intelligent Systems, 1-9.
switch_ripple_n <- function(x, constrained = FALSE){
  
  d <- length(x) # number of variables
  nr <- length(x) - 4 # n in the paper
  
  ## Define parameters
  Ell <- 400
  Prated <- 65000
  Vdc <- 800
  omega0 <- 100 * pi
  omegasw <- 32000 * pi
  Iref <- 141
  Rlk <- rep(0.005, nr)
  Rb <- Ell^2/Prated
  
  L1 <- x[1]
  L2 <- x[2]
  L3 <- x[3]
  Ck <- x[4:(d-1)]
  Cf <- x[d]
  
  Lk <- 1 / (Ck * (1:nr * omegasw)^2)
  
  GLC <- function(s){
    res <- (L2*s * (L3 * Cf * s^2 + 1) + L3 * s) / (Lk * s + Rlk + 1 / (Ck * s))
    sum(res)
  }
  
  G <- function(s){
    1 / (L1 * s * GLC(s) + L2 * s * (L3 * Cf * s^2 + 1) + L3 * s)
  }
  
  fi <- function(i) 20 * log(Mod(G(i * omegasw * 1i)))
  
  fs <- sapply(1:nr, fi)
  fM <- L1 + L2 + L3 + sum(Lk)
  
  if(constrained){
    g1 <- sum(Ck) + Cf - 0.05/(Rb * omega0)
    g2 <- L1 + L2 + L3 - 0.1 * Rb / omega0
    # g31 <- 0.2 - 2 * pi * Vdc / (8 * L1 * omegasw * Iref) # inactive constraints
    # g32 <- 2 * pi* Vdc / (8 * L1 * omegasw * Iref) - 0.4
    g4 <- L2 + L3 - L1
    g51 <- omegasw/2 - sqrt((L1 + L2 + L3) / (L1 * (L2 + L3) * (sum(Ck) + Cf)))
    g52 <- sqrt((L1 + L2 + L3)/(L1 * (L2 + L3) * (sum(Ck) + Cf))) - 3*omegasw / 4
    return(c(fs, fM, g1, g2, g4, g51, g52))
  }else{
    return(c(fs, fM))
  }
  
}


################################################################################
## Try game solving
library(GPGame)


set.seed(42)
d <- 8
nobj <- d-3 # unconstrained
lower <- c(0.000110815602836879, 7.83532027529331e-06, 1.29313391262165e-06, rep(1.29313391262165e-06, d-3))
upper <- c(0.000221631205673759, 0.000783532027529331, 0.000783532027529331, rep(6.46566956310825e-05, d-3))


fn <- function(x){
  res <- switch_ripple_n(x * (upper - lower) + lower)
  # ll <- length(res)
  # return(res[1:(ll-5)])
}

## Nash version

gridtype <- "lhs"
n.init <- 50
n.ite <- 50
n.s <- c(26, rep(11, nobj - 2), 51)

x.to.obj <- c(nobj, nobj, nobj, 1:(nobj - 1), 1)
filtercontrol <- list(nsimPoints=1000, ncandPoints=200,
                      filter=c("window", "Pnash"))

opt <- solve_game(fun = fn, nobj = nobj, d = d, n.init = n.init, n.ite = n.ite, x.to.obj = x.to.obj,
                  integcontrol=list(n.s=n.s, gridtype=gridtype), kmcontrol = list(model.trend=~.), 
                  filtercontrol=filtercontrol)

if(FALSE) pdf("RippleNash.pdf", width = 8, height = 8)
pairs(rbind(cbind(opt$model[[1]]@y, opt$model[[2]]@y, opt$model[[3]]@y, opt$model[[4]]@y, opt$model[[5]]@y),
            opt$Eq.poff), col = c(rep(1, n.init), rep(3, n.ite), rep(2, nrow(opt$Eq.design))), pch = c(rep(3, n.init), rep(4, n.ite), rep(20, nrow(opt$Eq.design))),
      labels = c(expression(f[1]), expression(f[2]), expression(f[3]), expression(f[4]), expression(f[5])))
if(FALSE) dev.off()

## identify best design: 
which(opt$model[[1]]@X[,1] ==opt$Eq.design[1])

plot(opt$model[[2]]@y, opt$model[[3]]@y)
points(matrix(opt$Eq.poff[,c(2,3)], ncol = 2), col = "red", pch = 20)

## Kalai version



