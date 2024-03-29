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
library(parallel)

parallel <- FALSE #TRUE #
if(parallel) ncores <- 5 else ncores <- 1

d <- 8
nobj <- d - 3 # unconstrained
lower <- c(0.000110815602836879, 7.83532027529331e-06, 1.29313391262165e-06, rep(1.29313391262165e-06, d-3))
upper <- c(0.000221631205673759, 0.000783532027529331, 0.000783532027529331, rep(6.46566956310825e-05, d-3))


fn <- function(x){
  res <- switch_ripple_n(x * (upper - lower) + lower)
  # ll <- length(res)
  # return(res[1:(ll-5)])
}



## Nash version

gridtype <- "lhs"
n.init <- 80
n.ite <- 70
n.s <- c(26, rep(11, nobj - 2), 51)

x.to.obj <- c(nobj, nobj, nobj, 1:(nobj - 1), 1)
filtercontrol <- list(nsimPoints=1000, ncandPoints=200,
                      filter=c("window", "Pnash"))

set.seed(42)
optNE <- solve_game(fun = fn, nobj = nobj, d = d, n.init = n.init, n.ite = n.ite, x.to.obj = x.to.obj, equilibrium = "NE",
                  integcontrol=list(n.s=n.s, gridtype=gridtype), kmcontrol = list(model.trend=~.), 
                  filtercontrol=filtercontrol, ncores = ncores)

if(FALSE) pdf("RippleNash.pdf", width = 8, height = 8)
pairs(rbind(cbind(optNE$model[[1]]@y, optNE$model[[2]]@y, optNE$model[[3]]@y, optNE$model[[4]]@y, optNE$model[[5]]@y),
            optNE$Eq.poff), col = c(rep(1, n.init), rep(3, n.ite), 2), pch = c(rep(3, n.init), rep(4, n.ite), 20),
      labels = c(expression(f[1]), expression(f[2]), expression(f[3]), expression(f[4]), expression(f[5])))
if(FALSE) dev.off()
save.image(paste0("RippleNE_", n.init, "_", n.ite, "_8d.RData"))

## identify best design: 
which(optNE$model[[1]]@X[,1] == optNE$Eq.design[1])

plot(optNE$model[[2]]@y, optNE$model[[3]]@y)
points(matrix(optNE$Eq.poff[,c(2,3)], ncol = 2), col = "red", pch = 20)

## Kalai version
set.seed(42)
gridtype <- "lhs"
n.init <- 80
n.ite <- 70
n.s2 <- 1e6
filtercontrol <- list(nsimPoints=1000, ncandPoints=200, filter=c("window", "window"))
optKS <- solve_game(fun = fn, nobj = nobj, d = d, n.init = n.init, n.ite = n.ite, x.to.obj = NULL, equilibrium = "KSE",
                  integcontrol=list(n.s=n.s2, gridtype=gridtype), kmcontrol = list(model.trend=~.), 
                  filtercontrol=filtercontrol, ncores = ncores)
save.image(paste0("RippleKS_", n.init, "_", n.ite, "_8d.RData"))

if(FALSE) pdf("RippleKS.pdf", width = 8, height = 8)
Y_KS <- cbind(optKS$model[[1]]@y, optKS$model[[2]]@y, optKS$model[[3]]@y, optKS$model[[4]]@y, optKS$model[[5]]@y)
isdPF_KS <- which(emoa::is_dominated(t(Y_KS)))
cols <- c(rep(1, n.init), rep(3, n.ite), 1, rep(NA, 150))
pchs <- c(rep(3, n.init), rep(4, n.ite), 24, rep(NA, 150))
bgs <- c(rep(NA, n.init+n.ite), 2, rep(NA, 150))
pchs[isdPF_KS] <- 1
cexs <- c(rep(1.5, nrow(Y_KS) + 1), rep(NA, 150))
cexs[isdPF_KS] <- 1
pairs(rbind(Y_KS, optKS$Eq.poff, Y_NKS), col = cols, pch = pchs, bg = bgs, cex = cexs,
      labels = c(expression(f[1]), expression(f[2]), expression(f[3]), expression(f[4]), expression(f[5])), main = "KS")
if(FALSE) dev.off()



## Nash Kalai version
set.seed(42)
gridtype <- "lhs"
n.init <- 80
n.ite <- 70
n.s <- c(26, rep(11, nobj - 2), 51)

x.to.obj <- c(nobj, nobj, nobj, 1:(nobj - 1), 1)
filtercontrol <- list(nsimPoints=1000, ncandPoints=200,
                      filter=c("window", "Pnash"))

optNKS <- solve_game(fun = fn, nobj = nobj, d = d, n.init = n.init, n.ite = n.ite, x.to.obj = x.to.obj, equilibrium = "NKSE",
                  integcontrol=list(n.s=n.s, gridtype=gridtype), kmcontrol = list(model.trend=~.), 
                  filtercontrol=filtercontrol, ncores = ncores)

save.image(paste0("RippleNKS_", n.init, "_", n.ite, "_8d.RData"))

# Get poff from Nash search: 
Nash_poff <- optNE$Eq.poff

if(FALSE) pdf("RippleNKS.pdf", width = 8, height = 8)
Y_NKS <- cbind(optNKS$model[[1]]@y, optNKS$model[[2]]@y, optNKS$model[[3]]@y, optNKS$model[[4]]@y, optNKS$model[[5]]@y)
isdPF_NKS <- which(emoa::is_dominated(t(Y_NKS)))
cols <- c(rep(1, n.init), rep(3, n.ite), 1, 1, rep(NA, 150))
pchs <- c(rep(3, n.init), rep(4, n.ite), 23, 25, rep(NA, 150))
bgs <- c(rep(NA, n.init+n.ite), 2, 4, rep(NA, 150))
pchs[isdPF_NKS] <- 1
cexs <- c(rep(1, nrow(Y_NKS) + 2), rep(NA, 150))
cexs[isdPF_NKS] <- 1
pairs(rbind(Y_NKS, optNKS$Eq.poff,  Nash_poff,Y_KS), col = cols, pch = pchs, cex = cexs, bg = bgs,
      labels = c(expression(f[1]), expression(f[2]), expression(f[3]), expression(f[4]), expression(f[5])), main = "NKS")
if(FALSE) dev.off()

################################################################################
## Added at MACODA revision

load("./Ripple_results/RippleKS_80_70_8d.RData")
load("./Ripple_results/RippleNE_80_70_8d.RData")
load("./Ripple_results/RippleNKS_80_70_8d.RData")


# Get poff from Nash search: 
Nash_poff <- optNE$Eq.poff

# New version
if(FALSE) pdf("RippleNKS2.pdf", width = 8, height = 8)
Y_NKS <- cbind(optNKS$model[[1]]@y, optNKS$model[[2]]@y, optNKS$model[[3]]@y, optNKS$model[[4]]@y, optNKS$model[[5]]@y)
isdPF_NKS <- which(emoa::is_dominated(t(Y_NKS)))
cols <- c(rep(1, n.init), rep(3, n.ite), 1, 1, 1, rep(NA, 150))
pchs <- c(rep(3, n.init), rep(4, n.ite), 23, 22, 25,rep(NA, 150))
bgs <- c(rep(NA, n.init+n.ite), 4, 7, 2, rep(NA, 150))
pchs[isdPF_NKS] <- 1
cexs <- c(rep(1.5, nrow(Y_NKS) + 3), rep(NA, 150))
cexs[isdPF_NKS] <- 1
cexs[151:153] <- 1.5
pairs(rbind(Y_NKS, Nash_poff, optKS$Eq.poff, optNKS$Eq.poff, Y_KS), col = cols, pch = pchs, cex = cexs, bg = bgs,
      labels = c(expression(f[1]), expression(f[2]), expression(f[3]), expression(f[4]), expression(f[5])), main = "NKS")
if(FALSE) dev.off()


if(export) pdf("parcoords_NKS.pdf", width = 12, height = 6)
cols2 <- cols
cols2[151:152] <- c(2,4)
ltys2 <- c(rep(1,150), 1,2, rep(5, 150))
ltys2[isdPF_NKS] <- 3
lwds2 <- c(rep(1,150), 3, 3, rep(0,150))
dataNKS <- rbind(Y_NKS, optNKS$Eq.poff,  Nash_poff,Y_KS)
colnames(dataNKS) <- c("f1", "f2", "f3", "f4", "f5")
MASS::parcoord(dataNKS, col=cols2, lty = ltys2, lwd = lwds2)
if(export) dev.off()


if(FALSE) pdf("RippleKS2.pdf", width = 8, height = 8)
Y_KS <- cbind(optKS$model[[1]]@y, optKS$model[[2]]@y, optKS$model[[3]]@y, optKS$model[[4]]@y, optKS$model[[5]]@y)
isdPF_KS <- which(emoa::is_dominated(t(Y_KS)))
cols <- c(rep(1, n.init), rep(3, n.ite), 1, 1, 1, rep(NA, 150))
pchs <- c(rep(3, n.init), rep(4, n.ite), 23, 22, 25, rep(NA, 150))
bgs <- c(rep(NA, n.init+n.ite), 4, 7, 2, rep(NA, 150))
pchs[isdPF_KS] <- 1
cexs <- c(rep(1.5, nrow(Y_KS) + 3), rep(NA, 150))
cexs[isdPF_KS] <- 1
pairs(rbind(Y_KS, Nash_poff, optKS$Eq.poff, optNKS$Eq.poff, Y_NKS), col = cols, pch = pchs, bg = bgs, cex = cexs,
      labels = c(expression(f[1]), expression(f[2]), expression(f[3]), expression(f[4]), expression(f[5])), main = "KS")
if(FALSE) dev.off()

if(export) pdf("parcoords_KS.pdf", width = 12, height = 6)
cols2 <- cols
cols2[151] <- c(7)
ltys2 <- c(rep(1,150), 1, rep(5, 150))
ltys2[isdPF_KS] <- 3
lwds2 <- c(rep(1,150), 3, rep(0,150))
dataKS <- rbind(Y_KS, optKS$Eq.poff, Y_NKS)
colnames(dataKS) <- c("f1", "f2", "f3", "f4", "f5")
MASS::parcoord(dataKS, col=cols2, lty = ltys2, lwd = lwds2)
if(export) dev.off()

# If needed: empirical comparison
empPF <- mco::nsga2(fn = fn, idim = d, odim = 5, lower.bounds = rep(0,d), upper.bounds = rep(1,d), popsize = 5000, generations = 300)
pairs(rbind(empPF$value, Y_KS,Y_NKS), col = c(rep(1,1000), rep(2,300)))
