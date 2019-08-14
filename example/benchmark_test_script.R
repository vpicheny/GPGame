library(GPGame)
library(GPareto)
library(plot3D)
library(rgl)
library(parallel)
library(DiceOptim)
# library(Rfast)
library(matrixStats)
set.seed(42)

## For use with a job
args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  print("No argument provided, use S by default")
  config <- "S"
}else{
  print(args)
  config <- args[1] # "S" or "M" or "L"
}

config <- "SMS"

testfun <- "DTLZ2"  # "hartman" "DTLZ2"
# config <- "S" # "S", "M", "L", "XL", "baseline", "RS", "SMS"
pb_type <- "discrete" # "discrete", "continuous"
equilibrium <- "CKS" #"KS", "CKS"
compute_actual <- TRUE

if (testfun == "DTLZ2"){
  directory <- "./Test_results/dtlz2/" # "~/Code/GPGame/example/Test_results/dtlz2/"
  fun <- DTLZ2
  dim <- 5
  nobj <- 4
  n.init <- 10
  n.ite <- 60
} else {
  directory <- "./Test_results/hartman/" # "~/Code/GPGame/example/Test_results/hartman/"
  fun <- function(x){
    I <- matrix(c(rep(0, 6), rep(1, 6), c(0,1,0,0,1,1), c(1,0,1,1,0,0), c(0,0,0,1,1,1), c(1,1,1,0,0,0)), byrow=TRUE, nrow=6)
    J <- matrix(c(1:6, 6:1, c(2,4,6,1,3,5), c(5,3,1,2,4,6), c(3,6,1,4,2,5), c(4,2,6,5,3,1)), byrow=TRUE, nrow=6)
    y <- rep(NA, 6)
    for (i in 1:6) {
      xj <- x[J[i,]]
      xi <- xj*(.5*I[i,] +.5*(1 - I[i,])) + .5*I[i,]
      y[i] <- -log(-hartman6(xi))
    }
    return(y)
  }
  
  dim <- 6
  nobj <- 6
  n.init <- 12
  n.ite <- 78
}

if (config == "S") {
  nsimPoints <- 250
  ncandPoints <- 125
  n.sim <- 25
  n.ynew <- 10
} else if (config == "M") {
  nsimPoints <- 500
  ncandPoints <- 250
  n.sim <- 50
  n.ynew <- 20
} else if (config == "L") {
  nsimPoints <- 1000
  ncandPoints <- 500
  n.sim <- 100
  n.ynew <- 20
} else {
  nsimPoints <- 2000
  ncandPoints <- 1000
  n.sim <- 200
  n.ynew <- 40
}

config_number <-  paste0(config, "_", pb_type)

ntests <- 10
n.s <- 2e5
n.s.large <- 1e7

dir.create(file.path(directory), showWarnings = FALSE)

exp_name <- paste0(directory, "config_", config_number, "_")

ncores <-  min(20, detectCores())

formals(fun)$nobj <- nobj

if (pb_type == "discrete") {
  integcontrol <- generate_integ_pts(n.s=n.s, d=dim, nobj=nobj, equilibrium="KSE",
                                     gridtype="lhs", lb = rep(0,dim), ub = rep(1,dim))
  integcontrol$renew <- FALSE
  integ.pts <- integcontrol$integ.pts
} else {
  integcontrol_large <- generate_integ_pts(n.s=n.s.large, d=dim, nobj=nobj, equilibrium="KSE",
                                           gridtype="lhs", lb = rep(0,dim), ub = rep(1,dim))
  
  integcontrol <- generate_integ_pts(n.s=n.s, d=dim, nobj=nobj, equilibrium="KSE",
                                     gridtype="subset", lb = rep(0,dim), ub = rep(1,dim),
                                     init_set = integcontrol_large$integ.pts)
  integcontrol$renew <- TRUE
  integ.pts <- integcontrol_large$integ.pts
}

# Set controls
filtercontrol <- list(nsimPoints=nsimPoints, ncandPoints=ncandPoints, 
                      filter=c("window", "window"))
simucontrol <- list(n.sim = n.sim, n.ynew=n.ynew)
kmcontrol <- list(lb=rep(.1, dim), ub = rep(3, dim), model.trend=~1, covtype="matern5_2")

if (compute_actual) {
  # Actual function values, Pareto front and KS
  if (testfun == "DTLZ2"){
    fun.grid <- fun(integ.pts)
  } else {
    fun.grid <- t(apply(integ.pts, 1, fun))
  }
  
  pf <- nonDom(fun.grid)
  Shadow <- apply(pf, 2, min)
  Nadir <- apply(pf, 2, max)
}

if (equilibrium == "KS") {
  if (compute_actual) {
    KS_act <- getEquilibrium(pf, equilibrium = "KSE", nobj = nobj)
    
    save(list = c("ntests", "n.init", "Shadow", "Nadir", "pf", "fun.grid", "integcontrol", "KS_act"),
         file=paste0(exp_name, "config_and_solution.RData"))
  }
  
  # Run solver
  KS_mat <- NULL
  # for(ii in 5){
  for(ii in 1:ntests){  #ntests
    if (config == "baseline") {
      res <- solve_game_baseline(fun, equilibrium = "KSE", crit = "sur", n.init=n.init, n.ite=n.ite,
                                 d = dim, nobj=nobj, x.to.obj = NULL,
                                 integcontrol=integcontrol, simucontrol = simucontrol,
                                 filtercontrol=filtercontrol, kmcontrol=kmcontrol,
                                 ncores = ncores, trace=2, seed=ii, freq.exploit=5)
    } else if (config == "RS") {
      res <- solve_game_baseline(fun, equilibrium = "KSE", crit = "sur", n.init=n.init, n.ite=n.ite,
                                 d = dim, nobj=nobj, x.to.obj = NULL,
                                 integcontrol=integcontrol, simucontrol = simucontrol,
                                 filtercontrol=filtercontrol, kmcontrol=kmcontrol,
                                 ncores = ncores, trace=2, seed=ii, baseline_type="RS")
    } else if(config == "SMS"){
      initialDoE <- DiceDesign::lhsDesign(n.init, dim, seed = ii)$design
      res <- easyGParetoptim(fn = fun, budget = n.ite + n.init, lower = rep(0, dim), upper = rep(1, dim), par = initialDoE,
                             control = list(inneroptim = "discrete", candidate.points = integcontrol$integ.pts), ncores = ncores)
      res$Eq.poff <- getEquilibrium(Z = res$value, nobj = nobj, equilibrium = "KSE", Nadir = Nadir, Shadow = Shadow)

    }else{
      res <- solve_game(fun, equilibrium = "KSE", crit = "sur", n.init=n.init, n.ite=n.ite,
                        d = dim, nobj=nobj, x.to.obj = NULL,
                        integcontrol=integcontrol, simucontrol = simucontrol,
                        filtercontrol=filtercontrol, kmcontrol=kmcontrol,
                        ncores = ncores, trace=2, seed=ii, freq.exploit=5)
    }
    KS_mat <- rbind(KS_mat, res$Eq.poff)
    save(list = "res", file=paste0(exp_name, "KSE_run_", ii, ".RData"))
  } 
} else {
  n.samp <- 1e4 #sampled points for Oakley-style interpolation of simulations
  integcontrol <- c(integcontrol, list(kweights = TRUE, nsamp = n.samp))
  
  if (compute_actual) {
    CKS_act <- getEquilibrium(fun.grid, equilibrium = "CKSE", nobj = nobj)
    
    save(list = c("ntests", "n.init", "Shadow", "Nadir", "pf", "fun.grid", "integcontrol", "CKS_act"),
         file=paste0(exp_name, "config_and_solution.RData"))
  }
  
  CKS_mat <- NULL
  
  for(ii in 1:ntests){
    if (config == "baseline") {
      res <- solve_game_baseline(fun, equilibrium = "CKSE", crit = "sur", n.init=n.init, n.ite=n.ite,
                                 d = dim, nobj=nobj, x.to.obj = NULL,
                                 integcontrol=integcontrol, simucontrol = simucontrol,
                                 filtercontrol=filtercontrol, kmcontrol=kmcontrol,
                                 ncores = ncores, trace=2, seed=ii, freq.exploit=5)
    } else if (config == "RS") {
      res <- solve_game_baseline(fun, equilibrium = "CKSE", crit = "sur", n.init=n.init, n.ite=n.ite,
                                 d = dim, nobj=nobj, x.to.obj = NULL,
                                 integcontrol=integcontrol, simucontrol = simucontrol,
                                 filtercontrol=filtercontrol, kmcontrol=kmcontrol,
                                 ncores = ncores, trace=2, seed=ii, baseline_type="RS")
    } else{
      res <- solve_game(fun, equilibrium = "CKSE", crit = "sur", n.init=n.init, n.ite=n.ite,
                        d = dim, nobj=nobj, x.to.obj = NULL,
                        integcontrol=integcontrol, simucontrol = simucontrol,
                        filtercontrol=filtercontrol, kmcontrol=kmcontrol,
                        ncores = ncores, trace=3, seed=ii, freq.exploit=5)
    }
    
    CKS_mat <- rbind(CKS_mat, res$Eq.poff)
    save(list = "res", file=paste0(exp_name, "CKSE_run_", ii, ".RData"))
  }
}
