context('Equilibrium computation')
library(GPareto)

test_that("KS equilibrium computation is correct", {
  
  ## 2D
  xgrid <- seq(0,1,length.out = 101)
  Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
  Z <- P2(Xgrid)
  
  PF <- nonDom(Z)
  Shadow <- apply(PF, 2, min)
  Nadir <- apply(PF, 2, max)
  KSeq <- GPGame:::getKS(Z = Z, Nadir = Nadir, Shadow = Shadow)
  
  ratios <- (matrix(Nadir, nrow = nrow(PF), ncol = ncol(PF), byrow = TRUE) - PF) %*% diag(1/(Nadir - Shadow))
  KS_ref <- PF[which.max(apply(ratios, 1, min)),, drop = FALSE]
  
  expect_equal(KSeq$KS, KS_ref)
  
  # plot(Z)
  # plotParetoEmp(PF, col ='red')
  # points(KSeq$KS, col = "green", pch = 20)
  
  ## 5D (Pareto front of a ball)
  nvar <- 5
  n <- 1e6
  Z <- matrix(rnorm(n * nvar), n)
  Z <- t(apply(Z, 1, function(x) return(-abs(x)/sqrt(sum(x^2)))))
  Shadow <- rep(-1, nvar)
  Nadir <- rep(0, nvar)
  KS_ref <- matrix(rep(- 1/sqrt(nvar), nvar), nrow = 1)
  KSeq <- GPGame:::getKS(Z = Z, Nadir = Nadir, Shadow = Shadow)
  
  expect_equal(KSeq$KS, KS_ref, tol = 1e-2)
})