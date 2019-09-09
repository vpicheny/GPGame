#------------------------------------------------
# Transform observations in the calibration case
transform_obs <- function(observations, calibcontrol) {
  if (!is.null(calibcontrol$target)) {
    observations <- (observations - matrix(rep(calibcontrol$target, nrow(observations)), byrow=TRUE, nrow=nrow(observations)))^2
    if (calibcontrol$log) observations <- log(observations + calibcontrol$offset)
  }
  return(observations)
}

#------------------------------------------------
# Aggregate sets of observations
aggregate_obs <- function(model_list, n.init=0) {
  if (n.init > 0) {
    observations <- Reduce(cbind, lapply(model_list[[1]], slot, "y"))[1:n.init,]
    X <- model_list[[1]][[1]]@X[1:n.init,]
    noise.var <- Reduce(cbind, lapply(model_list[[1]], slot, "noise.var"))[1:n.init,]
  } else {
    X <- noise.var <- observations <- c()
  }
  
  for (i in 1:length(model_list)){
    obs <- Reduce(cbind, lapply(model_list[[i]], slot, "y"))
    x <- model_list[[i]][[1]]@X
    noise <- Reduce(cbind, lapply(model_list[[i]], slot, "noise.var"))
    I <- (n.init+1):model_list[[i]][[1]]@n
    X <- rbind(X, x[I,])
    noise.var <- rbind(noise.var, noise[I,])
    observations <- rbind(observations, obs[I,])
  }
  
  return(list(X=X, observations=observations, noise.var=noise.var))
}

#------------------------------------------------
# Build GP models
build_models <- function(X, observations, noise.var, ncores=1) {
  d <- ncol(X)
  nobj <- ncol(observations)
  my.km <- function(i) {
    km(~., design = X, response = observations[, i], optim.method="BFGS",
       noise.var=noise.var[,i], multistart=5,
       control=list(trace=FALSE), lower=rep(.1,d), upper=rep(5,d))
  }
  model <- mclapply(1:nobj, my.km, mc.cores=ncores)
  return(model)
}

#------------------------------------------------
# Get Nadir and Shadow
get_nadir_and_shadow <- function(observations, model, n.s, nsimPoints, calibcontrol, Nadir=NULL, Shadow=NULL, nsim=1, ncores=1) {
  nobj <- ncol(observations)
  
  if (is.null(Nadir)) Nadir <- rep(Inf, nobj)
  if (is.null(Shadow)) Shadow <- rep(-Inf, nobj)
  
  # First create a new set of integration points
  integcontrol <- list(gridtype='lhs', renew=TRUE, n.s=n.s)
  res <- generate_integ_pts(n.s=n.s, d=model[[1]]@d, nobj=nobj, x.to.obj=NULL, equilibrium="KSE",
                            gridtype="lhs", include.obs=FALSE)
  integcontrol$integ.pts <- integ.pts <- res$integ.pts

  # Compute prediction at those points  

  pred <- mclapply(model, FUN=predict, newdata = integ.pts, checkNames = FALSE, type = "UK", light.return = TRUE, mc.cores=ncores)
  pred <- lapply(pred, function(alist) list(mean=alist$mean, sd=alist$sd))
  predmean <- Reduce(cbind, lapply(pred, function(alist) alist$mean))
  
  # Retain only the most interesting ones for nadir and shadow
  integ_pts_ns <- get_nadir_and_shadow_integpoints(model, pred, integ.pts, Nadir, Shadow, nsimPoints, nobj, target=calibcontrol$target)
  
  if (nsim == 1) {
    # If nsim is one: compute the nadir and shadow of the GP means
    
    # Compute predictions at new points and concatenate
    pred_ns <- predict_kms(model=model, newdata = integ_pts_ns, checkNames = FALSE, type = "UK",
                           light.return = TRUE, mc.cores=ncores)
    predmean <- rbind(predmean, t(pred_ns$mean), observations)
    
    # Calibration mode
    predmean <- transform_obs(pred_mean, calibcontrol)
    
    # Get non-dominated predictions
    PFmean <- nonDom(predmean)
    
    # Get Shadow and Nadir
    Shadow <- apply(PFmean, 2, min)
    Nadir <- apply(PFmean, 2, max)
    
  } else {
    # Conditional simulations at new points
    Z <- try(t(Reduce(rbind, mclapply(model, simulate, nsim=nsim, newdata=integ_pts_ns, 
                                                 cond=TRUE, checkNames=FALSE, nugget.sim = 10^-6, mc.cores=ncores))))

    # Compute N and S for each sample
    Nadir <- Shadow <- matrix(NA, nsim, nobj)
    for (u in 1:nsim) {
      # Retain only the non-dominated points of the sample u
      J <- seq(u, ncol(Z), nsim)
      ZZ <- rbind(Z[,J, drop = FALSE], observations)
      
      # Calibration mode
      ZZ <- transform_obs(ZZ, calibcontrol)
      Zred <- nonDom(ZZ)
      
      # Get Shadow and Nadir
      Shadow[u, ] <- apply(Zred, 2, min)
      Nadir[u, ]  <- apply(Zred, 2, max)
    }
  }
  return(list(Nadir=Nadir, Shadow=Shadow))
}

#---------------------------------
get_nadir_and_shadow_integpoints <- function(model, predictions, integ.pts, Nadir, Shadow, nsimPoints, nobj, target) {
  source("R/prob.of.non.domination.R")
  ## Add potential minimas of each objective based on EI (utopia), and EI x Pdom (nadir)
  IKS <- NULL # points of interest for KS (i.e., related to KS and Nadir)
  
  # Number of points to keep:
  n_each <- round(nsimPoints / nobj / 2)
  
  if (max(Nadir) == Inf) {
    PNDom <- prob.of.non.domination(model=model, integration.points=integ.pts, predictions=predictions, nsamp=NULL, target=target)
  }
  if (is.null(target)) {
    PFobs <- nonDom(Reduce(cbind, lapply(model, slot, "y")))
    for (jj in 1:length(model)) {
      discard <- which(predictions[[jj]]$sd/sqrt(model[[jj]]@covariance@sd2) < 1e-06)
      if (Shadow[jj] == -Inf) {
        # EI(min) on each objective (to find potential Utopia points)
        xcr <-  (min(model[[jj]]@y) - predictions[[jj]]$mean)/predictions[[jj]]$sd
        test <- (min(model[[jj]]@y) - predictions[[jj]]$mean)*pnorm(xcr) + predictions[[jj]]$sd * dnorm(xcr)
        test[discard] <- NA
        I <- sort(-test, index.return=TRUE)$ix
        IKS <- c(IKS, I[1:n_each])
      }
      # EI(max) x Pdom on each objective (to find potential Nadir points) unless Nadir is provided
      if (Nadir[jj] == Inf) {
        xcr <-  -(max(PFobs[,jj]) - predictions[[jj]]$mean)/predictions[[jj]]$sd
        test <- (-max(PFobs[,jj]) + predictions[[jj]]$mean)*pnorm(xcr) + predictions[[jj]]$sd * dnorm(xcr)
        test[discard] <- NA
        test <- test * PNDom
        I <- sort(-test, index.return=TRUE)$ix
        IKS <- c(IKS, I[1:n_each])
      }
    }
  } else {
    # Calibration mode
    observations <- Reduce(cbind, lapply(model, slot, "y"))
    PFobs <- nonDom((observations - matrix(rep(target, nrow(observations)), byrow=TRUE, nrow=nrow(observations)))^2)
    for (jj in 1:length(model)) {
      discard <- which(predictions[[jj]]$sd/sqrt(model[[jj]]@covariance@sd2) < 1e-06)
      
      # EI(min) on each objective (to find potential Utopia points)
      zmin <- min((model[[jj]]@y - target[jj])^2)
      mu <- predictions[[jj]]$mean - target[jj]
      sigma   <- predictions[[jj]]$sd
      a2 <- (sqrt(zmin) - mu)/sigma
      a1 <- (-sqrt(zmin) - mu)/sigma
      da2 <- dnorm(a2);   da1 <- dnorm(a1)
      pa2 <- pnorm(a2);   pa1 <- pnorm(a1)
      test <- (zmin - sigma^2 - mu^2) * (pa2 - pa1) + sigma^2 * (a2*da2 - a1*da1) + 2*mu*sigma*(da2 - da1)
      test[discard] <- NA
      I <- sort(-test, index.return=TRUE)$ix
      IKS <- c(IKS, I[1:n_each])
      
      # EI(max) x Pdom on each objective (to find potential Nadir points) unless Nadir is provided
      if (Nadir[jj] == Inf) {
        zmax <- max(PFobs[,jj])
        b2 <- (sqrt(zmax) - mu)/sigma
        b1 <- (-sqrt(zmax) - mu)/sigma
        db2 <- dnorm(b2);   db1 <- dnorm(b1)
        pb2 <- pnorm(b2);   pb1 <- pnorm(b1)
        test <- (sigma^2 + mu^2 - zmax) * (1 - pb2 + pb1) + sigma^2 * (b2*db2 - b1*db1) + 2*mu*sigma*(db2 - db1)
        test[discard] <- NA
        test <- test * PNDom
        I <- sort(-test, index.return=TRUE)$ix
        IKS <- c(IKS, I[1:n_each])
      }
    }
  }
  
  integ.pts <- integ.pts[unique(IKS),]
  return(integ.pts)
}

#------------------------------------------------
# Get reference set for copula calculation
get_reference_set <- function(model, n.s, nsimPoints=NULL, calibcontrol, nsim=1, ncores=1) {
  nobj <- length(model)
  
  # First create a new set of integration points
  integcontrol <- list(gridtype='lhs', renew=TRUE, n.s=n.s)
  res <- generate_integ_pts(n.s=n.s, d=model[[1]]@d, nobj=nobj, x.to.obj=NULL, equilibrium="KSE",
                            gridtype="lhs", include.obs=FALSE)
  integcontrol$integ.pts <- integ.pts <- res$integ.pts
  
  # Compute prediction at those points  
  pred <- mclapply(model, FUN=predict, newdata = integ.pts, checkNames = FALSE, type = "UK", light.return = TRUE, mc.cores=ncores)
  pred <- lapply(pred, function(alist) list(mean=alist$mean, sd=alist$sd))
  
  if (nsim ==1) {
    
    predmean <- Reduce(cbind, lapply(pred, function(alist) alist$mean))
    predmean <- transform_obs(pred_mean, calibcontrol)
    
    return(list(X=integ.pts, observations=predmean))
    
  } else {
    # Use the pseudo-observation trick from Oakley
    
    # Gather points with highest variance
    predvar <- apply(Reduce(cbind, lapply(pred, function(alist) log(alist$sd^2))), 1, mean)
    I <- sort(predvar, index.return=TRUE, decreasing=TRUE)$ix[1:nsimPoints]
    my.integ.pts <- integ.pts[I, ]
    
    # Compute pseudo points and kweights
    Xsamp <- matrix(runif(d * n.s), n.s)
    kweights <- list()
    for (u in 1:nobj) {
      kn <- covMat1Mat2(model[[u]]@covariance, Xsamp, my.integ.pts)
      Knn <- try(chol2inv(chol(covMat1Mat2(model[[u]]@covariance, my.integ.pts, my.integ.pts) + diag(1e-6, nrow(my.integ.pts)))))
      kweights <- c(kweights, list(list(kn = kn, Knn = Knn)))
    }

    # Conditional simulations at those points
    Simu <- try(t(Reduce(rbind, mclapply(model, simulate, nsim=nsim, newdata=my.integ.pts, 
                                                 cond=TRUE, checkNames=FALSE, nugget.sim = 10^-6, mc.cores=ncores))))

    # Pseudo-random trajectories using kweights
    ZZ <- list()
    for (u in 1:nsim) {
      J <- seq(u, ncol(Simu), nsim)
      Z <- matrix(NA, nrow = nrow(kweights[[1]]$kn), ncol = length(J))
      for(jj in 1:length(J)) {
        Z[,jj] <- kweights[[jj]]$kn %*% (kweights[[jj]]$Knn %*% Simu[,J[jj]])
      }
      Z <- transform_obs(Z, calibcontrol)
      ZZ <- c(ZZ, list(Z))
    }
    return(list(X=Xsamp, observations=ZZ))
  }
}

#------------------------------------------------
# Transforms observations into rank with respect to a set 
# of reference observations
convert_to_ranks <- function(observations, observations_ref) {
  
  obsrank <- matrix(NA, nrow(observations), ncol(observations))
  
  nrow_ref <- nrow(observations_ref)
  
  for (i in 1:nrow(observations)) {
    diff <- matrix(rep(observations[i,], nrow_ref), nrow=nrow_ref, byrow=TRUE) - observations_ref
    obsrank[i,] <- apply(diff>0, 2, sum)
  }
  return(obsrank / nrow(observations_ref))
}

#------------------------------------------------
# Centrality plot
plot_centrality <- function(observations, observations_ref=NULL, Nadir=NULL, Shadow=NULL, n.init=NULL, add=FALSE, copula=FALSE){
  if (copula) {
    nobj <- ncol(observations)
    if (is.null(observations_ref)) observations_ref <- list(observations[1:n.init,])
    obs <- list()
    for (i in 1:length(observations_ref)) {
      obs <- c(obs, list(convert_to_ranks(observations, observations_ref[[i]])))
    }
    observations <- obs
    # observations <- convert_to_ranks(observations, observations_ref)/n.init
    Nadir <- rep(1, nobj)
    Shadow <- rep(0, nobj)
    
  }
  if (is.null(dim(Nadir))) {
    if (length(observations_ref) == 1) {
      ratios <- sweep(sweep(observations, 2, Nadir, "-"), 2, Shadow - Nadir, "/")
      dks <- apply(ratios, 1, min)
      DKS <- apply(ratios, 1, max)
      
      plot(x, dks, xlab="Iteration #", ylab="Marginal gain", ylim=range(c(dks, DKS)))
      arrows(x, dks, x, DKS, length=0.05, angle=90, code=3)
    } else {
      dks <- dks2 <- matrix(NA, length(observations_ref), nrow(observations[[1]]))
      for (i in 1:length(observations_ref)) {
        ratios <- sweep(sweep(observations[[i]], 2, Nadir, "-"), 2, Shadow - Nadir, "/")
        dks[i,] <- apply(ratios, 1, min)
        dks2[i,] <- apply(ratios, 1, mean)
      }
      x = 1:nrow(ratios)
      plot(x, apply(dks, 2, mean), xlab="Iteration #", ylab="Marginal gain", ylim=range(c(dks, dks2)))
      arrows(x, apply(dks, 2, min), x, apply(dks, 2, max), length=0.05, angle=90, code=3)
      dks <- apply(dks, 2, mean)
      dks2 <- apply(dks2, 2, mean)
    }
  } else {
    dks <- dks2 <- matrix(NA, nrow(Nadir), nrow(observations))
    for (i in 1:nrow(Nadir)) {
      ratios <- sweep(sweep(observations, 2, Nadir[i,], "-"), 2, Shadow[i,] - Nadir[i,], "/")
      dks[i,] <- apply(ratios, 1, min)
      dks2[i,] <- apply(ratios, 1, mean)
    }
    x = 1:nrow(ratios)
    plot(x, apply(dks, 2, mean), xlab="Iteration #", ylab="Marginal gain", ylim=range(c(dks, dks2)))
    arrows(x, apply(dks, 2, min), x, apply(dks, 2, max), length=0.05, angle=90, code=3)
    dks <- apply(dks, 2, mean)
    dks2 <- apply(dks2, 2, mean)
  }
  
  if (!add){
    lines(cummax(dks), col="blue")
    lines(cummax(dks2), col="red")
    if (!is.null(n.init)) abline(v=n.init+.5, lty=2)
  }
}

#------------------------------------------------
# Pairs plot of observations and KS
# Different colors for dominated observations
# Equilibrium found in red
# Actual (if provided) in green
plot_pairs <- function(observations, Eq.poff, n.init, Nadir=NULL, Shadow=NULL){
  if (is.null(nrow(Eq.poff))) Eq.poff <- matrix(Eq.poff, nrow=1)
  nobs <- nrow(observations)
  
  objnames2 <- c("Heat_W","ECS_W","Cook_W","Other_W",
                 "Winter_T","T_cook","T_relax",
                 "T_sleep","T_out")
  # objnames2 <- c("Heat_W","ECS_W","Cook_W","Other_W",
  #                "Winter_T","T_cook","T_relax",
  #                "T_sleep")
  colnames(observations) <- objnames2
  
  n.ite <- nobs - n.init
  ind <- nonDom(observations, return.idx = TRUE)
  colpf <- rep("red", nobs)
  colpf[ind] <- "blue"
  
  colpf <- c(colpf, rep("green", nrow(Eq.poff)), "violet", "violet")
  pch <- c(rep('.', n.init), rep('.', n.ite), rep('.', nrow(Eq.poff)+2))
  cex <- c(rep(4, nobs), 10, 10, 10)
  
  pairs(rbind(observations, Eq.poff, Nadir, Shadow), pch=pch, cex=cex, col=colpf, upper.panel=NULL)
}

#------------------------------------------------
# Plot observations
plot_observations <- function(observations, cols, n.init){
  par(mfrow=c(1,2), mar=c(4,4,4,4))
  plot(NA, xlim=c(0,nobs), ylim=range(observations), xlab="Nb iterations", ylab="logMSE", main="All observed values")
  for (i in 1:ncol(observations)) points(1:nobs, observations[,i], col=cols[i])
  abline(v=n.init, col="black", lty=2)
}

#------------------------------------------------
# Plot all added X
plot_X <- function(observations, cols, n.init, Eq.design){
  par(mfrow=c(3,5), mar=c(2,2,0,1))
  colX <- sample(rainbow(19))
  nobs <- nrow(observations)
  
  for (i in 1:ncol(X)) {
    plot((n.init+1):nobs, X[(n.init+1):nobs,i], col=colX[i], pch=19, xlim=c(n.init, nobs), 
         ylim=c(0,1), ylab="", xlab="", xaxt='n', yaxt='n', ann=FALSE)
    abline(h=Eq.design[i], col=colX[i])
    title(ylab=paste0("X",i), mgp=c(1,1,0), cex.lab=1.2) 
  }
}

#------------------------------------------------
# Plot marginal distribution and KS
plot_X <- function(observations, cols, n.init, objnames, calibcontrol, Eq.poff){
  par(mfrow=c(3,3), mar=c(4,4,4,4))
  for (i in 1:nobj) {
    hist(observations[1:nobs,i], main=objnames[i], xlab="", freq=TRUE)
    abline(v=Eq.poff[i], col=cols[i], lwd=2)
  }
  if (!is.null(calibcontrol$target)) {
    # Calibration mode
    par(mfrow=c(3,3), mar=c(4,4,4,4))
    for (i in 1:nobj) {
      hist(sqrt(exp(observations[1:nobs,i] - calibcontrol$offset))*100, main=objnames[i], xlab="RMSE (%)", freq=TRUE)
      abline(v=sqrt(exp(Eq.poff[i]- calibcontrol$offset))*100, col=cols[i], lwd=2)
    }
  }
}

#------------------------------------------------
# Plot marginal distribution and KS but with a different color 
# for initial points and added ones
plot_histogram_bicolor <- function(observations, n.init, Eq.poff){
  nobs <- nrow(observations)
  nobj <- ncol(observations)
  n.ite <- nobs - n.init
  plots <- list()
  par(mfrow=c(3,3))
  for (i in 1:nobj) {
    dat <- data.frame(value=observations[,i], label=c(rep("init", n.init), rep("optim", n.ite)))
    plots[[i]] <- ggplot(dat, aes(x=value, fill=label)) +  geom_histogram() + 
      ggtitle(colnames(observations)[i]) + geom_vline(xintercept = Eq.poff[i])
  }
  multiplot(plotlist=plots, cols=3)
}

#------------------------------------------------
# Plots everything
plot_bundle <- function(observations, Nadir, Shadow, n.init, Eq.poff, Eq.design, cols, objnames, calibcontrol){
  cols <- rainbow(ncol(observations))
  par(mfrow=c(1, 1))
  plot_centrality(observations, Nadir, Shadow, n.init, add=FALSE)
  plot_pairs(observations, Eq.poff, n.init)
  par(mfrow=c(1, 1))
  plot_observations(observations, cols, n.init)
  plot_X(observations, cols, n.init, objnames, calibcontrol, Eq.poff)
  # plot_histogram_bicolor(observations, n.init, Eq.poff)
}


#------------------------------------------------
# Convergence plots
convergence_plots <- function(models, solution, Nadir=NULL, Shadow=NULL, title=NULL,
                              observations_ref=NULL, copula=FALSE, dist_ratio=TRUE) {
  n_runs <- length(models)
  n_ite <- models[[1]][[1]]@n
  
  # results stores all the minratios, results2 the cumulative min of the minratios
  results <- results2 <- matrix(NA, n_runs, n_ite)
  
  # Compute first the minratio at the solution (if N/S are provided)
  if (!is.null(Nadir)) ratio_sol <- min((solution - Nadir) / (Shadow - Nadir))
  
  for(i in 1:n_runs){
    # model is a list, models is a list of lists
    model <- models[[i]]
    observations <- Reduce(cbind, lapply(model, slot, "y"))
    
    # if CKS, transform observations
    if (copula) {
      observations <- convert_to_ranks(observations, observations_ref)
      Nadir <- rep(1, ncol(observations))
      Shadow <- rep(0, ncol(observations))
      ratio_sol <- min((convert_to_ranks(solution, observations_ref) - Nadir) / (Shadow - Nadir))
    }
    
    if (!dist_ratio) {
      # If dist_ratio=FALSE: convergence is calculated using the L2 norm
      results[i,] <- as.matrix(dist(rbind(solution, observations)))[2:(model[[1]]@n+1),1]
    } else {
      # Regular case: compute the distance to solution in terms of minratio
      ratios <- sweep(sweep(observations, 2, Nadir, "-"), 2, Shadow - Nadir, "/")
      results[i,] <- ratio_sol - apply(ratios, 1, min)
    }
    # results2 for the convergence plots
    results2[i,] <- cummin(results[i,])
  }
  require(ggplot2)
  require(reshape2)
  require(gridExtra)

  df <- data.frame(t(results))
  df$time = 1:n_ite
  df <- melt(df ,  id.vars = 'time', variable.name = 'run')
  
  df2 <- data.frame(t(results2))
  df2$time = 1:n_ite  # to reduce the curves to anything less that n_ite if wanted
  df2 <- melt(df2,  id.vars = 'time', variable.name = 'run')
  
  # p1 <- ggplot(df, aes(time, value)) + geom_line(aes(colour = run))
  p2 <- ggplot(df2, aes(time, value)) + geom_line(aes(colour = run)) + ylim(0, NA) + ggtitle(title)
  grid.arrange(p2, nrow = 1)

  return(list(cummin=results2, alldist=results))
}