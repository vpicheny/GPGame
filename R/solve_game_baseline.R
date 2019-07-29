#' @noRd
#' @export
solve_game_baseline <- function(
  fun, ..., equilibrium="NE", crit="sur", model=NULL, n.init=NULL, n.ite, d, nobj, x.to.obj=NULL, noise.var = NULL,
  Nadir=NULL, Shadow=NULL, integcontrol=NULL, simucontrol=NULL, filtercontrol=NULL, kmcontrol=NULL, returncontrol=NULL,
  ncores=1, trace=1, seed=NULL, calibcontrol=NULL, freq.exploit=1e3, baseline_type="1d") {
  
  t1 <- Sys.time()
  set.seed(seed)
  if(ncores > 1) parallel <- TRUE
  
  ####################################################################################################
  #### CHECK INPUTS ##################################################################################
  ####################################################################################################
  # Check integcontrol
  integ.pts <- integcontrol$integ.pts
  expanded.indices <- integcontrol$expanded.indices
  n.s <- integcontrol$n.s
  init_set <- integcontrol$init_set
  gridtype <- switch(1+is.null(integcontrol$gridtype), integcontrol$gridtype, "lhs")
  lb <- switch(1+is.null(integcontrol$lb), integcontrol$lb, rep(0, d))
  ub <- switch(1+is.null(integcontrol$ub), integcontrol$ub, rep(1, d))
  integcontrol$renew <- switch(1+is.null(integcontrol$renew), integcontrol$renew, equilibrium %in% c("KSE", "CKSE"))
  
  # Check simucontrol
  n.ynew <- switch(1+is.null(simucontrol$n.ynew), simucontrol$n.ynew, 10)
  n.sim <- switch(1+is.null(simucontrol$n.sim), simucontrol$n.sim, 10)
  IS <- switch(1+is.null(simucontrol$IS), simucontrol$IS, TRUE)
  cross <- FALSE
  
  # Check kmcontrol
  cov.reestim <- switch(1+is.null(kmcontrol$cov.reestim), kmcontrol$cov.reestim, TRUE) 
  model.trend <- switch(1+is.null(kmcontrol$model.trend), kmcontrol$model.trend, ~1)
  kmlb <- switch(1+is.null(kmcontrol$lb), kmcontrol$lb, rep(.1,d))
  kmub <- switch(1+is.null(kmcontrol$ub), kmcontrol$ub, rep(1,d))
  kmnugget <- switch(1+is.null(kmcontrol$nugget), kmcontrol$nugget, 1e-8)
  control <- switch(1+is.null(kmcontrol$control), kmcontrol$control, list(trace=FALSE))
  covtype <- switch(1+is.null(kmcontrol$covtype), kmcontrol$covtype, "matern5_2")
  
  # Check Noise
  if (!is.null(noise.var)) {
    if (typeof(noise.var) == "closure") noise.var <- match.fun(noise.var)
    else if (typeof(noise.var) == "double" && length(noise.var==1)) noise.var <- rep(noise.var, nobj)
    kmnugget <- NULL
  }
  
  # Check calibcontrol
  if (!is.null(calibcontrol)) {
    if (is.null(calibcontrol$target)){
      cat("calibcontrol should contain a target vector \n")
      return(NA)
    } 
    if (is.null(calibcontrol$log)) calibcontrol$log <- FALSE
    if (is.null(calibcontrol$offset)) calibcontrol$offset <- 0
  }
  
  if (trace>0) cat("--------------------------\n Starting", equilibrium, "\n",
                   "among (", n.s, ") strategies \n --------------------------\n " )
  
  ####################################################################################################
  #### INITIALIZE VARIABLES AND MODELS ###############################################################
  ####################################################################################################
  if (baseline_type == "RS") {
    if (is.null(n.init))     n.init <- 0
    n.init <- n.init + n.ite
    n.ite <- 0
  }
  
  #### Initial design and models ####################
  if (is.null(model)){
    if (is.null(n.init))     n.init <- 5*d
    design <- lhsDesign(n.init, d, seed = seed)$design
    response <- t(apply(design, 1, fun, ... = ...))
    
    if (!is.null(noise.var)) {
      if (typeof(noise.var) == "closure") {
        newnoise.var <- apply(design, 1, noise.var, ...)
      } else if (typeof(noise.var) == "double") {
        newnoise.var <- matrix(noise.var, nrow=n.init, ncol=ncol(response), byrow=TRUE)
      } else {
        tmp <- newnoise.var <- NULL # initialization
        for (i in 1:length(response)) {
          tmp <- rbind(tmp, response[[i]][[1]])
          newnoise.var <- rbind(newnoise.var, response[[i]][[2]])
        }
        response <- tmp
      }
    } else {
      newnoise.var <- NULL
    }
    nobj <- ncol(response)
    
    my.km <- function(i) {
      km(model.trend, design = design, response = response[, i], covtype=covtype,
         control=control, lower=kmlb, upper=kmub, nugget=kmnugget, noise.var=newnoise.var[,i])
    }
    model <- mclapply(1:nobj, my.km, mc.cores=ncores)
  }
  
  #### Integration points ###########################
  if (is.null(integcontrol$integ.pts) || is.null(integcontrol$expanded.indices)) {
    if (equilibrium %in% c("KSE", "CKSE")) include.obs <- TRUE else include.obs <- FALSE
    res <- generate_integ_pts(n.s=n.s, d=d, nobj=nobj, x.to.obj=x.to.obj, equilibrium=equilibrium,
                              gridtype=gridtype, lb = lb, ub = ub, include.obs=include.obs, model=model)
    integcontrol$integ.pts <- integ.pts <- res$integ.pts
    integcontrol$expanded.indices <- expanded.indices <- res$expanded.indices
  }
  if (is.null(expanded.indices)) sorted <- FALSE
  else sorted <- !is.unsorted(expanded.indices[, ncol(expanded.indices)])
  
  ####################################################################################################
  #### DEFINE SEQUENCE OF ACTIONS ####################################################################
  ####################################################################################################
  # J is the individual objective to consider at time ii (0 means all => equilibrium)
  J <- rep(c(1, 1:nobj, 1:nobj), ceiling((n.ite+1) / (2*nobj-1)))
  
  # Task is either 0 (equilibrium), 1 (shadow), 2 (nadir), or 3 (variance)
  if (equilibrium == "KSE") {
    task <- rep(c(0, rep(1, nobj), rep(2, nobj)), ceiling((n.ite+1) / (2*nobj-1)))
  } else {
    task <- rep(c(0, rep(3, nobj), rep(3, nobj)), ceiling((n.ite+1) / (2*nobj-1)))
  }
  
  ####################################################################################################
  #### MAIN LOOP STARTS HERE #########################################################################
  ####################################################################################################
  
  for (ii in 1:(n.ite+1)){
    
    ##################################################################
    # RENEW INTEGRATION POINTS
    ##################################################################
    if (ii > 1 && integcontrol$renew) {
      include.obs <- equilibrium %in% c("KSE", "CKSE")
      res <- generate_integ_pts(n.s=n.s, d=d, nobj=nobj, x.to.obj=x.to.obj, equilibrium=equilibrium,
                                gridtype=gridtype, lb = lb, ub = ub, include.obs=include.obs, model=model, 
                                init_set=init_set, seed=ii)
      integcontrol$integ.pts <- integ.pts <- res$integ.pts
      integcontrol$expanded.indices <- expanded.indices <- res$expanded.indices
    }
    
    ##################################################################
    # MODELS UPDATE
    ##################################################################
    if (ii > 1) {
      newmodel <- model
      X.new <- matrix(xnew,nrow=1)
      my.update <- function(u) {
        try(update(object = model[[u]], newX = X.new, newy=ynew[u], newX.alreadyExist=FALSE, newnoise.var = newnoise.var[u],
                   cov.reestim = cov.reestim, kmcontrol = list(control = list(trace = FALSE))), silent = TRUE)
      }
      newmodel <- mclapply(1:nobj, my.update, mc.cores=ncores)
      
      for (u in 1:nobj){
        if (typeof(newmodel[[u]]) == "character" && cov.reestim) {
          cat("Error in hyperparameter estimation - old hyperparameter values used instead for model ", u, "\n")
          newmodel[[u]] <- try(update(object = model[[u]], newX = X.new, newy=ynew[u], newnoise.var = newnoise.var[u],
                                      newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
        }
        if (typeof(newmodel[[u]]) == "character") {
          cat("Unable to udpate kriging model ", u, " at iteration", ii-1, "- optimization stopped \n")
          cat("lastmodel ", u, " is the model at iteration", ii-1, "\n")
          cat("par and values contain the ",ii , "th observation \n \n")
          return(list(par = X.new, values = ynew, nsteps = ii, model = model, integcontrol=integcontrol))
        } else {
          model[[u]] <- newmodel[[u]]
        }
      }
    }
    observations <- Reduce(cbind, lapply(model, slot, "y"))
    
    ##################################################################
    # PREDICTION
    ##################################################################
    pred <- mclapply(model, FUN=predict, newdata = integ.pts, checkNames = FALSE, type = "UK", light.return = TRUE, mc.cores=ncores)
    
    ##################################################################
    # ACQUISITION
    ##################################################################
    jj <- J[ii]
    get_shadow <- task[ii] == 1
    get_nadir <- task[ii] == 2
    max_predvar <- task[ii] == 3
    
    discard <- which(pred[[jj]]$sd/sqrt(model[[jj]]@covariance@sd2) < 1e-06)
    
    if (max_predvar) {
      ################################################################
      # Maximise prediction variance
      i <- which.max(pred[[jj]]$sd)
      
    } else if (get_shadow) {
      ################################################################
      # Search for Shadow
      if (trace>1) cat("Looking for individual minimum of objective: ", jj, "\n")
      if (is.null(calibcontrol$target)) {
        # EI(min) on each objective (to find potential Utopia points)
        xcr <-  (min(model[[jj]]@y) - pred[[jj]]$mean)/pred[[jj]]$sd
        test <- (min(model[[jj]]@y) - pred[[jj]]$mean)*pnorm(xcr) + pred[[jj]]$sd * dnorm(xcr)
      } else {
        # EI(min) on each objective (to find potential Utopia points)
        zmin <- min((model[[jj]]@y - target[jj])^2)
        mu <- pred[[jj]]$mean - target[jj]
        sigma   <- pred[[jj]]$sd
        a2 <- (sqrt(zmin) - mu)/sigma
        a1 <- (-sqrt(zmin) - mu)/sigma
        da2 <- dnorm(a2);   da1 <- dnorm(a1)
        pa2 <- pnorm(a2);   pa1 <- pnorm(a1)
        test <- (zmin - sigma^2 - mu^2) * (pa2 - pa1) + sigma^2 * (a2*da2 - a1*da1) + 2*mu*sigma*(da2 - da1)
      }
      test[discard] <- NA
      i <- which.max(test)
      
    } else if (get_nadir) {
      ################################################################
      # Search for Nadir
      if (trace>1) cat("Looking for nadir of objective: ", jj, "\n")
      PFobs <- nonDom(Reduce(cbind, lapply(model, slot, "y")))
      PNDom <- prob.of.non.domination(model=model, integration.points=integ.pts, predictions=pred, target=calibcontrol$target)
      
      if (is.null(calibcontrol$target)) {
        # EI(max) x Pdom on each objective (to find potential Nadir points) unless Nadir is provided
        xcr <-  -(max(PFobs[,jj]) - pred[[jj]]$mean)/pred[[jj]]$sd
        test <- (-max(PFobs[,jj]) + pred[[jj]]$mean)*pnorm(xcr) + pred[[jj]]$sd * dnorm(xcr)
      } else {
        zmax <- max(PFobs[,jj])
        b2 <- (sqrt(zmax) - mu)/sigma
        b1 <- (-sqrt(zmax) - mu)/sigma
        db2 <- dnorm(b2);   db1 <- dnorm(b1)
        pb2 <- pnorm(b2);   pb1 <- pnorm(b1)
        test <- (sigma^2 + mu^2 - zmax) * (1 - pb2 + pb1) + sigma^2 * (b2*db2 - b1*db1) + 2*mu*sigma*(db2 - db1)
      }
      test[discard] <- NA
      test <- test * PNDom
      i <- which.max(test)
      
      ################################################################
      # Search for KS
    } else {
      if (trace>1) cat("Looking for equilibrium \n")
      predmean <- Reduce(cbind, lapply(pred, function(alist) alist$mean))
      if (!is.null(calibcontrol$target)) {
        # Calibration mode
        predmean <- (predmean - matrix(rep(calibcontrol$target, nrow(predmean)), byrow=TRUE, nrow=nrow(predmean)))^2
        if (calibcontrol$log) {
          predmean <- log(predmean + calibcontrol$offset)
        }
      }
      currentEq <- try(getEquilibrium(Z = predmean,  equilibrium = equilibrium, nobj = nobj,
                                      expanded.indices=expanded.indices, n.s=n.s, kweights = NULL,
                                      sorted = sorted, cross = cross, return.design = TRUE, Nadir=Nadir, Shadow=Shadow))
      i <- currentEq$NE
    }
    
    ##################################################################
    # GET NEW OBSERVATION
    ##################################################################
    xnew <- integ.pts[i,]
    while (duplicated(rbind(xnew, model[[1]]@X), fromLast = TRUE)[1]) xnew <- integ.pts[sample.int(nrow(integ.pts), 1),]
    
    ## Compute new observation
    ynew <- try(fun(xnew, ...))
    
    if (typeof(ynew) == "character" ) {
      cat("Unable to compute objective function at iteration ", i, "- optimization stopped \n")
      cat("Problem occured for the design: ", xnew, "\n")
      cat("Last model and problematic design (xnew) returned \n")
      return(list(model=model, predEq=predEq, integcontrol=integcontrol, xnew=xnew))
    }
    
    if (!is.null(noise.var)) {
      if (typeof(noise.var) == "closure") {
        newnoise.var <- noise.var(xnew)
      } else if (typeof(noise.var) == "double") {
        newnoise.var <- noise.var
      } else {#noise.var ="given_by_fn"
        newnoise.var <- ynew[[2]]
        ynew <- ynew[[1]]
      }
    } else {
      newnoise.var <- NULL
    }
  }
  ####################################################################################################
  #### MAIN LOOP ENDS HERE ###########################################################################
  ####################################################################################################
  return(list(model = model, integcontrol=integcontrol))
}