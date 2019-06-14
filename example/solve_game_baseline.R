solve_game_baseline <- function(
  fun, ..., equilibrium="NE", crit="sur", model=NULL, n.init=NULL, n.ite, d, nobj, x.to.obj=NULL, noise.var = NULL,
  Nadir=NULL, Shadow=NULL, integcontrol=NULL, simucontrol=NULL, filtercontrol=NULL, kmcontrol=NULL, returncontrol=NULL,
  ncores=1, trace=1, seed=NULL, calibcontrol=NULL, freq.exploit=1e3) {
  
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
  
  if (trace>0) cat("--------------------------\n Starting", equilibrium, "search with:", "\n",
                   crit, "strategy,","\n",
                   simu.filter, "simulation point filter,", "\n",
                   cand.filter, "candidate point filter,", "\n",
                   "among (", n.s, ") strategies \n --------------------------\n " )
  
  ####################################################################################################
  #### INITIALIZE VARIABLES AND MODELS ###############################################################
  ####################################################################################################

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
  #### MAIN LOOP STARTS HERE #########################################################################
  ####################################################################################################
  Jnres <- rep(NA, n.ite + 1)
  
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
    
    #--------------------------------------------------------#
    # Regular loop : find infill point
    #--------------------------------------------------------#
    pred <- mclapply(model, FUN=predict, newdata = integ.pts, checkNames = FALSE, type = "UK", light.return = TRUE, mc.cores=ncores)
    
    if (((ii == n.ite) && returncontrol$force.exploit.last) || ((ii %% freq.exploit == 0) && (ii > 1))) {
      if (equilibrium %in% c("KSE", "CKSE")) {
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
        
        if (typeof(currentEq)== "character") {
          cat("Unable to compute exploitation step, last model returned \n")
          return(list(model=model, predEq=predEq, integcontrol=integcontrol))
        }
        
        i <- currentEq$NE
        xnew <- integ.pts[i,]
        
        if (duplicated(rbind(integ.pts[currentEq$NE,], model[[1]]@X), fromLast = TRUE)[1]) FailToExploit <- TRUE
      }
      
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
  }
  ####################################################################################################
  #### MAIN LOOP ENDS HERE ###########################################################################
  ####################################################################################################
  return(list(model = model, integcontrol=integcontrol))
}