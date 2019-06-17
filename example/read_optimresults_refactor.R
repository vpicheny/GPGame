library(GPareto)
library(GPGame)
library(plot3D)
library(rgl)
library(GGally)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)
library(DiceDesign)
setwd("~/Code/GPGame")
source("../Gama/multiplot.R")
source('~/Code/GPGame/example/plot_utils.R')

set.seed(42)
n.init <- 100
ncores = 9

# Set directory
load("../Gama/v2.0/optim_v2.0_KSE.RData")
mod_ks <- res$model
Eq.design_ks <- res$Eq.design
Eq.poff_ks <- res$Eq.poff

load("../Gama/v2.0/optim_v2.0_CKSE.RData")
mod_cks <- res$model
Eq.design_cks <- res$Eq.design
Eq.poff_cks <- res$Eq.poff

load("../Gama/v2.0/optim_v2.0_PKSE.RData")
mod_pks <- res$model
Eq.design_pks <- res$Eq.design
Eq.poff_pks <- res$Eq.poff

res <- aggregate_obs(list(mod_ks, mod_cks, mod_pks), n.init=n.init)
observations <- res$observations
X <- res$X
noise.var <- res$noise.var
model <- build_models(X, observations, noise.var, ncores)

calibcontrol <- list(target=rep(1, nobj), log=TRUE, offset=0.01)

NS <- get_nadir_and_shadow(observations, model, n.s=1e5, nsimPoints=500, calibcontrol=calibcontrol, nsim=10, ncores=ncores)

ref_set <- get_reference_set(model, n.s=1e5, nsimPoints=500, calibcontrol, nsim=10, ncores=ncores)

obs_ks <- Reduce(cbind, lapply(mod_ks, slot, "y"))
obs_ks <- transform_obs(obs_ks, calibcontrol)

obs_cks <- Reduce(cbind, lapply(mod_cks, slot, "y"))
obs_cks <- transform_obs(obs_cks, calibcontrol)

obs_pks <- Reduce(cbind, lapply(mod_pks, slot, "y"))
obs_pks <- transform_obs(obs_pks, calibcontrol)

plot_centrality(rbind(obs_ks, Eq.poff_ks), Nadir=NS$Nadir, Shadow=NS$Shadow, n.init=100, add=FALSE, copula=FALSE)

plot_centrality(rbind(obs_cks, Eq.poff_cks), observations_ref=ref_set$observations, n.init=100, add=FALSE, copula=TRUE)

A = convert_to_ranks(obs_cks, ref_set$observations[[1]])
apply(ref_set$observations[[1]], 2, range)

plot_pairs(transform_obs(observations, calibcontrol), Eq.poff_ks, n.init)
