library(GPareto)
library(GPGame)
library(ggplot2)
library(parallel)
library(gridExtra)
require(ggplot2)
require(reshape2)

setwd("~/Code/GPGame")
source('~/Code/GPGame/example/plot_utils.R')

testfun <- "DTLZ2" # "hartman" "DTLZ2"
# config <- "L" # "S", "M", "L", "XL"
# pb_type <- "continuous" # "discrete", "continuous"
equilibrium <- "KS" #"KS", "CKS"



config_number1 <-  paste0("S", "_", "discrete")
config_number2 <-  paste0("L", "_", "discrete")
config_number3 <-  paste0("baseline", "_", "discrete")
config_number4 <-  paste0("S", "_", "continuous")
config_number5 <-  paste0("L", "_", "continuous")

config_number <- c(config_number1, config_number2, config_number3, config_number4, config_number5)

if (testfun == "DTLZ2") { directory <- "~/Code/GPGame/example/Test_results/dtlz2/" 
} else {                  directory <- "~/Code/GPGame/example/Test_results/hartman/"}

average_perf <- all_perf <- min_perf <- max_perf <- c()
for (i in 1:length(config_number)) {
  exp_name <- paste0(directory, "config_", config_number[i], "_")
  
  if (i != 3) load(paste0(exp_name, "config_and_solution.RData"))
  
  ntests = 3
  list_of_model <- vector("list", ntests) 
  
  for (ii in 1:ntests) {
    load(paste0(exp_name, "KSE_run_", ii, ".RData"))
    list_of_model[[ii]] <- res$model
  }
  solution <- KS_act
  print(KS_act)
  
  perf <- convergence_plots(list_of_model, solution)
  average_perf <- rbind(average_perf, apply(perf$cummin, 2, mean))
  min_perf <- rbind(min_perf, apply(perf$cummin, 2, min))
  max_perf <- rbind(max_perf, apply(perf$cummin, 2, max))
  all_perf <- rbind(all_perf, perf$cummin)
}

t_max <- ncol(average_perf)
# t_max <- 42

df <- data.frame(t(average_perf[,1:t_max, drop=FALSE]))
df$time = 1:t_max
names(df) <- c(config_number, "time")
df <- melt(df ,  id.vars = 'time', variable.name = 'config')

dfm <- data.frame(t(min_perf[,1:t_max, drop=FALSE]))
dfm$time = 1:t_max
names(dfm) <- c(config_number, "time")
dfm <- melt(dfm,  id.vars = 'time', variable.name = 'lo')

dfM <- data.frame(t(max_perf[,1:t_max, drop=FALSE]))
dfM$time = 1:t_max
names(dfM) <- c(config_number, "time")
dfM <- melt(dfM,  id.vars = 'time', variable.name = 'up')

ddf = cbind(df, dfm, dfM)[-c(4,7)]

p <- ggplot(df, aes(time, value, fill=config)) + geom_line(aes(colour = config)) + ylim(0, max(average_perf))
p <- p + geom_vline(xintercept=n.init+1/2, linetype="dashed", color = "black")
grid.arrange(p, nrow = 1)

pp <- ggplot(ddf, aes(time, value, ymin=value.1, ymax=value.2, fill=config)) + geom_line(aes(colour = config)) + ylim(0, max(average_perf))
pp <- pp + geom_vline(xintercept=n.init+1/2, linetype="dashed", color = "black")
pp <- pp + geom_ribbon(alpha=0.5)
grid.arrange(pp, nrow = 1)
