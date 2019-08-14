library(GPareto)
library(GPGame)
library(ggplot2)
library(parallel)
library(gridExtra)
require(ggplot2)
require(reshape2)

# setwd("~/Code/GPGame")
source('example/plot_utils.R')

testfun <- "DTLZ2" # "hartman" "DTLZ2"
# config <- "L" # "S", "M", "L", "XL"
pb_type <- "discrete" # "discrete", "continuous"
equilibrium <- "KSE" #"KSE", "CKSE"
ratios <- TRUE
copula <- FALSE

config_number1 <-  paste0("S", "_", "discrete")
# config_number2 <-  paste0("M", "_", "discrete")
config_number3 <-  paste0("L", "_", "discrete")
config_number4 <-  paste0("baseline", "_", "discrete")
config_number5 <-  paste0("RS", "_", "discrete")
config_number6 <-  paste0("SMS", "_", "discrete")

config_number7 <-  paste0("S", "_", "continuous")
# config_number8 <-  paste0("M", "_", "continuous")
config_number9 <-  paste0("L", "_", "continuous")
config_number10 <-  paste0("baseline", "_", "continuous")

# config_number4 <-  paste0("RS", "_", "discrete")
# config_number5 <-  paste0("baseline", "_", "discrete")
# config_number6 <-  paste0("SMS", "_", "discrete")

# config_number4 <- NULL
# config_number3 <-  paste0("baseline", "_", "discrete")
# config_number4 <-  paste0("RS", "_", "discrete")
# config_number5 <-  paste0("S", "_", "continuous")
# config_number6 <-  paste0("L", "_", "continuous")
# config_number7 <-  paste0("baseline", "_", "continuous")
# config_number8 <-  paste0("RS", "_", "discrete")

# config_number <- c(config_number1, config_number2, config_number4, config_number5)
# config_number <- c(config_number1, config_number2, config_number3, config_number4, 
#                    config_number5, config_number6, config_number7)
if (pb_type == "continuous") {
  config_number <- c(config_number7, config_number9, config_number10)
} else {
  config_number <- c(config_number1, config_number3, config_number4, config_number5, config_number6)
}

if (testfun == "DTLZ2") { directory <- "example/Test_results/dtlz2/" 
} else {                  directory <- "example/Test_results/hartman/"}

average_perf <- all_perf <- median_perf <- min_perf <- max_perf <- c()

load(paste0(directory, pb_type, "_config_and_solution.RData"))

for (i in 1:length(config_number)) {
  exp_name <- paste0(directory, "config_", config_number[i], "_")
  
  # if (i %in% c(1, 5)) load(paste0(exp_name, "config_and_solution.RData"))
  
  ntests = 10
  list_of_model <- vector("list", ntests) 
  
  for (ii in 1:ntests) {
    load(paste0(exp_name, equilibrium, "_run_", ii, ".RData"))
    list_of_model[[ii]] <- res$model
  }
  if(equilibrium == "KSE") solution <- KS_act else solution <- CKS_act
  print(solution)
  
  # if (i == 8) exp_name <- paste0("RS", "_", "continuous")
  
  if (ratios) {
    perf <- convergence_plots(models = list_of_model, solution = solution, Nadir = Nadir, Shadow = Shadow, title=config_number[i],
                              observations_ref = fun.grid, copula = copula)
  } else {
    perf <- convergence_plots(models = list_of_model, solution = solution, title=config_number[i], copula = copula)
  }
 
  average_perf <- rbind(average_perf, apply(perf$cummin, 2, mean))
  median_perf <- rbind(median_perf, apply(perf$cummin, 2, median))
  min_perf <- rbind(min_perf, apply(perf$cummin, 2, min))
  max_perf <- rbind(max_perf, apply(perf$cummin, 2, max))
  all_perf <- rbind(all_perf, perf$cummin)
}

ymin <- 5e-3
myperf <- average_perf
myperf <- pmax(myperf, ymin)
min_perf <- pmax(min_perf, ymin)
t_max <- ncol(myperf)
# t_max <- 42

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
col= gg_color_hue(6)

# config_number[8] <- paste0("RS", "_", "continuous")
# config_number <- c("S", "L", "base", "RS")
if (pb_type == "continuous") {
  config_number <- c("SUR coarse", "SUR fine", "Baseline")
} else {
  config_number <- c("SUR coarse", "SUR fine",  "Baseline", "Random", "Hypervolume")
}
cols = col[1:length(config_number)]

df <- data.frame(t(myperf[,1:t_max, drop=FALSE]))
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

names(df)[names(df) == "value"] <- "Gap"

p <- ggplot(df, aes(time, Gap, fill=config)) + geom_line(aes(colour = config)) + ylim(0, NA)
p <- p + geom_vline(xintercept=n.init+1/2, linetype="dashed", color = "black")
p <- p + ggtitle(paste0(testfun, " ", pb_type))
grid.arrange(p, nrow = 1)

if (equilibrium == "KSE") equilibrium <- "KS" else equilibrium <- "CKS" 

names(ddf)[names(ddf) == "value"] <- "Gap"
# names(ddf)[names(ddf) == "config"] <- "Algorithm"

# dddf = ddf[-which(ddf$config %in% c("RS_discrete", "RS_continuous")),]
pp <- ggplot(ddf, aes(time, Gap, ymin=value.1, ymax=value.2, fill=config))
pp <- pp + geom_vline(xintercept=n.init+1/2, linetype="dashed", color = "black")
pp <- pp + geom_ribbon(alpha=0.2)
pp <- pp + geom_line(aes(colour = config), size=1) #+ ylim(0.001, NA)
pp <- pp + scale_y_continuous(trans='log10', limits=c(ymin, .5)) #+ ylim(0.001, NA)
pp <- pp + ggtitle(paste0(testfun, " - ", pb_type, " - ", equilibrium)) + labs(x = "Nb of observations") 
pp <- pp + scale_color_manual(values = cols) + scale_fill_manual(values = cols) + theme(legend.title=element_blank())
# pp <- pp + guides(colour=guide_legend(title="Algorithm"))
# pp <- pp + coord_trans(y='log10')
pp <- pp #+ ylim(0.001, NA)
grid.arrange(pp, nrow = 1)
