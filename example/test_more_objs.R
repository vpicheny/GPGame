## Test of the interest of considering all objectives vs aggregations

library(GPGame)
library(DiceKriging)
set.seed(42)

# 6 objs function
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

# 2 objs function
fun2 <- function(x, group1 = c(2,5,6), group2 = c(1,3,4)){
  tmp <- fun(x)
  return(c(sum(tmp[group1]), sum(tmp[c(group2)])))
}

set.seed(42)

nset <- 5e5
Xset <- matrix(runif(nset * 6), nset) 
Yset <- t(apply(Xset, 1, fun))
group1 <- c(2,5,6)
group2 <- c(1,3,4)

library(mco) # to compute 2d (aggregated) and 6d reference PFs
solref <- nsga2(fun, 6, 6, lower.bounds = rep(0, 6), upper.bounds = rep(1,6), popsize = 6000)

solref2 <- nsga2(fun2, 6, 2, lower.bounds = rep(0, 6), upper.bounds = rep(1,6), popsize = 2000)
Yall <- rbind(Yset, solref$value, t(apply(solref2$par, 1, fun)))


# Compute (pure) 2d PFs
for(i in 1:5){
  for(j in (i+1):6){
    sol2d <- nsga2(fun2, 6, 2, lower.bounds = rep(0, 6), upper.bounds = rep(1,6), popsize = 200,
                   group1 = i, group2 = j)
    Yall <- rbind(Yall, t(apply(sol2d$par, 1, fun)))
  }
}

# 3 objs function
fun3 <- function(x, group1 = c(2,5), group2 = c(1,3), group3 = c(4,6)){
  tmp <- fun(x)
  return(c(sum(tmp[group1]), sum(tmp[c(group2)]), sum(tmp[c(group3)])))
}

# Compute (pure) 3d PFs
for(i in 1:4){
  for(j in (i+1):5){
    for(k in (j+1):6){
      sol3d <- nsga2(fun3, 6, 3, lower.bounds = rep(0, 6), upper.bounds = rep(1,6), popsize = 300,
                     group1 = i, group2 = j, group3 = k)
      Yall <- rbind(Yall, t(apply(sol3d$par, 1, fun)))
    }
  }
}

# 4 objs function
fun4 <- function(x, group1 = c(2,5), group2 = c(1,3), group3 = c(4), group4 = c(6)){
  tmp <- fun(x)
  return(c(sum(tmp[group1]), sum(tmp[c(group2)]), sum(tmp[c(group3)]), sum(tmp[c(group4)])))
}

# Compute (pure) 4d PFs
for(i in 1:3){
  for(j in (i+1):4){
    for(k in (j+1):5){
      for(l in (k+1):6){
        sol4d <- nsga2(fun3, 6, 4, lower.bounds = rep(0, 6), upper.bounds = rep(1,6), popsize = 400,
                       group1 = i, group2 = j, group3 = k, group4 = l)
        Yall <- rbind(Yall, t(apply(sol4d$par, 1, fun)))
      }
    }
  }
}

# 5 objs function
fun5 <- function(x, group1 = c(2,5), group2 = c(1), group3 = c(4), group4 = c(6), group5 = c(3)){
  tmp <- fun(x)
  return(c(sum(tmp[group1]), sum(tmp[c(group2)]), sum(tmp[c(group3)]), sum(tmp[c(group4)]), sum(tmp[c(group5)])))
}

# Compute (pure) 5d PFs
for(i in 1:2){
  for(j in (i+1):3){
    for(k in (j+1):4){
      for(l in (k+1):5){
        for(m in (l+1):6){
          sol5d <- nsga2(fun3, 6, 5, lower.bounds = rep(0, 6), upper.bounds = rep(1,6), popsize = 500,
                         group1 = i, group2 = j, group3 = k, group4 = l)
          Yall <- rbind(Yall, t(apply(sol5d$par, 1, fun)))
        }
      }
    }
  }
}


Pset <- nonDom(Yall, return.idx = TRUE)

KSeq <- getEquilibrium(Z = Yall, equilibrium = "KSE", nobj = 6)
Shadow <- apply(Yall, 2, min)
Nadir <- apply(Yall[Pset,], 2, max)
iKS <- which(Yall[,1] == KSeq[1])

cols <- rep(1, length(Pset))
pairs(rbind(Yall[Pset,], Yall[iKS,]), col = c(cols,2), pch = c(rep(1, length(Pset)),20))

## Now aggregate/select objectives
cor(Yset) # highest correlations: 5-6, 2-4, 2-5, 3-5, 4-5
Zall <- cbind(rowSums(Yall[,group1]), rowSums(Yall[,group2]))

Pset2 <- nonDom(Zall, return.idx = TRUE)
cols2 <- rep(1, nrow(Zall))
cols2[Pset2] <- 3
plot(Zall, col = cols2)
points(Zall[iKS,1], Zall[iKS,2], col = 'red', pch = 20)

## Calcul de la distance sur couple d'objectifs

par(mfrow = c(5, 6))
diff_all <- diff_12 <- c()

for (id1 in 1:5){
  for(id2 in (id1+1):6){
    
    Pset_12 <- nonDom(Yall[,c(id1,id2)], return.idx = TRUE)
    ratios <- (matrix(Nadir, nrow = length(Pset_12), ncol = 6, byrow = T) - Yall[Pset_12,]) %*% diag(1/(Nadir - Shadow))
    
    ratio_all <- apply(ratios, 1, min)
    ratio_KS <- min((Nadir - KSeq)/(Nadir - Shadow))
    diff_all <- c(diff_all, max(ratio_all) - ratio_KS)
    
    plot(ratio_all, ylim = c(0,1), main = paste('6 obj results, PF obj_1:', id1, 'obj_2:', id2, ' diff:', signif(max(ratio_all) - ratio_KS, 3)))
    abline(h = ratio_KS, col = 'red')
    
    
    ratio_12 <- (matrix(Nadir, nrow = length(Pset_12), ncol = 6, byrow = T) - Yall[Pset_12,]) %*% diag(1/(Nadir - Shadow))
    ratio_12 <- apply(ratio_12[,c(id1, id2)], 1, min)
    ratio_KS_12 <- min(((Nadir - KSeq)/(Nadir - Shadow))[c(id1, id2)])
    
    plot(ratio_12, ylim = c(0,1), main = paste('2 obj results, PF obj_1:', id1, 'obj_2:', id2, ' diff:', signif(max(ratio_12) - ratio_KS_12, 3)))
    abline(h = ratio_KS_12, col = 'blue')
    
    diff_12 <- c(diff_12, max(ratio_12) - ratio_KS_12)
  }
}
par(mfrow = c(1,1))

## Putting things on the same graphs
if(save) pdf("InterestMaO.pdf", width = 12, height = 12)
par(mfrow = c(4, 4))
for (id1 in 1:5){
  for(id2 in (id1+1):6){
    
    Pset_12 <- nonDom(Yall[,c(id1,id2)], return.idx = TRUE)
    ratios <- (matrix(Nadir, nrow = length(Pset_12), ncol = 6, byrow = T) - Yall[Pset_12,]) %*% diag(1/(Nadir - Shadow))
    
    ratio_all <- apply(ratios, 1, min)
    ratio_KS <- min((Nadir - KSeq)/(Nadir - Shadow))
    
    ratio_12 <- (matrix(Nadir, nrow = length(Pset_12), ncol = 6, byrow = T) - Yall[Pset_12,]) %*% diag(1/(Nadir - Shadow))
    ratio_12 <- apply(ratio_12[,c(id1, id2)], 1, min)
    ratio_KS_12 <- min(((Nadir - KSeq)/(Nadir - Shadow))[c(id1, id2)])
    
    plot(ratio_all, ylim = c(0,1), pch = 3, col = "black", ylab = "ratio", xlab = "Index of solution",
         main = substitute(paste("PF (", f[a], ",", f[b], ")"
                                 #, ', diff KS/PF 6objs:', c,
                                 # ' diff KS/PF 2objs:', d
         ), 
         list(a = id1, b = id2
              #, c = signif(max(ratio_all) - ratio_KS, 3), d = signif(max(ratio_12) - ratio_KS_12, 3)
         )))
    abline(h = ratio_KS, col = 'red', lty = 2)
    
    
    points(ratio_12, col = "green", pch = 1) #main = paste('2 obj results, PF obj_1:', id1, 'obj_2:', id2, ' diff:', signif(max(ratio_12) - ratio_KS_12, 3)))
    abline(h = ratio_KS_12, col = 'blue', lty = 3)
    
    
  }
}
plot(NA, NA, ylab = "", xlab = "")
legend("bottom", col = c(1, 3, 2, 4), pch = c(3, 1, NA, NA), lty = c(NA, NA, 2, 3), 
       legend = c("6 objectives", "2 objectives", "KS 6 objectives", "KS 2 objectives"))
par(mfrow = c(1,1))
if(save) dev.off()

if(save) pdf("BoxplotMaO.pdf", width = 6, height = 6)
boxplot(-diff_all, diff_12, names = c("-Diff on 6 objs", "Diff on 2 objs"))
if(save) dev.off()
