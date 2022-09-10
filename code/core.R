### PROGETTO PER IL CORSO DI DATI FUNZIONALI
# Francesco Martella, Sara Meneghetti, Alessio Piraccini

rm(list = ls())
set.seed(42)
opar <- par()


# libraries --------------------------------------------------------------------

library(fda)
library(refund)
library(tidyverse)


# utils ------------------------------------------------------------------------

proj_dir  <-  "~/R_Jobs/DatiFunzionali/progetto_funzionali"
if (normalizePath(getwd()) != normalizePath(proj_dir)) setwd(proj_dir)

load("data/finaldata.RData")
ls()

### function fo extacting MISE
get_err <- function(Y, X){
  iae <- ise <- numeric(NROW(Y))
  for(i in 1:NROW(Y)){
    cat(i, "... ")
    test <- i
    y <- Y[-test,]
    x <- X[-test,]
    d <- data.frame(I(y), I(x))
    m <- pffr(y ~ x,
              yind = t,
              bs.yindex = list(bs = "ps", k = 40, m = c(2, 1)),
              method="REML", data = d)
    d2 <- data.frame(x = I(t(X[test,])))
    pi <- predict(m, newdata = d2)
    iae[i] <- sum(abs(Y[test,] - drop(pi)))
    ise[i] <- sum((Y[test,] - drop(pi))^2)
  }
  out <- c(mean(iae), mean(ise))
  names(out) <- c("MIAE", "MISE")
  out
}

# r2
r2 <- function(yhat, yoss = Y){
  num <- apply((yoss - yhat)^2, 2, sum)
  den <- apply((yoss - apply(yoss, 2, mean))^2, 2, sum)
  r2 <- 1 - num/den
}


# deaths on lagged cases -------------------------------------------------------

# indici temporali
grid1 
s2 <- 1:115
s <- 2:116
t <- 3:117
# 117 settimane, per fare un lag ne tolgo 1!

# dati
Y <- t(as.matrix(deaths$mat.sm))
Y <- Y[,t]

X1 <- t(as.matrix(cases$mat.sm))
X1 <- X1[,s]

set.seed(42)
dat<-data.frame(I(Y), I(X1))

# mod
m1 <- pffr(Y ~ X1,
           yind = t,
           bs.yindex = list(bs = "ps", k = 40, m = c(2, 1)),
           method="REML", data = dat)
summary(m1) #edf effective degrees of freedom of the fit
plot(m1, scheme=T, pages=1, scale = 0)
#plot(m1, scheme=T, pages=1, scale = -1) # stessa scala

# x report
pdf("plots/beta0m1.pdf", width = 7, height = 5)
plot(m1, scheme=1, select = 1, rug = F, xlab = "", ylab = bquote(beta[0]), scale = 0)
dev.off()

pdf("plots/beta1m1.pdf", width = 7, height = 5)
plot(m1, scheme=1, select = 2, rug = F, xlab = "", ylab = bquote(beta[1]), scale = 0)
dev.off()

err1 <- get_err(Y, X1)
err1

# extra: prova modello non simultaneo
X1 <- t(as.matrix(cases$mat.sm))
X1 <- X1[,t]

m11 <- pffr(Y ~ ff(X1, xind = t, limits = "s<t"),
           yind = t,
           bs.yindex = list(bs = "ps", k = 40, m = c(2, 1)),
           method="REML", data = dat) 
summary(m11)
plot(m11, scheme = 1, select = 1)
plot(m11, scheme = 2, select = 2)
#m11$pffr$ff$`ff(X1, xind = t, limits = "s<t")`

# x report
pdf("plots/beta0m11.pdf", width = 7, height = 5)
plot(m11, scheme=1, select = 1, rug = F, xlab = "", ylab = bquote(beta[0]), scale = 0)
dev.off()

pdf("plots/beta1m11.pdf", width = 7, height = 5)
plot(m11, scheme=2, select = 2, rug = F, xlab = "", ylab = bquote(beta[1]), scale = 0, main = "")
dev.off()

# x report
pdf("plots/res_m1.pdf", width = 7, height = 5)
matplot(grid1[t], t(residuals(m1)), type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "residui")
abline(h = 0, lty = 2)
matlines(grid1[t], apply(t(residuals(m1)), 1, mean), lty = 1, lwd = 2, col = 2)
dev.off()

pdf("plots/r2_m1.pdf", width = 7, height = 5)
plot(grid1[t], r2(fitted(m1), yoss = Y), type = "l", xlab = "", ylab = expression(R^2), col = "blue", ylim = c(0, 1))
abline(h = 0, lty = 2)
dev.off()

pdf("plots/res_m11.pdf", width = 7, height = 5)
matplot(grid1[t], t(residuals(m11)), type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "residui")
abline(h = 0, lty = 2)
matlines(grid1[t], apply(t(residuals(m11)), 1, mean), lty = 1, lwd = 2, col = 2)
dev.off()

pdf("plots/r2_m11.pdf", width = 7, height = 5)
plot(grid1[t], r2(fitted(m11), yoss = Y), type = "l", xlab = "", ylab = expression(R^2), col = "blue", ylim = c(0, 1))
abline(h = 0, lty = 2)
dev.off()

dim(fitted(m1))


# cases on other covariates ----------------------------------------------------

# modello  simultaneo laggato
# dati
Y <- t(as.matrix(cases$mat.sm))
Y <- Y[,s]

X1 <- t(as.matrix(residential$mat.sm))
X1 <- X1[,s2]
X2 <- t(as.matrix(grocery_and_pharmacy$mat.sm))
X2 <- X2[,s2]
X3 <- t(as.matrix(transit$mat.sm))
X3 <- X3[,s2]
X4 <- t(as.matrix(retail_and_recreation$mat.sm))
X4 <- X4[,s2]
X5 <- t(as.matrix(parks$mat.sm))
X5 <- X5[,s2]
X6 <- t(as.matrix(workplaces$mat.sm))
X6 <- X6[,s2]

set.seed(42)
dat <- data.frame(I(Y), I(X1), I(X2), I(X3), I(X4), I(X5), I(X6))

# mod
m2  <-  pffr(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
             yind = t,
             bs.yindex = list(bs = "ps", k = 40, m = c(2, 1)),
             method="REML", data = dat)
summary(m2)
plot(m2, scheme=T, pages=1, scale = 0)

dim(fitted(m2))

# residui e r2
matplot(grid1[s], t(residuals(m2)), type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "residui")
abline(h = 0, lty = 2)
matlines(grid1[s], apply(t(residuals(m2)), 1, mean), lty = 1, lwd = 2, col = 2)

r2 <- function(yhat, yoss = Y){
  num <- apply((yoss - yhat)^2, 2, sum)
  den <- apply((yoss - apply(yoss, 2, mean))^2, 2, sum)
  r2 <- 1 - num/den
}
plot(grid1[s], r2(fitted(m2)), type = "l", xlab = "", ylab = expression(R^2), col = "blue", ylim = c(0, 1))
abline(h = 0, lty = 2)

# x report
pdf("plots/res_m2.pdf", width = 7, height = 5)
matplot(grid1[s], t(residuals(m2)), type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "residui")
abline(h = 0, lty = 2)
matlines(grid1[s], apply(t(residuals(m2)), 1, mean), lty = 1, lwd = 2, col = 2)
dev.off()

pdf("plots/r2_m2.pdf", width = 7, height = 5)
plot(grid1[s], r2(fitted(m2)), type = "l", xlab = "", ylab = expression(R^2), col = "blue", ylim = c(0, 1))
abline(h = 0, lty = 2)
dev.off()

# modello non simultaneo
# dati
Y <- t(as.matrix(cases$mat.sm))
Y <- Y[,s]

X1 <- t(as.matrix(residential$mat.sm))
X1 <- X1[,s]
X2 <- t(as.matrix(grocery_and_pharmacy$mat.sm))
X2 <- X2[,s]
X3 <- t(as.matrix(transit$mat.sm))
X3 <- X3[,s]
X4 <- t(as.matrix(retail_and_recreation$mat.sm))
X4 <- X4[,s]
X5 <- t(as.matrix(parks$mat.sm))
X5 <- X5[,s]
X6 <- t(as.matrix(workplaces$mat.sm))
X6 <- X6[,s]

set.seed(42)
dat <- data.frame(I(Y), I(X1), I(X2), I(X3), I(X4), I(X5), I(X6))

# mod
m22 <- pffr(Y ~ ff(X1, xind = t, limits = "s<t") + ff(X2, xind = t, limits = "s<t") + 
              ff(X3, xind = t, limits = "s<t") + ff(X4, xind = t, limits = "s<t") + 
              ff(X5, xind = t, limits = "s<t") + ff(X6, xind = t, limits = "s<t"),
           yind = t,
           bs.yindex = list(bs = "ps", k = 40, m = c(2, 1)),
           method="REML", data = dat)
summary(m22)
plot(m22, scheme=2, pages=1, scale = 0)

dim(fitted(m22))

pdf("plots/res_m22.pdf", width = 7, height = 5)
matplot(grid1[s], t(residuals(m22)), type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "residui")
abline(h = 0, lty = 2)
matlines(grid1[s], apply(t(residuals(m22)), 1, mean), lty = 1, lwd = 2, col = 2)
dev.off()

pdf("plots/r2_m22.pdf", width = 7, height = 5)
plot(grid1[s], r2(fitted(m22), yoss = Y), type = "l", xlab = "", ylab = expression(R^2), col = "blue", ylim = c(0, 1))
abline(h = 0, lty = 2)
dev.off()


# deaths on lagged predicted cases ---------------------------------------------

# dati
Y <- t(as.matrix(deaths$mat.sm))
Y <- Y[,t]

X1 <- fitted(m2)
X2 <- fitted(m22)

set.seed(42)
dat <- data.frame(I(Y), I(X1), I(X2))

# mod
m3  <-  pffr(Y ~ X1,
             yind = t,
             bs.yindex = list(bs = "ps", k = 40, m = c(2, 1)),
             method="REML", data = dat)
summary(m3)
plot(m3, scheme=T, pages=1, scale = 0)

# x report
pdf("plots/beta0m3.pdf", width = 7, height = 5)
plot(m3, scheme=1, select = 1, rug = F, xlab = "", ylab = bquote(beta[0]), scale = 0)
dev.off()

pdf("plots/beta1m3.pdf", width = 7, height = 5)
plot(m3, scheme=1, select = 2, rug = F, xlab = "", ylab = bquote(beta[1]), scale = 0, main = "")
dev.off()

# mod
m33  <-  pffr(Y ~ X2,
             yind = t,
             bs.yindex = list(bs = "ps", k = 40, m = c(2, 1)),
             method="REML", data = dat)
summary(m33)
plot(m33, scheme=T, pages=1, scale = 0)

# x report
pdf("plots/beta0m33.pdf", width = 7, height = 5)
plot(m33, scheme=1, select = 1, rug = F, xlab = "", ylab = bquote(beta[0]), scale = 0)
dev.off()

pdf("plots/beta1m33.pdf", width = 7, height = 5)
plot(m33, scheme=1, select = 2, rug = F, xlab = "", ylab = bquote(beta[1]), scale = 0, main = "")
dev.off()

# x report
pdf("plots/res_m3.pdf", width = 7, height = 5)
matplot(grid1[t], t(residuals(m3)), type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "residui")
abline(h = 0, lty = 2)
matlines(grid1[t], apply(t(residuals(m3)), 1, mean), lty = 1, lwd = 2, col = 2)
dev.off()

pdf("plots/r2_m3.pdf", width = 7, height = 5)
plot(grid1[t], r2(fitted(m3), yoss = Y), type = "l", xlab = "", ylab = expression(R^2), col = "blue", ylim = c(0, 1))
abline(h = 0, lty = 2)
dev.off()

pdf("plots/res_m33.pdf", width = 7, height = 5)
matplot(grid1[t], t(residuals(m33)), type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "residui")
abline(h = 0, lty = 2)
matlines(grid1[t], apply(t(residuals(m33)), 1, mean), lty = 1, lwd = 2, col = 2)
dev.off()

pdf("plots/r2_m33.pdf", width = 7, height = 5)
plot(grid1[t], r2(fitted(m33), yoss = Y), type = "l", xlab = "", ylab = expression(R^2), col = "blue", ylim = c(0, 1))
abline(h = 0, lty = 2)
dev.off()

# errori
err2 <- get_err(Y, X1)
err22 <- get_err(Y, X2)
err1; err2; err22

# - fare cv leave one out per ottenere errore
# - altern. si potrebbe usare fRegress ma non sembra banale fare il lag come 
#   vogliamo noi, da quello che ho capito dovremmo fare un magheggio facendogli
#   rilisciare le cose partendo da mat.sm con il lag giusto