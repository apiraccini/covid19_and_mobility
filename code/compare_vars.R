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

# indici temporali
grid1 
s2 <- 1:115
s <- 2:116
t <- 3:117
# 117 settimane, per fare un lag ne tolgo 1!


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

plot(grid1[s], r2(fitted(m2)), type = "l", xlab = "", ylab = expression(R^2), col = "blue", ylim = c(0, 1))
abline(h = 0, lty = 2)

matplot(grid1[s], t(residuals(m2)), type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "residui")
abline(h = 0, lty = 2)
matlines(grid1[s], apply(t(residuals(m2)), 1, mean), lty = 1, lwd = 2, col = 2)

plot(grid1[s], r2(fitted(m2)), type = "l", xlab = "", ylab = expression(R^2), col = "blue", ylim = c(0, 1))
abline(h = 0, lty = 2)


plot(m2, scheme=T, pages=1, scale = 0)


######
metrics <- matrix(NA, 6, 2)
rownames(metrics) <- c("residential", "grocery_and_pharmacy", "transit",
                       "retail_and_recreation", "parks", "workplaces")
colnames(metrics) <- c("MIAE", "MISE")
for(i in 1:6){
  x <- get(paste0("X", i, collapse = NULL))
  metrics[i,] <- get_err(Y, x)
}
knitr::kable(metrics, digits = 4)


######
pdf("plots/m2_b2.pdf", width = 7, height = 5)
plot(m2, scheme = 1, scale = 0, select = 3, rug = 0, xlab = "", ylab = bquote(beta[2]))
dev.off()

pdf("plots/m2_b5.pdf", width = 7, height = 5)
plot(m2, scheme = 1, scale = 0, select = 6, rug = 0, xlab = "", ylab = bquote(beta[5]))
dev.off()

pdf("plots/m2_b6.pdf", width = 7, height = 5)
plot(m2, scheme = 1, scale = 0, select = 7, rug = 0, xlab = "", ylab = bquote(beta[6]))
dev.off()
