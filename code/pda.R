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


# PDA --------------------------------------------------------------------------

range <- cases$fd$basis$rangeval
nbasis <- 20
pdabasisfd <- create.bspline.basis(range, nbasis)
betafdPar <- fdPar(fd(matrix(0,nbasis,1),pdabasisfd))

# fissiamo la lista di oggetti per i parametri funzionali per i pesi
bwtlist <- vector("list", 2)
bwtlist[[1]] <- betafdPar
bwtlist[[2]] <- betafdPar

# fissiamo la lista con i dati in formato fd
xfdlist <- list(cases$fd)

# analisi differenziale principale
pdaList <- pda.fd(xfdlist, bwtlist)

# parametri stimati
bwtestlist <- pdaList$bwtlist

# grafici dei coefficienti
day <- grid1

bvals <- eval.fd(day,bwtestlist[[1]]$fd)
plot(day, bvals, type = 'l', xlab = "tempo", ylab = bquote(beta[0]), lwd = 2)
abline(h = 0, lty = 2, col = 2)

bvals <- eval.fd(day,bwtestlist[[2]]$fd)
plot(day, bvals, type = 'l', xlab = "tempo", ylab = bquote(beta[1]), lwd = 2)
abline(h = 0, lty = 2, col = 2)

# x report
pdf("plots/pda_beta0.pdf", width = 7, height = 5)
bvals <- eval.fd(day,bwtestlist[[1]]$fd)
plot(day, bvals, type = 'l', xlab = "tempo", ylab = bquote(beta[0]), lwd = 2)
abline(h = 0, lty = 2, col = 2)
dev.off()

pdf("plots/pda_beta1.pdf", width = 7, height = 5)
bvals <- eval.fd(day,bwtestlist[[2]]$fd)
plot(day, bvals, type = 'l', xlab = "tempo", ylab = bquote(beta[1]), lwd = 2)
abline(h = 0, lty = 2, col = 2)
dev.off()


# confronto tra accelerazione e funzionale stimato

Lfdest <- Lfd(2, bwtestlist)
funz <- eval.fd(day, cases$fd, Lfdest)

accel  <- eval.fd(day, cases$fd, 2)

matplot(grid1, accel, type = "l", col = "lightgrey", lty = 1, ylim = c(-0.2, 0.2), ylab = expression(D^2~(x)), xlab = "tempo")
abline(h = 0, col = 1, lty = 2)
matlines(grid1, apply(accel, 1, mean), col = "red", lwd = 2)

matplot(grid1, funz, type = "l", col = "lightgrey", lty = 1, ylim = c(-0.2, 0.2), ylab = bquote(L(x)), xlab = "tempo")
abline(h = 0, col = 1, lty = 2)
matlines(grid1, apply(funz, 1, mean), col = "red", lwd = 2)

# x report
pdf("plots/pda_accel.pdf", width = 7, height = 5)
matplot(grid1, accel, type = "l", col = "lightgrey", lty = 1, ylim = c(-0.2, 0.2), ylab = bquote(D^2(x)), xlab = "tempo")
abline(h = 0, col = 1, lty = 2)
matlines(grid1, apply(accel, 1, mean), col = "red", lwd = 2)
dev.off()

pdf("plots/pda_funz.pdf", width = 7, height = 5)
matplot(grid1, funz, type = "l", col = "lightgrey", lty = 1, ylim = c(-0.2, 0.2), ylab = bquote(L(x)), xlab = "tempo")
abline(h = 0, col = 1, lty = 2)
matlines(grid1, apply(funz, 1, mean), col = "red", lwd = 2)
dev.off()