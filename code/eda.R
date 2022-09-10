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

proj_dir <- "~/R_Jobs/DatiFunzionali/progetto_funzionali"
if (normalizePath(getwd()) != normalizePath(proj_dir)) setwd(proj_dir)


# setup ------------------------------------------------------------------------

load("data/smoothed.RData")

# griglie
grid1 <- seq(1, length(unique(d$date)), by = 7)
grid2 <- seq(1, length(unique(d$date)), by = 30)

# dataset per refund
dd <- data.frame(stato = unique(d$stato))


# cases ------------------------------------------------------------------------

#GCV
pdf("plots/gcv_casi.pdf", width =7, height =5 )
plot(cases$metrics[,1], cases$metrics[,3], xlab = expression(log(lambda)), ylab = "GCV", type = "b")
abline(v = cases$metrics[which.min(cases$metrics[,3]) ,1], col = 2, lty = 2)
dev.off()

pdf("plots/fda_casi.pdf", width =7, height =5 )
plot(cases$fd, xlab="tempo",ylab="log casi")
dev.off()
cases$mat.sm <- eval.fd(grid1, cases$fd)

# media
pdf("plots/mean_casi.pdf", width =7, height =5 )
matplot(grid1, cases$mat.sm, type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "")
matlines(grid1, apply(cases$mat.sm, 1, mean), col = "red", lwd = 2)
dev.off()

# covarianza
cov.fd <- var.fd(cases$fd)

g <- seq(grid1[1], grid1[length(grid1)], len = 100)
cov.sm <- eval.bifd(g, g, cov.fd)
#persp(g, g, cov.sm, theta = -45, phi = 25, r = 3, expand = 0.5, ticktype = "detailed",
#      xlab = "", ylab = "", zlab = "", main = "Covarianza cases")

pdf("plots/cov_casi.pdf", width =7, height =5 )
image(g, g, cov.sm, xlab = "", ylab = "", main = "", col=hcl.colors(50,palette="Blues", rev=T))
contour(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza casi", add=T)
dev.off()

#boxplot
pdf("plots/boxplot_casi.pdf", width =7, height =5 )
boxplot(cases$fd, xlab="tempo", ylab="")
dev.off()

# fpca
ncomp <- 4
cases$pca <- pca.fd(cases$fd, ncomp)
#print(cases$pca)$harmonics

#par(mfrow = c(2,2), mar = rep(2, 4))
#plot(cases$pca, cex.main = 0.9)
#par(opar)

pdf("plots/pca_casi.pdf", width =7, height =5 )
plot(cases$pca$harmonics, xlab="tempo", ylab="")
dev.off()

names(cases)
dd$cases <- t(cases$mat.sm)


# deaths -----------------------------------------------------------------------

#GCV
pdf("plots/gcv_decessi.pdf", width =7, height =5 )
plot(deaths$metrics[,1], deaths$metrics[,3], xlab = expression(log(lambda)), ylab = "GCV", type = "b")
abline(v = deaths$metrics[which.min(deaths$metrics[,3]) ,1], col = 2, lty = 2)
dev.off()

pdf("plots/fda_decessi.pdf", width =7, height =5 )
plot(deaths$fd, xlab="tempo", ylab="log decessi")
dev.off()

deaths$mat.sm <- eval.fd(grid1, deaths$fd)


# media
pdf("plots/mean_decessi.pdf", width =7, height =5 )
matplot(grid1, deaths$mat.sm, type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "")
matlines(grid1, apply(deaths$mat.sm, 1, mean), col = "red", lwd = 2)
dev.off()

# covarianza
cov.fd <- var.fd(deaths$fd)

cov.sm <- eval.bifd(g, g, cov.fd)
#persp(g, g, cov.sm, theta = -45, phi = 25, r = 3, expand = 0.5, ticktype = "detailed",
#     xlab = "", ylab = "", zlab = "", main = "Covarianza deaths")

pdf("plots/cov_decessi.pdf", width =7, height =5 )
image(g, g, cov.sm, xlab = "", ylab = "", main = "", col=hcl.colors(50,palette="Blues", rev=T))
contour(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza decessi", add=T)
dev.off()

#boxplot
pdf("plots/boxplot_decessi.pdf", width =7, height =5 )
boxplot(deaths$fd, xlab="tempo", ylab="")
dev.off()

# fpca
ncomp <- 4
deaths$pca <- pca.fd(deaths$fd, ncomp)
#print(deaths$pca)$harmonics

#par(mfrow = c(2,2), mar = rep(2, 4))
#plot(deaths$pca, cex.main = 0.9)
par(opar)

pdf("plots/pca_decessi.pdf", width =7, height =5 )
plot(deaths$pca$harmonics,xlab="tempo", ylab="")
dev.off()

names(deaths)
dd$deaths <- t(deaths$mat.sm)

# transit --------------------------------------------------------

plot(transit$fd)
transit$mat.sm <- eval.fd(grid1, transit$fd)

# media
matplot(transit$mat.sm, type = "l", col = "lightgrey", main = "transit (smoothed)")
matlines(apply(transit$mat.sm, 1, mean), col = "red", lwd = 2)

# covarianza
cov.fd <- var.fd(transit$fd)

cov.sm <- eval.bifd(g, g, cov.fd)
persp(g, g, cov.sm, theta = -45, phi = 25, r = 3, expand = 0.5, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "", main = "Covarianza transit")
contour(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza transit")
image(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza transit")

# fpca
ncomp <- 4
transit$pca <- pca.fd(transit$fd, ncomp)
#print(retail_and_recreation$pca)$harmonics

par(mfrow = c(2,2), mar = rep(2, 4))
plot(transit$pca, cex.main = 0.9)
par(opar)
plot(transit$pca$harmonics)

names(transit)
dd$transit <- t(transit$mat.sm)

# x report
pdf("plots/transit1.pdf", width = 7, height = 5)
plot(transit$metrics[,1], transit$metrics[,3], xlab = expression(log(lambda)), ylab = "GCV", type = "b")
abline(v =transit$metrics[which.min(transit$metrics[,3]) ,1], col = 2, lty = 2)
dev.off()

pdf("plots/transit2.pdf", width = 7, height = 5)  
matplot(grid1,transit$mat.sm, type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "")
matlines(grid1, apply(transit$mat.sm, 1, mean), col = "red", lwd = 2)
dev.off()

# retail_and_recreation --------------------------------------------------------

plot(retail_and_recreation$fd)
retail_and_recreation$mat.sm <- eval.fd(grid1, retail_and_recreation$fd)

# media
matplot(retail_and_recreation$mat.sm, type = "l", col = "lightgrey", main = "retail_and_recreation (smoothed)")
matlines(apply(retail_and_recreation$mat.sm, 1, mean), col = "red", lwd = 2)

# covarianza
cov.fd <- var.fd(retail_and_recreation$fd)

cov.sm <- eval.bifd(g, g, cov.fd)
persp(g, g, cov.sm, theta = -45, phi = 25, r = 3, expand = 0.5, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "", main = "Covarianza retail_and_recreation")
contour(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza retail_and_recreation")
image(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza retail_and_recreation")

# fpca
ncomp <- 4
retail_and_recreation$pca <- pca.fd(retail_and_recreation$fd, ncomp)
#print(retail_and_recreation$pca)$harmonics

par(mfrow = c(2,2), mar = rep(2, 4))
plot(retail_and_recreation$pca, cex.main = 0.9)
par(opar)
plot(retail_and_recreation$pca$harmonics)

names(retail_and_recreation)
dd$retail_and_recreation <- t(retail_and_recreation$mat.sm)

# x report
pdf("plots/retail_and_recreation1.pdf", width = 7, height = 5)
plot(retail_and_recreation$metrics[,1], retail_and_recreation$metrics[,3], xlab = expression(log(lambda)), ylab = "GCV", type = "b")
abline(v = retail_and_recreation$metrics[which.min(retail_and_recreation$metrics[,3]) ,1], col = 2, lty = 2)
dev.off()

pdf("plots/retail_and_recreation2.pdf", width = 7, height = 5)  
matplot(grid1, retail_and_recreation$mat.sm, type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "")
matlines(grid1, apply(retail_and_recreation$mat.sm, 1, mean), col = "red", lwd = 2)
dev.off()


# grocery_and_pharmacy ---------------------------------------------------------

plot(grocery_and_pharmacy$fd)
grocery_and_pharmacy$mat.sm <- eval.fd(grid1, grocery_and_pharmacy$fd)

# media
matplot(grocery_and_pharmacy$mat.sm, type = "l", col = "lightgrey", main = "grocery_and_pharmacy (smoothed)")
matlines(apply(grocery_and_pharmacy$mat.sm, 1, mean), col = "red", lwd = 2)

# covarianza
cov.fd <- var.fd(grocery_and_pharmacy$fd)

cov.sm <- eval.bifd(g, g, cov.fd)
persp(g, g, cov.sm, theta = -45, phi = 25, r = 3, expand = 0.5, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "", main = "Covarianza grocery_and_pharmacy")
contour(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza grocery_and_pharmacy")
image(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza grocery_and_pharmacy")

# fpca
ncomp <- 4
grocery_and_pharmacy$pca <- pca.fd(grocery_and_pharmacy$fd, ncomp)
#print(grocery_and_pharmacy$pca)$harmonics

par(mfrow = c(2,2), mar = rep(2, 4))
plot(grocery_and_pharmacy$pca, cex.main = 0.9)
par(opar)
plot(grocery_and_pharmacy$pca$harmonics)

names(grocery_and_pharmacy)
dd$grocery_and_pharmacy <- t(grocery_and_pharmacy$mat.sm)

# x report
pdf("plots/grocery_and_pharmacy1.pdf", width = 7, height = 5)
plot(grocery_and_pharmacy$metrics[,1], grocery_and_pharmacy$metrics[,3], xlab = expression(log(lambda)), ylab = "GCV", type = "b")
abline(v = grocery_and_pharmacy$metrics[which.min(grocery_and_pharmacy$metrics[,3]) ,1], col = 2, lty = 2)
dev.off()

pdf("plots/grocery_and_pharmacy2.pdf", width = 7, height = 5)  
matplot(grid1, grocery_and_pharmacy$mat.sm, type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "")
matlines(grid1, apply(grocery_and_pharmacy$mat.sm, 1, mean), col = "red", lwd = 2)
dev.off()


# parks ------------------------------------------------------------------------

plot(parks$fd)
parks$mat.sm <- eval.fd(grid1, parks$fd)

# media
matplot(parks$mat.sm, type = "l", col = "lightgrey", main = "parks (smoothed)")
matlines(apply(parks$mat.sm, 1, mean), col = "red", lwd = 2)

# covarianza
cov.fd <- var.fd(parks$fd)

cov.sm <- eval.bifd(g, g, cov.fd)
persp(g, g, cov.sm, theta = -45, phi = 25, r = 3, expand = 0.5, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "", main = "Covarianza parks")
contour(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza parks")
image(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza parks")

# fpca
ncomp <- 4
parks$pca <- pca.fd(parks$fd, ncomp)
#print(parks$pca)$harmonics

par(mfrow = c(2,2), mar = rep(2, 4))
plot(parks$pca, cex.main = 0.9)
par(opar)
plot(parks$pca$harmonics)

names(parks)
dd$parks <- t(parks$mat.sm)

# x report
pdf("plots/parks1.pdf", width = 7, height = 5)
plot(parks$metrics[,1], parks$metrics[,3], xlab = expression(log(lambda)), ylab = "GCV", type = "b")
abline(v = parks$metrics[which.min(parks$metrics[,3]) ,1], col = 2, lty = 2)
dev.off()

pdf("plots/parks2.pdf", width = 7, height = 5)  
matplot(grid1, parks$mat.sm, type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "")
matlines(grid1, apply(parks$mat.sm, 1, mean), col = "red", lwd = 2)
dev.off()


# workplaces -------------------------------------------------------------------

plot(workplaces$fd)
workplaces$mat.sm <- eval.fd(grid1, workplaces$fd)

# media
matplot(workplaces$mat.sm, type = "l", col = "lightgrey", main = "workplaces (smoothed)")
matlines(apply(workplaces$mat.sm, 1, mean), col = "red", lwd = 2)

# covarianza
cov.fd <- var.fd(workplaces$fd)

cov.sm <- eval.bifd(g, g, cov.fd)
persp(g, g, cov.sm, theta = -45, phi = 25, r = 3, expand = 0.5, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "", main = "Covarianza workplaces")
contour(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza workplaces")
image(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza workplaces")

# fpca
ncomp <- 4
workplaces$pca <- pca.fd(workplaces$fd, ncomp)
#print(workplaces$pca)$harmonics

par(mfrow = c(2,2), mar = rep(2, 4))
plot(workplaces$pca, cex.main = 0.9)
par(opar)
plot(workplaces$pca$harmonics)

names(workplaces)
dd$workplaces <- t(workplaces$mat.sm)

# x report
pdf("plots/workplaces1.pdf", width = 7, height = 5)
plot(workplaces$metrics[,1], workplaces$metrics[,3], xlab = expression(log(lambda)), ylab = "GCV", type = "b")
abline(v = workplaces$metrics[which.min(workplaces$metrics[,3]) ,1], col = 2, lty = 2)
dev.off()

pdf("plots/workplaces2.pdf", width = 7, height = 5)  
matplot(grid1, workplaces$mat.sm, type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "")
matlines(grid1, apply(workplaces$mat.sm, 1, mean), col = "red", lwd = 2)
dev.off()


# residential ------------------------------------------------------------------

plot(residential$fd)
residential$mat.sm <- eval.fd(grid1, residential$fd)

# media
matplot(residential$mat.sm, type = "l", col = "lightgrey", main = "residential (smoothed)")
matlines(apply(residential$mat.sm, 1, mean), col = "red", lwd = 2)

# covarianza
cov.fd <- var.fd(residential$fd)

cov.sm <- eval.bifd(g, g, cov.fd)
persp(g, g, cov.sm, theta = -45, phi = 25, r = 3, expand = 0.5, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "", main = "Covarianza residential")
contour(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza residential")
image(g, g, cov.sm, xlab = "", ylab = "", main = "Covarianza residential")

# fpca
ncomp <- 4
residential$pca <- pca.fd(residential$fd, ncomp)
#print(residential$pca)$harmonics

par(mfrow = c(2,2), mar = rep(2, 4))
plot(residential$pca, cex.main = 0.9)
par(opar)
plot(residential$pca$harmonics)

names(residential)
dd$residential <- t(residential$mat.sm)

# x report
pdf("plots/residential1.pdf", width = 7, height = 5)
plot(residential$metrics[,1], residential$metrics[,3], xlab = expression(log(lambda)), ylab = "GCV", type = "b")
abline(v = residential$metrics[which.min(residential$metrics[,3]) ,1], col = 2, lty = 2)
dev.off()

pdf("plots/residential2.pdf", width = 7, height = 5)  
matplot(grid1, residential$mat.sm, type = "l", col = "lightgrey", lty = 1, xlab = "tempo", ylab = "")
matlines(grid1, apply(residential$mat.sm, 1, mean), col = "red", lwd = 2)
dev.off()


# salvataggio ------------------------------------------------------------------

save(d, dd, grid1, grid2, cases, deaths, grocery_and_pharmacy, parks,
     residential, retail_and_recreation, transit, workplaces, file = "data/finaldata.RData")