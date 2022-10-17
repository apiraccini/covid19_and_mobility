rm(list = ls())
set.seed(42)
opar <- par()


# libraries --------------------------------------------------------------------

library(fda)
library(refund)
library(tidyverse)


# utils ------------------------------------------------------------------------

load("data/finaldata.RData")


# old clustering ---------------------------------------------------------------

# based on betas from basis repr.

dim(cases$fd$coefs)
betas <- t(cases$fd$coefs)
dimnames(betas)[[1]] <- unique(d$stato)

# la variabilita'/caratt.individuali delle curve sono racchiuse nei beta, 
# che moltiplicano funz tutte uguali
# -> si puo' fare clustering con un semplice k-means sui beta
# l'unico problema di questo approccio e' non considerare la variabilita' dei 
# coeff stimati (pazienza)

set.seed(42)
tryk <- 1:26
withinss <- numeric(length(tryk))
for (i in seq_along(tryk)){
  cat("Trying with", tryk[i], "centers...\n")
  km <- kmeans(betas, centers = tryk[i], nstart = 200, iter.max = 1000)
  withinss[i] <- km$withinss
}
plot(tryk, withinss, main = "Within SS", xlab = "k", ylab = "", type = "b")

set.seed(42)
km1 <- kmeans(betas, centers = 4, nstart = 200, iter.max = 1000)
matplot(cases$mat.sm, type = "l", col = km1$cluster, lty = 1, xlab = "")
legend("bottomright", lty = 1, col = 1:4, legend = paste0("group ",1:4), bty = "n", cex = 0.8)

for (j in 1:4) cat("Cluster ", j, ": ", levels(d$stato)[km1$cluster == j], "\n")

hc <- hclust(dist(betas), method = "complete")
plot(hc)
rect.hclust(hc, k = 4)
hc1 <- cutree(hc, k = 4)

kmeans1 <- data.frame(name = names(km1$cluster), cluster = km1$cluster, row.names = NULL)

plot_kmeans <- eu %>% 
  merge(kmeans1, by = "name") %>% 
  mutate(cluster = as.factor(cluster)) %>% 
  ggplot(aes(fill = cluster)) + 
  geom_sf() +
  scale_x_continuous(limits = c(-10, 30)) +
  scale_y_continuous(limits = c(35, 65)) +
  scale_fill_brewer(palette = "Set3") + 
  labs(title = "K-means")
plot_kmeans  

hclust1 <- data.frame(name = names(hc1), cluster2 = hc1, row.names = NULL)
hclust1$cluster[hclust1$cluster2 == 3] <- 4
hclust1$cluster[hclust1$cluster2 == 4] <- 2
hclust1$cluster[hclust1$cluster2 == 2] <- 1
hclust1$cluster[hclust1$cluster2 == 1] <- 3

plot_hclust <- eu %>% 
  merge(hclust1, by = "name") %>% 
  mutate(cluster = as.factor(cluster)) %>% 
  ggplot(aes(fill = cluster)) + 
  geom_sf() +
  scale_x_continuous(limits = c(-10, 30)) +
  scale_y_continuous(limits = c(35, 65)) +
  scale_fill_brewer(palette = "Set3") + 
  labs(title = "Hclust")
plot_hclust 

ggpubr::ggarrange(plot_kmeans, plot_hclust)


# PDA --------------------------------------------------------------------------


range<-cases$fd$basis$rangeval
nbasis<-20
pdabasisfd <- create.bspline.basis(range, nbasis)
betafdPar  <- fdPar(fd(matrix(0,nbasis,1),pdabasisfd))

#  fissiamo la lista di oggetti per i parametri funzionali per i pesi
bwtlist = vector("list", 2)
bwtlist[[1]] <- betafdPar
bwtlist[[2]] <- betafdPar
#
# fissiamo la lista con i dati in formato fd
xfdlist <- list(cases$fd)

#analisi differenziale principale
pdaList <- pda.fd(xfdlist, bwtlist)
#
#parametri stimati
bwtestlist <- pdaList$bwtlist

#grafici
par(mfrow=c(2,1),pty="m")
day<-seq(1,813, len=813)
#
for (j in 1:2) {
  bfdParj <- bwtestlist[[j]]
  bvals = eval.fd(day,bwtestlist[[j]]$fd)
  plot(day,bvals,type='l', xlab="tempo", ylab=paste("beta",j-1))
  abline(h=0, lty=2)
}

par(mfrow=c(1,1))

#calcolo discriminante
dfd = 0.25*pdaList$bwtlist[[2]]$fd^2 - pdaList$bwtlist[[1]]$fd 
plot(dfd, ylab="discriminante", xlab="tempo")


#diagramma di biforcazione
#
pda.overlay.new<-function (pdaList, nfine = 501, ncoarse = 11, ...) 
{
  fdrange = pdaList$bwtlist[[1]]$fd$basis$range
  tfine = seq(fdrange[1], fdrange[2], len = nfine)
  beta0vals = eval.fd(tfine, pdaList$bwtlist[[1]]$fd, 0)
  beta1vals = eval.fd(tfine, pdaList$bwtlist[[2]]$fd, 0)
  plot(-beta1vals, -beta0vals, type = "l", col = 4,xlab=expression(beta[1]), ylab=expression(beta[0]) , ...)
  abline(h = 0, col = 2)
  abline(v = 0, col = 2)
  bv = seq(min(-beta1vals), max(-beta1vals), len = nfine)
  lines(bv, -(bv/2)^2, col = 2, lty = 2, ...)
  tcoarse = seq(fdrange[1], fdrange[2], len = ncoarse)
  beta0valsc = eval.fd(tcoarse, pdaList$bwtlist[[1]]$fd, 0)
  beta1valsc = eval.fd(tcoarse, pdaList$bwtlist[[2]]$fd, 0)
  text(x = -beta1valsc, y = -beta0valsc, labels = tcoarse, col = 4)
}

pda.overlay.new(pdaList,lwd=2,cex.lab=1.5,cex.axis=1.5)

#Analisi dei residui

#  calcoliamo le funzioni derivate 
#
Lfdest <- Lfd(2, bwtestlist)
#
#valutazione di Lx_i
force <- eval.fd(day, cases$fd, Lfdest)
#
#accelerazione, derivata seconda, per i dati originali (lisciati)
cases_accel     <- eval.fd(day, cases$fd, 2)
#accelerazione media
cases_meanaccel <- apply(cases_accel, 1, mean)
#
par(mfrow=c(2,1))
matplot(day, force, type="l", lty=1)
matplot(day, cases_accel, type="l", lty=1)
#
#oppure conforntiamo con la accelerazione media
par(mfrow=c(1,1))
yrange <- c(min(min(cases_meanaccel),min(force)),
            max(max(cases_meanaccel),max(force)))
matplot(day, force, type="l", lty=1, ylim=yrange)
#
lines(day, cases_meanaccel, lty=4, lwd=2) #non si vede manco il cazzo

# confrontiamo le due funzioni medie
forcemean <- apply(force, 1, mean)
#
plot(day, forcemean, type="l", lty=1, ylim=c(-0.05,0.05))
#
lines(day, cases_meanaccel, lty=4, col=2)

#calcolo di R2
resmat1 = eval.fd(day,pdaList$resfdlist[[1]])
resmat0 = cases_accel

SSE0 = apply((resmat0)^2, 1, sum)
SSE1 = apply(resmat1^2, 1, sum)
cases.R2 = (SSE0-SSE1)/SSE0
plot(day, cases.R2, lwd=2, xlab="tempo", type="l", ylab=expression(paste(R^2,"(t)")), col=4)

#  solve equation
result <- odesolv(bwtestlist)
xp <- result[[1]]
yp <- result[[2]]

#  plot the two solutions
par(mfrow=c(2,1),pty="m")
pltrng <- c(min(yp[1,,]), max(yp[1,,]))
matplot(xp,t(yp[1,,]), type="l", lty=1, ylim=pltrng, main="Funzione")
abline(h=0, lty=2)
#
pltrng <- c(min(yp[2,,]), max(yp[2,,]))
matplot(xp,t(yp[2,,]), type="l", lty=1, ylim=pltrng, main="Derivata")
abline(h=0, lty=2)
#


#  plot the two solutions
par(mfrow=c(2,1),pty="m")
plot(xp,t(yp[1,1,]), type="l", lty=1, lwd=2, main="Funzione", ylab=expression(xi[1](t)), xlab="tempo")
abline(h=0, lty=2)
plot(xp,t(yp[1,2,]), type="l", lty=1, lwd=2, main="Funzione", ylab=expression(xi[2](t)), xlab="tempo")
abline(h=0, lty=2)
#
#
par(mfrow=c(2,1),pty="m")
plot(xp,t(yp[2,1,]), type="l", lty=1, main="Derivata", ylab=expression(paste("D", xi[1](t))), xlab="tempo (secondi)", lwd=2)
abline(h=0, lty=2)
#
plot(xp,t(yp[2,2,]), type="l", lty=1, main="Derivata", ylab=expression(paste("D", xi[2](t))), xlab="tempo (secondi)", lwd=2)
abline(h=0, lty=2)


# SIR model --------------------------------------------------------------------

casi <- as.numeric(unlist(d[d$stato=="Italy","cases"]))
morti <- as.numeric(unlist(d[d$stato=="Italy","deaths"]))
length(casi)
n <- 6*10^7

#title: "Il modello SIR e sue estensioni - Studi di simulazione"
#subtitle: "Statistica medica ed epidemiologia Progredito (a.a. 2020/21)"
#date: "5 maggio 2022"
#author: "Prof.sse A.R. Brazzale & G. Boccuzzo"

#Load required packages
if("deSolve" %in% rownames(installed.packages()) == FALSE) {install.packages("deSolve", repos='http://cran.stat.unipd.it')}
library(deSolve)

set.seed(21)

# Section 1: First-hand exploration of a SIR model
# transmission rate 
beta <- 1/7 # tasso di infezione, tempo di latenza (circa 5gg prima che la malattia si presenti)
# recovery rate
gamma <- 1/21 #tasso di guarigione, tempo medio di attesa per la guarigione-valore atteso di una esponenziale
# basic reproduction number
R0 <- beta/gamma
# numero medio di soggetti che vengono contagiati ad inizio epidemia (circa 2)
R0

# 1. Ordinary differential equations (ODE)
# vogliamo risolvere per diversi valori di t 
# l'equazione differenziale=consistenza numerica dei comparti

# Scegliamo t=day
# time grid

time <- 1:length(casi)

# Modello deterministico 
# SIR model

# funzione che restituisce la proporzione per comparto per ogni t
# data la composizione iniziale della popolazione e i parametri
# che caratterizzano l'epidemia
n <- 6*10^5
casi <- cases$mat.raw[1,1]#dovrei prendere n e casi corrispondenti con n popolazione della nazione
pop0 <- c(S=1-casi[1]/n, I=casi[1]/n) 

odeSIR <- function(t, states, param)  
{
  with(as.list(c(states, param)),
       { 
         dS = - beta * S * I
         dI = beta * S * I - gamma * I 
         list(c(dS, dI))
       })
}

# serve il vettore dei tempi per i quali cerchiamo la soluzione
# ODE solver
solSIR <- ode(y = pop0, times = 1:190, method="euler",
              func = odeSIR, parms = c(beta, gamma)) 

?ode

str(solSIR) 
head(solSIR)

# sample paths of SIR model
St_ode <- solSIR[,2]
It_ode <- solSIR[,3]
Rt_ode <- 1 - St_ode - It_ode 

# epidemic model
plot(St_ode, type="l", col="green", xlab="time", ylab="population",
     ylim=c(0,1), main="sample paths of SIR model, with ODE")

lines(It_ode, col="red")
lines(Rt_ode, col="blue")

# questa seconda funzione è presa da internet

sir_1 <- function(beta, gamma, S0, I0, R0, times) {
  require(deSolve) # for the "ode" function
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS <- -beta * I * S
      dI <-  beta * I * S - gamma * I
      dR <-  gamma * I
      return(list(c(dS, dI, dR)))
    })
  }
  
  # the parameters values:
  parameters_values <- c(beta  = beta, gamma = gamma)
  
  # the initial values of variables:
  initial_values <- c(S = S0, I = I0, R = R0)
  
  # solving
  out <- ode(initial_values, times, sir_equations, parameters_values)
  ?ode.2D
  
  # returning the output:
  as.data.frame(out)
}

sir_1(beta = 0.004, gamma = 0.5, S0 = 999, I0 = 1, R0 = 0, times = seq(0, 10))


# positive smoothing ------------------------------------------------------

dim(cases$mat.raw)

library(fda)

m <- exp(deaths$mat.raw)-1
loglambdaseq =seq(3, 500, len = 20)[1:3]
#positive smooth and evaluate
get_fd <- function(m, loglambdaseq = NULL){
  
  # imputazione na per le date iniziali (si assume assenza di rilevazioni)
  for (i in 1:nrow(m)) m[i, 1:min(which(!is.na(m[i,])))] <- 0
  
  # parametri per il lisciamento
  day <- seq(1, ncol(m), by = 7)
  M <- 6
  k <- length(day) + M - 2
  base <- create.bspline.basis(range(day), k, M, day)
  lfd.degree <- 2
  
  # penalita'
  loglambda <- seq(3, 5, len = 20)
  if(!is.null(loglambdaseq)) loglambda <- loglambdaseq
  metrics <- matrix(NA, length(loglambda), 3)
  colnames(metrics) <- c("lambda", "df", "gcv")
  betas <- array(NA, dim = c(nrow(m), k, length(loglambda)))
  # ciclo di smoothing
  for (j in seq_along(loglambda)){
    
    cat("Using lambda ", j, "/", length(loglambda), "...\n")
    fdPar <- fdPar(base, lfd.degree, lambda = 10)
    
    metrics2 <- matrix(NA, nrow(m), 2)
    
    # una serie per volta
    for(i in 1:nrow(m)){
      
      #if(i %% 5 == 0) cat("\tSmoothing series", i, "/", nrow(m), "...\n")
      s <- as.numeric(na.omit(m[i,]))
      ids <- which(!is.na(m[i,]))
      # le successive due righe calcolano la media mobile della serie
      s <- as.numeric(na.omit(forecast::ma(s, 6)))
      ids <- ids[4:(length(ids)-3)]
      
      sb <- smooth.pos(ids, s, fdPar)
      
      metrics2[i,] <- c(10, sb$Flist$f)
      betas[i,,j] <- sb$Wfdobj$coefs
      
    }
    
    # calcolo metriche
    metrics[j,] <- c(loglambda[j], apply(metrics2, 2, mean))
    cat(j)
    cat(i)
  }
  str(sb)
  # visualizzazione metriche di errore
  ##plot(metrics[,1], metrics[,3], xlab = expression(log(lambda)), ylab = "GCV", type = "b")
  #abline(v = metrics[which.min(metrics[,3]) ,1], col = 2, lty = 2)
  #plot(metrics[,1], metrics[,2], xlab = expression(log(lambda)), ylab = "df", type = "b")
  
  # output
  betabest <- betas[,,which.min(metrics[,3])]
  myfd <- fd(coef = t(betabest), basisobj = base)
  out <- list(fd = myfd, metrics = metrics)
  
  out
}

out <- get_fd(m)
y_out = exp(eval.fd(grid1, out$fd))-1
plot(grid1, y_out, type="l", ylim = c(-1,9*10^5), ylab="value")


# model based clustering --------------------------------------------------


library(fdapace)

d <- data.frame(d)
d$date <- as.Date(d$date, format = "%Y-%m-%d")
d11 <- d
d <- d11
d <- d[1:6000,]
Flies <- MakeFPCAInputs(d$stato, d$date, d$cases) 

data(medfly25)
str(medfly25)
Flies <- MakeFPCAInputs(medfly25$ID, medfly25$Days, medfly25$nEggs)
library(EMCluster)
newClust <- FClust(Flies$Ly, Flies$Lt, k = 4, optnsFPCA = 
                     list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.70))

newClust <- FClust(Flies$Ly, Flies$Lt, k = 4, optnsFPCA = 
                     list(methodMuCovEst = 'smooth', userBwCov = 2, FVEthreshold = 0.70))

str(newClust)
