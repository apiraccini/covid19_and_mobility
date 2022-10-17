### PROGETTO PER IL CORSO DI DATI FUNZIONALI
# Francesco Martella, Sara Meneghetti, Alessio Piraccini

rm(list = ls())
set.seed(42)


# libraries --------------------------------------------------------------------

library(fda)
library(refund)
library(tidyverse)


# data loading -----------------------------------------------------------------

mobility <- read.csv("data/Global_Mobility_Report.csv", header = T) # was too big to upload
covid_eu <- read.csv("data/51DCR-Ln.csv")


# preprocessing ----------------------------------------------------------------

# check
str(mobility)
str(covid_eu)

# diversi punti di rilevazione per ogni stato (a volte con valori mancanti)
mobility %>% 
  group_by(country_region) %>% summarise(k = length(unique(place_id)))

# nazioni europee
eu_nations <- unique(covid_eu$countriesAndTerritories)

# manipolazione mobility
mobility_eu <- mobility %>% 
  filter(country_region %in% eu_nations) %>% 
  select(date, country_region, place_id, retail_and_recreation_percent_change_from_baseline,
         grocery_and_pharmacy_percent_change_from_baseline, parks_percent_change_from_baseline,
         transit_stations_percent_change_from_baseline, workplaces_percent_change_from_baseline,
         residential_percent_change_from_baseline) %>% 
  rename(stato = country_region,
         retail_and_recreation_pcfb = retail_and_recreation_percent_change_from_baseline,
         grocery_and_pharmacy_pcfb = grocery_and_pharmacy_percent_change_from_baseline,
         parks_pcfb = parks_percent_change_from_baseline,
         transit_pcfb = transit_stations_percent_change_from_baseline,
         workplaces_pcfb = workplaces_percent_change_from_baseline,
         residential_pcfb = residential_percent_change_from_baseline) %>% 
  group_by(stato, date) %>% 
  summarise(retail_and_recreation_pcfb = mean(retail_and_recreation_pcfb, na.rm = T),
            grocery_and_pharmacy_pcfb = mean(grocery_and_pharmacy_pcfb, na.rm = T),
            parks_pcfb = mean(parks_pcfb, na.rm = T),
            transit_pcfb = mean(transit_pcfb, na.rm = T),
            workplaces_pcfb = mean(workplaces_pcfb, na.rm = T),
            residential_pcfb = mean(residential_pcfb, na.rm = T)) %>% 
  arrange(stato, date)

# manipolazione covid
covid_eu <- covid_eu %>% 
  select(dateRep, countriesAndTerritories, cases, deaths, popData2020) %>% 
  rename(date = dateRep, stato = countriesAndTerritories, pop2020 = popData2020)

# check
str(mobility_eu)
str(covid_eu)

summary(mobility_eu)
summary(covid_eu)

# per ora
mobility_eu <- na.omit(mobility_eu)
covid_eu <- na.omit(covid_eu)

str(mobility_eu)
str(covid_eu)


# merging ----------------------------------------------------------------------

mobility_eu$date <- as.Date(mobility_eu$date, format = "%Y-%m-%d")
covid_eu$date <- as.Date(covid_eu$date, format = "%d/%m/%Y")

range(covid_eu$date)
range(mobility_eu$date)

d <- merge(mobility_eu, covid_eu, by = c("date", "stato") , all.x = T) %>% as_tibble()
str(d)
d$stato <- as.factor(d$stato)
summary(d) # na sono giorni non corrispondenti
d <- na.omit(d)

table(d$cases<0)
d$cases[d$cases<0] <- 0

table(d$deaths<0)
d$deaths[d$deaths<0] <- 0

save(d, file = "data/d.RData")
rm(list = ls())
gc()

load("data/d.RData")


# visualisations ---------------------------------------------------------------

plotf <- function(var = "deaths"){
  p <- ggplot(d, aes(x = date, y = .data[[var]], colour = stato)) +
    geom_line(alpha = 0.8) +
    labs(x = "", y = "", title = var) +
    theme(legend.position = "none")
  p
}

p <- lapply(setdiff(names(d), c("date", "stato", "pop2020")), plotf)
ggpubr::ggarrange(plotlist = p, nrow = 2)


# functional data --------------------------------------------------------------

# idea: ogni covariata dentro una lista, dove teniamo matrice grezza, matrice 
# lisciata, dato in formato fd, ecc.. (cosi non creiamo mille oggetti)

# dati funzionali in formato matriciale
mat.raw <- d %>% select(date, stato, deaths) %>% 
  pivot_wider(names_from = "date", values_from = "deaths") %>% select(-stato) %>% as.matrix()
mat.raw[!is.na(mat.raw)] <- log(mat.raw[!is.na(mat.raw)] + 1) # log(1+x) per deaths e cases
deaths <- list(mat.raw = mat.raw)

mat.raw <- d %>% select(date, stato, cases) %>% 
  pivot_wider(names_from = "date", values_from = "cases") %>% select(-stato) %>% as.matrix()
mat.raw[!is.na(mat.raw)] <- log(mat.raw[!is.na(mat.raw)] + 1)
cases <- list(mat.raw = mat.raw)

mat.raw <- d %>% select(date, stato, transit_pcfb) %>% 
  pivot_wider(names_from = "date", values_from = "transit_pcfb") %>% select(-stato) %>% as.matrix()
transit <- list(mat.raw = mat.raw)

mat.raw <- d %>% select(date, stato, retail_and_recreation_pcfb) %>% 
  pivot_wider(names_from = "date", values_from = "retail_and_recreation_pcfb") %>% select(-stato) %>% as.matrix()
retail_and_recreation <- list(mat.raw = mat.raw)

mat.raw <- d %>% select(date, stato, grocery_and_pharmacy_pcfb) %>% 
  pivot_wider(names_from = "date", values_from = "grocery_and_pharmacy_pcfb") %>% select(-stato) %>% as.matrix()
grocery_and_pharmacy <- list(mat.raw = mat.raw)

mat.raw <- d %>% select(date, stato, parks_pcfb) %>% 
  pivot_wider(names_from = "date", values_from = "parks_pcfb") %>% select(-stato) %>% as.matrix()
parks <- list(mat.raw = mat.raw)

mat.raw <- d %>% select(date, stato, workplaces_pcfb) %>% 
  pivot_wider(names_from = "date", values_from = "workplaces_pcfb") %>% select(-stato) %>% as.matrix()
workplaces <- list(mat.raw = mat.raw)

mat.raw <- d %>% select(date, stato, residential_pcfb) %>% 
  pivot_wider(names_from = "date", values_from = "residential_pcfb") %>% select(-stato) %>% as.matrix()
residential <- list(mat.raw = mat.raw)

m <- transit$mat.raw; loglambdaseq = seq(-6, -2, len = 20); usema = F; i <- j <- 1
# smooth and evaluate
get_fd <- function(m, loglambdaseq = NULL, usema = F, penderiv = 2){
  
  # imputazione na per le date iniziali (si assume assenza di rilevazioni)
  for (i in 1:nrow(m)) m[i, 1:min(which(!is.na(m[i,])))] <- 0
  
  # parametri per il lisciamento
  day <- seq(1, ncol(m), by = 7)
  M <- 6
  k <- length(day) + M - 2
  base <- create.bspline.basis(range(day), k, M, day)
  lfd.degree <- penderiv
  
  # penalita'
  loglambda <- seq(3, 5, len = 20)
  if(!is.null(loglambdaseq)) loglambda <- loglambdaseq
  metrics <- matrix(NA, length(loglambda), 3)
  colnames(metrics) <- c("lambda", "df", "gcv")
  betas <- array(NA, dim = c(nrow(m), k, length(loglambda)))
  
  # ciclo di smoothing
  for (j in seq_along(loglambda)){
    
    cat("Using lambda ", j, "/", length(loglambda), "...\n")
    fdPar <- fdPar(base, lfd.degree, lambda = exp(loglambda[j]))
    
    metrics2 <- matrix(NA, nrow(m), 2)
    
    # una serie per volta
    for(i in 1:nrow(m)){
      
      #if(i %% 5 == 0) cat("\tSmoothing series", i, "/", nrow(m), "...\n")
      s <- as.numeric(na.omit(m[i,]))
      ids <- which(!is.na(m[i,]))
      # le successive due righe calcolano la media mobile della serie
      ids <- ids[4:(length(ids)-3)]
      s <- ifelse(rep(usema,length(ids)), as.numeric(na.omit(forecast::ma(s, 6))), s[4:(length(s)-3)])
      sb <- smooth.basis(ids, s, fdPar)
      metrics2[i,] <- c(sb$df, sb$gcv)
      betas[i,,j] <- sb$fd$coefs
      
    }
    
    # calcolo metriche
    metrics[j,] <- c(loglambda[j], apply(metrics2, 2, mean))
    
  }
  
  # visualizzazione metriche di errore
  plot(metrics[,1], metrics[,3], xlab = expression(log(lambda)), ylab = "GCV", type = "b")
  abline(v = metrics[which.min(metrics[,3]) ,1], col = 2, lty = 2)
  #plot(metrics[,1], metrics[,2], xlab = expression(log(lambda)), ylab = "df", type = "b")
  
  # output
  betabest <- betas[,,which.min(metrics[,3])]
  myfd <- fd(coef = t(betabest), basisobj = base)
  out <- list(fd = myfd, metrics = metrics)
  
  out
}


# deaths
system.time(out <- get_fd(deaths$mat.raw, loglambdaseq = seq(4, 4.7, len = 40), usema = T, penderiv = 4))
deaths$fd <- out$fd
deaths$metrics <- out$metrics
plot(deaths$fd)

# cases
out <- get_fd(cases$mat.raw, loglambdaseq = seq(4, 4.7, len = 40), usema = T, penderiv = 4)
cases$fd <- out$fd
cases$metrics <- out$metrics
plot(cases$fd)

# transit
out <- get_fd(transit$mat.raw, loglambdaseq = seq(4, 6, len = 40))
transit$fd <- out$fd
transit$metrics <- out$metrics
plot(transit$fd)

# retail_and_recreation
out <- get_fd(retail_and_recreation$mat.raw, loglambdaseq = seq(5, 6, len = 40))
retail_and_recreation$fd <- out$fd
retail_and_recreation$metrics <- out$metrics
plot(retail_and_recreation$fd)

# grocery_and_pharmacy
out <- get_fd(grocery_and_pharmacy$mat.raw, loglambdaseq = seq(6, 10, len = 40))
grocery_and_pharmacy$fd <- out$fd
grocery_and_pharmacy$metrics <- out$metrics
plot(grocery_and_pharmacy$fd)

# parks
out <- get_fd(parks$mat.raw, loglambdaseq = seq(2, 5, len = 40))
parks$fd <- out$fd
parks$metrics <- out$metrics
plot(parks$fd)

# workplaces
out <- get_fd(workplaces$mat.raw, loglambdaseq = seq(5, 6.5, len = 40))
workplaces$fd <- out$fd
workplaces$metrics <- out$metrics
plot(workplaces$fd)

# residential
out <- get_fd(residential$mat.raw, loglambdaseq = seq(5, 7, len = 40))
residential$fd <- out$fd
residential$metrics <- out$metrics
plot(residential$fd)


# salvataggio ------------------------------------------------------------------

save(d, deaths, cases, transit, retail_and_recreation, 
     grocery_and_pharmacy, parks, workplaces, residential, file = "data/smoothed.RData")
