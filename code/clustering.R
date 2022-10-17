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

library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
world$name[world$name == "Czech Rep."] <- "Czechia"
class(world)

eu <- world %>% 
  filter(name %in% dd$stato) %>% 
  select(name, geometry)

ggplot(eu, aes(fill = "blue")) + 
  geom_sf() + 
  scale_x_continuous(limits = c(-10, 30)) +
  scale_y_continuous(limits = c(35, 65))


# extra: clustering of cases ---------------------------------------------------

#dim(cases$fd$coefs)
#betas <- t(cases$fd$coefs)
#dimnames(betas)[[1]] <- unique(d$stato)

# basato su scores e pca
xmat <- cases$pca$scores
rownames(xmat) <- unique(d$stato)
xmat

# la variabilita'/caratt.individuali delle curve sono racchiuse nei beta, 
# che moltiplicano funz tutte uguali
# -> si puo' fare clustering con un semplice k-means sui beta
# l'unico problema di questo approccio e' non considerare la variabilita' dei 
# coeff stimani, amen ce ne faremo una ragione

set.seed(42)
tryk <- 1:26
withinss <- numeric(length(tryk))
for (i in seq_along(tryk)){
  cat("Trying with", tryk[i], "centers...\n")
  km <- kmeans(xmat, centers = tryk[i], nstart = 200, iter.max = 1000)
  withinss[i] <- km$withinss
}
plot(tryk, withinss, main = "Within SS", xlab = "k", ylab = "", type = "b")
#plot(tryk[-1], diff(withinss), main = "Within SS gain", xlab = "k", ylab = "", type = "b")
# possiamo provare 5, 6, 7

set.seed(42)
km1 <- kmeans(xmat, centers = 4, nstart = 200, iter.max = 1000)
matplot(cases$mat.sm, type = "l", col = km1$cluster, lty = 1, xlab = "")
legend("bottomright", lty = 1, col = 1:4, legend = paste0("group ",1:4), bty = "n", cex = 0.8)

# x report
pdf("plots/kmeans1.pdf", width = 7, height = 5)
matplot(cases$mat.sm, type = "l", col = km1$cluster, lty = 1, xlab = "")
legend("bottomright", lty = 1, col = 1:4, legend = paste0("group ",1:4), bty = "n", cex = 0.8)
dev.off()

for (j in 1:4) cat("Cluster ", j, ": ", levels(d$stato)[km1$cluster == j], "\n")

hc <- hclust(dist(xmat), method = "complete")
plot(hc)
#plot(hc, hang = -1)
rect.hclust(hc, k = 4)
hc1 <- cutree(hc, k = 4)

# x report
pdf("plots/hclust1.pdf", width = 7, height = 5)
plot(hc, xlab = "", ylab = "", axes = F, frame.plot = F, ann = F)
rect.hclust(hc, k = 4)
dev.off()


# visualisations ---------------------------------------------------------------

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
#ggsave("plots/kmeans_map.pdf", plot = plot_kmeans, width = 7, height = 5, units = "in")


hclust1 <- data.frame(name = names(hc1), cluster2 = hc1, row.names = NULL)
#hclust1$cluster[hclust1$cluster2 == 3] <- 4
#hclust1$cluster[hclust1$cluster2 == 4] <- 2
hclust1$cluster[hclust1$cluster2 == 2] <- 1
hclust1$cluster[hclust1$cluster2 == 1] <- 2

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
#ggsave("plots/hclust_map.pdf", plot = plot_hclust, width = 7, height = 5, units = "in")


pp <- ggpubr::ggarrange(plot_kmeans, plot_hclust, common.legend = T, legend = "bottom")
ggsave("plots/clustering_maps.pdf", plot = pp, width = 10, height = 5, units = "in")
