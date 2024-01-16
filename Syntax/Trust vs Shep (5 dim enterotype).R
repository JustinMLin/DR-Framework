library(ggplot2)
library(MASS)
library(dplyr)
library(Rtsne)
library(viridis)
library(gridExtra)
library(tsne)
library(phyloseq)

# source("../../../Algorithms/DR metrics.R")
source("Algorithms/DR metrics.R")

data(enterotype)

data = as.data.frame(t(enterotype@otu_table))

data = data %>%
  select(where(function(x) mean(x < 0.1e-05) < 0.6)) %>%
  mutate_all(function(x) log(1+x))

###

pca = prcomp(data, center = TRUE, scale. = TRUE, rank. = 66)
Z = pca$x

###

n = dim(Z)[1]
p = dim(Z)[2]
r = 5

Y = Z[,(1:r)]

###

perplexity = seq(from = 5, to = 90, by = 5)
b = 20
k = 10

trusts = sheps = vector(length = length(perplexity)*b)
trusts_noise = sheps_noise = vector(length = length(perplexity)*b)

print("Starting loop!")

count = 1
for (i in 1:length(perplexity)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity[i], theta = 0.0)$Y
    
    trusts[count] = trustworthiness_full(Y, X_tsne, k)
    sheps[count] = dist_cor_full(Y, X_tsne)
    
    trusts_noise[count] = trustworthiness_full(Z, X_tsne, k)
    sheps_noise[count] = dist_cor_full(Z, X_tsne)
    
    print(paste0(count, "/", length(perplexity)*b, " loops complete!"))
    
    count = count + 1
  }
}

###

df = data.frame(perplexity = rep(perplexity, each = b),
                trust = trusts,
                shep = sheps,
                trust_noise = trusts_noise,
                shep_noise = sheps_noise)

p1 = ggplot(df, aes(x = trust, y = shep, col = perplexity)) + 
  geom_point(size = 2) + 
  scale_color_viridis() + 
  labs(x = "Trustworthiness", y = "Shepard Goodness", title = "Replicating Signal")


p2 = ggplot(df, aes(x = trust_noise, y = shep_noise, col = perplexity)) + 
  geom_point(size = 2) + 
  scale_color_viridis() + 
  labs(x = "Trustworthiness", y = "Shepard Goodness", title = "Replicating Signal + Noise")

# save(Y, df, p1, p2, file = "trust vs shep (5 dim enterotype).Rda")
save(Y, df, p1, p2, file = "~/DR noise/Output/trust vs shep (5 dim enterotype).Rda")
