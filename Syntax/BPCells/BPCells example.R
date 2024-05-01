library(dplyr)
library(ggplot2)
library(umap)
library(reticulate)

use_python("/N/soft/rhel8/python/gnu/3.10.5/bin/python")
py_config()
py_available()

source("Algorithms/DR metrics.R")

load("Data/BPCells clean data.Rda")

###
Z = unname(pca$x)
Z_dist = as.matrix(dist(Z))

n = dim(Z)[1]
p = dim(Z)[2]
r = 4

Y = unname(pca$x[,1:r])
Y_dist = as.matrix(dist(Y))

###

neighbors = seq(from=5, to=100, by=5)
b = 20
k = 50

trusts = vector(length = length(neighbors)*b)
trusts_noise = vector(length = length(neighbors)*b)

best_X = NA
best_n_neighbors = NA
best_trust = -Inf

best_X_noise = NA
best_n_neighbors_noise = NA
best_trust_noise = -Inf

print("Starting loop!")
count = 1
for (i in 1:length(neighbors)) {
  for (j in 1:b) {
    trans_umap = umap(Z, method = "umap-learn", n_neighbors = neighbors[i], n_components = 2)
    X_umap = trans_umap$layout
    
    X_dist = as.matrix(dist(X_umap))
    
    sample_indices = sample(1:n, 1000)
    
    trust = trustworthiness_full_approx2(Y_dist, X_dist, k, sample_indices)
    trust_noise = trustworthiness_full_approx2(Z_dist, X_dist, k, sample_indices)
    
    if (trust > best_trust) {
      best_X = X_umap
      best_trust = trust
      best_n_neighbors = neighbors[i]
    }
    
    if (trust_noise > best_trust_noise) {
      best_X_noise = X_umap
      best_trust_noise = trust_noise
      best_n_neighbors_noise = neighbors[i]
    }
    
    trusts[count] = trust
    trusts_noise[count] = trust_noise
    
    print(paste0(count, "/", length(neighbors)*b, " loops complete!"))
    
    count = count + 1
  }
}

df = data.frame(n_neighbors = rep(neighbors, each = b),
                trust = trusts,
                trust_noise = trusts_noise)

save(Y, df, best_X, best_n_neighbors, best_trust, best_X_noise, best_n_neighbors_noise, best_trust_noise, file = "~/DR noise/Output/BPCells example.Rda")