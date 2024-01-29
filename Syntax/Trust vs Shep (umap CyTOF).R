library(ggplot2)
library(MASS)
library(dplyr)
library(umap)
library(viridis)
library(gridExtra)
library(reticulate)

use_python("/N/soft/rhel8/python/gnu/3.10.5/bin/python")
py_config()
py_available()

# source("../../../Algorithms/DR metrics.R")
source("Algorithms/DR metrics.R")

# data = read.delim("../../../../../CyTOF data/Exp1_NK_uns_CD4_inf.txt")
data = read.delim("Data/Exp1_NK_uns_CD4_inf.txt")
data = data %>%
  select(-c(1,2)) %>%
  sample_n(5000) %>%
  mutate_all(function(x) log(1+x))

###

pca = prcomp(data, center = TRUE, scale. = TRUE, rank. = 30)
Z = pca$x

###

n = dim(Z)[1]
p = dim(Z)[2]
r = 5

Y = Z[,c(1,2)]

###

neighbors = seq(from=10, to=300, by=10)
b = 20
k = 50

trusts = sheps = vector(length = length(neighbors)*b)
trusts_noise = sheps_noise = vector(length = length(neighbors)*b)

print("Starting loop!")

count = 1
for (i in 1:length(neighbors)) {
  for (j in 1:b) {
    trans_umap = umap(Z, method = "umap-learn", n_neighbors = neighbors[i], n_components = 2)
    X_umap = trans_umap$layout
    
    sample_indices = sample(1:n, 200)
    
    trusts[count] = trustworthiness_full_approx(Y, X_umap, k, sample_indices)
    sheps[count] = dist_cor_full_approx(Y, X_umap, sample_indices)
    
    trusts_noise[count] = trustworthiness_full_approx(Z, X_umap, k, sample_indices)
    sheps_noise[count] = dist_cor_full_approx(Z, X_umap, sample_indices)
    
    print(paste0(count, "/", length(neighbors)*b, " loops complete!"))
    
    count = count + 1
  }
}

###

df = data.frame(n_neighbors = rep(neighbors, each = b),
                trust = trusts,
                shep = sheps,
                trust_noise = trusts_noise,
                shep_noise = sheps_noise)

p1 = ggplot(df, aes(x = trust, y = shep, col = n_neighbors)) + 
  geom_point(size = 2) + 
  scale_color_viridis() + 
  labs(x = "Trustworthiness", y = "Shepard Goodness", title = "Replicating Signal")


p2 = ggplot(df, aes(x = trust_noise, y = shep_noise, col = n_neighbors)) + 
  geom_point(size = 2) + 
  scale_color_viridis() + 
  labs(x = "Trustworthiness", y = "Shepard Goodness", title = "Replicating Signal + Noise")

# save(Y, df, p1, p2, file = trust vs shep (umap CyTOF).Rda")
save(Y, df, p1, p2, file = "~/DR noise/Output/trust vs shep (umap CyTOF).Rda")
