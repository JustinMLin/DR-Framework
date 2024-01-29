library(ggplot2)
library(MASS)
library(dplyr)
library(Rtsne)
library(viridis)
library(gridExtra)
library(EFAtools)

# source("../../../Algorithms/DR metrics.R")
source("Algorithms/DR metrics.R")

# data = read.delim("../../../../../CyTOF data/Exp1_NK_uns_CD4_inf.txt")
data = read.delim("Data/Exp1_NK_uns_CD4_inf.txt")
data = data %>%
  select(-c(1,2)) %>%
  sample_n(5000) %>%
  mutate_all(function(x) log(1+x))

# par = PARALLEL(data)
# data.frame(par$eigenvalues_PCA) %>%
#   janitor::clean_names() %>%
#   mutate(x = 1:length(real_eigenvalues)) %>%
#   ggplot(aes(x=x, y=real_eigenvalues)) +
#     geom_col() +
#     labs(y="Eigenvalue", title="Scree Plot for CyTOF Dataset") +
#     theme(axis.title.x=element_blank())

###

pca = prcomp(data, center = TRUE, scale. = TRUE, rank. = 30)
Z = pca$x

###

n = dim(Z)[1]
p = dim(Z)[2]
r = 5

Y = Z[,(1:r)]

###

perplexity = seq(from=10, to=300, by=10)
b = 20
k = 50

trusts = sheps = vector(length = length(perplexity)*b)
trusts_noise = sheps_noise = vector(length = length(perplexity)*b)

print("Starting loop!")

count = 1
for (i in 1:length(perplexity)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity[i])$Y
    
    sample_indices = sample(1:n, 200)
    
    trusts[count] = trustworthiness_full_approx(Y, X_tsne, k, sample_indices)
    sheps[count] = dist_cor_full_approx(Y, X_tsne, sample_indices)
    
    trusts_noise[count] = trustworthiness_full_approx(Z, X_tsne, k, sample_indices)
    sheps_noise[count] = dist_cor_full_approx(Z, X_tsne, sample_indices)
    
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

# save(Y, df, p1, p2, file = trust vs shep (5 dim CyTOF).Rda")
save(Y, df, p1, p2, file = "~/DR noise/Output/trust vs shep (5 dim CyTOF).Rda")
