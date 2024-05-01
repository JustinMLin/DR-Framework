library(dplyr)
library(datarix)
library(irlba)

load("BPCells matrix.Rda")

data = unname(t(as.matrix(mat_norm)))

total_var = sum(diag(t(data) %*% data))/nrow(data)
gc()

pca = prcomp_irlba(data, n=500)

sum(pca$sdev^2)/total_var
# contains 86% of variance

save(data, pca, file="BPCells clean data.Rda")