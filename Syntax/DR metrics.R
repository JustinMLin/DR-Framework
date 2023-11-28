library(dbscan)

################################################################################
## Helper Functions

# returns the k-nearest neighbors of pt in data
get_nns = function(pt, data, k) {
  kNN(data, k, matrix(pt, nrow = 1))$id
}

# returns the rank of point j among pt's neighbors in data
nn_rank = function(pt, j, data) {
  n = length(data[,1])
  rank = match(j, kNN(data, n-1, matrix(pt, nrow = 1))$id)
  
  if (is.na(rank)) n else rank
}

################################################################################
## DR

# computes trustworthiness using every point
trustworthiness_full = function(Z, X, k) {
  num_pts = length(Z[,1])
  
  total = 0
  for (i in 1:num_pts) {
    high_dim_neighbors = get_nns(Z[i,], Z, k+1)[-1]
    low_dim_neighbors = get_nns(X[i,], X, k+1)[-1]
    
    U = setdiff(low_dim_neighbors, high_dim_neighbors)
    
    if (length(U) != 0) {
      total = total + sum(sapply(U, function(j) nn_rank(Z[i,], j, Z) - 1 - k))
    }
  }
  
  1 - 2/(num_pts *k*(2*num_pts  - 3*k - 1))*total
}

# approximates trustworthiness by downsampling to points with the given indices
trustworthiness_full_approx = function(Z, X, k, indices) {
  num_pts = length(Z[,1])
  b = length(indices)
  
  total = 0
  for (i in indices) {
    high_dim_neighbors = get_nns(Z[i,], Z, k+1)[-1]
    low_dim_neighbors = get_nns(X[i,], X, k+1)[-1]
    
    U = setdiff(low_dim_neighbors, high_dim_neighbors)
    
    if (length(U) != 0) {
      total = total + sum(sapply(U, function(j) nn_rank(Z[i,], j, Z) - 1 - k))
    }
  }
  
  1 - 2/(b*k*(2*num_pts  - 3*k - 1))*total
}

# computes continuity using every point
continuity_full = function(Z, X, k) {
  n = length(Z[,1])
  
  total = 0
  for (i in 1:n) {
    high_dim_neighbors = get_nns(Z[i,], Z, k+1)[-1]
    low_dim_neighbors = get_nns(X[i,], X, k+1)[-1]
    
    V = setdiff(high_dim_neighbors, low_dim_neighbors)
    
    if (length(V != 0)) {
      total = total + sum(sapply(V, function(j) nn_rank(X[i,], j, X) - 1 - k))
    }
  }
  
  1 - 2/(n*k*(2*n - 3*k - 1))*total
}

# computes local stress
local_stress_full = function(Z, X, k) {
  n = length(Z[,1])
  
  total_stress = 0
  for (i in 1:n) {
    neighbors = get_nns(Z[i,], Z, k+1)[-1]
    
    z_mat = matrix(Z[i,], nrow = k, ncol = length(Z[i,]), byrow = TRUE)
    z_dists = sqrt(rowSums((Z[neighbors,] - z_mat)^2))
    
    x_mat = matrix(X[i,], nrow = k, ncol = length(X[i,]), byrow = TRUE)
    x_dists = sqrt(rowSums((X[neighbors,] - x_mat)^2))
    
    total_stress = total_stress + sum((z_dists - x_dists)^2) / sum((z_dists)^2)
  }
  
  total_stress/n
}

# computes Shepard goodness using every point
dist_cor_full = function(Z, X) {
  cor(dist(Z), dist(X), method="spearman")
}

# approximates Shepard goodness by first downsampling
dist_cor_full_approx = function(Z, X, indices) {
  cor(dist(Z[indices,]), dist(X[indices,]), method="spearman")
}