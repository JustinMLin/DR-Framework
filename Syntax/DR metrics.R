library(dbscan)

################################################################################
## Helper Functions

get_nns = function(pt, data, k) {
  kNN(data, k, matrix(pt, nrow = 1))$id
}

get_nns2 = function(id, data_dist, k) {
  order(data_dist[id,])[1:k]
}

nn_rank = function(pt, j, data) {
  n = length(data[,1])
  rank = match(j, kNN(data, n-1, matrix(pt, nrow = 1))$id)
  
  if (is.na(rank)) n else rank
}

nn_rank2 = function(id, j, data_dist) {
  n = length(data_dist[,1])
  rank = match(j, order(data_dist[id,]))
  
  if (is.na(rank)) n else rank
}

################################################################################
## DR

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

trustworthiness_full2 = function(Z_dist, X_dist, k) {
  num_pts = length(Z_dist[,1])
  
  total = 0
  for (i in 1:num_pts) {
    high_dim_neighbors = get_nns2(i, Z_dist, k+1)[-1]
    low_dim_neighbors = get_nns2(i, X_dist, k+1)[-1]
    
    U = setdiff(low_dim_neighbors, high_dim_neighbors)
    
    if (length(U) != 0) {
      total = total + sum(sapply(U, function(j) nn_rank2(i, j, Z_dist) - 1 - k))
    }
  }
  
  1 - 2/(num_pts *k*(2*num_pts  - 3*k - 1))*total
}

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

trustworthiness_full_approx2 = function(Z_dist, X_dist, k, indices) {
  num_pts = length(Z_dist[,1])
  b = length(indices)
  
  total = 0
  for (i in indices) {
    high_dim_neighbors = get_nns2(i, Z_dist, k+1)[-1]
    low_dim_neighbors = get_nns2(i, X_dist, k+1)[-1]
    
    U = setdiff(low_dim_neighbors, high_dim_neighbors)
    
    if (length(U) != 0) {
      total = total + sum(sapply(U, function(j) nn_rank2(i, j, Z_dist) - 1 - k))
    }
  }
  
  1 - 2/(b*k*(2*num_pts  - 3*k - 1))*total
}
  
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

dist_cor_full = function(Z, X) {
  cor(dist(Z), dist(X), method="spearman")
}

dist_cor_full_approx = function(Z, X, indices) {
  cor(dist(Z[indices,]), dist(X[indices,]), method="spearman")
}