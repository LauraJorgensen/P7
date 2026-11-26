# Function for tuning lambda tilde
tuning_lambda = function(V, k, seed = 42, grid = c(0.0000001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2), rel_tol = 0.05, tol = 1e-5, maxiter = 1000, delta = 1){
  grid <- sort(grid)
  gkl <- rep(0, length(grid))
  fit <- vector("list", length(grid))
  
  for(j in seq_along(grid)){
    fit[[j]] <- MINVOLNMF(V, k, seed = seed, tol = tol, maxiter = maxiter, lambda_tilde = grid[j], delta = delta)
    gkl[j] <- fit[[j]]$gkl
  }
  
  ref <- gkl[1]
  ok <- is.finite(gkl) & gkl <= (1 + rel_tol) * ref
  best <- max(which(ok))
  
  list(lambda_tilde = grid[best],
       fit = fit[[best]],
       table = data.frame(lambda_tilde = grid, gkl = gkl, accepted = ok))
}
