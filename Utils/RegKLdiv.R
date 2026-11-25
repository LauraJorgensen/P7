# Function for calculating the regularized generalised KL-divergence
gkl.dev.reg <- function(V, W, H, lambda, k, delta = 1){
  div <- gkl.dev(as.vector(V), as.vector(W %*% H))
  X <- t(W) %*% W + delta * diag(nrow = k)
  return(div + lambda * logdet(X))
}