# Function to calculate BIC
gkl.bic <- function(k, V, W, H){
  V <- as.matrix(V)
  m <- dim(V)[1]
  n <- dim(V)[2]
  p <- m * k + k * n - k
  t <- m * n
  gkl <- gkl.dev(as.vector(V), as.vector(W %*% H))
  model_complexity <- p * log(t)
  return(model_complexity + 2*gkl)
}