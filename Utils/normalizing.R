# Function for normalizing W and H
normalize <- function(W, H){
  s <- colSums(W)
  list(W = t(t(W) / s), H = s * H)
}