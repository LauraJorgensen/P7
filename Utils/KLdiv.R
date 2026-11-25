# Function for calculating the generalised KL-divergence
gkl.dev <- function(x, y){
  r <- y 
  p <- which(x > 0)   # 0 log(0) = 0
  r[p] <- (x * (log(x) - log(y)) - x + y)[p]   
  return(sum(r))
}