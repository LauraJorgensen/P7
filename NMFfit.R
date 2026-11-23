# Function for calculating the generalised KL-divergence
gkl.dev <- function(x, y){
  if(length(x) != length(y)){
    stop("The input has different lengths.")
  }
  if(any(x < 0 | y < 0)){
    stop("The input is negative.")
  }
  r <- y 
  p <- which(x > 0)   # 0 log(0) = 0
  r[p] <- (x * (log(x) - log(y)) - x + y)[p]   
  return(sum(r))
}

# KL-NMF algorithm as described in Section 2.3
KLNMF = function(V, k, seed = 42, arrange = TRUE, tol = 1e-5, maxiter = 1000){
  if(any(V < 0) || k < 0){
    stop("The input is negative.")
  }
  
  m <- dim(V)[1]  
  n <- dim(V)[2]  
  eps <- 1e-10   # to avoid underflow
  
  div <- rep(0, length(seed)) 
  Wlist <- list()
  Hlist <- list()
  reslist <- list() 
  
  for(i in 1:length(seed)){
    set.seed(seed[i])
    
    # Initialize W and H
    W <- matrix(stats::runif(m*k), nrow = m, ncol = k)  
    H <- matrix(stats::runif(k*n), nrow = k, ncol = n)  
    
    ones <- matrix(1, nrow = m, ncol = n) 
    
    obj <- rep(NA_real_, maxiter)
    niter <- maxiter
    
    for(t in 1:maxiter){
      # Update of H
      WH <- W %*% H
      WH[WH <= 0] <- eps
      H <- H * (t(W) %*% (V / WH)) / (t(W) %*% ones)
      H[H <= 0] <- eps
      
      # Update of W 
      WH <- W %*% H
      WH[WH <= 0] <- eps
      W <- W * ((V / WH) %*% t(H)) / (ones %*% t(H))
      W[W <= 0] <- eps
      
      # Stopping Criterion
      obj[t] <- gkl.dev(as.vector(V), as.vector(W %*% H))
      if(t > 1 && is.finite(obj[t]) && abs(obj[t-1] - obj[t]) / (abs(obj[t-1]) + eps) < tol){
        niter <- t
        break
      }
    }
    
    # Normalization of W
    H <- colSums(W) * H
    W <- t(t(W) * 1/colSums(W))
    
    Wlist[[i]] <- W
    Hlist[[i]] <- H
    div[i] <- gkl.dev(as.vector(V), as.vector(W %*% H))
    reslist[[i]] <- list(iter = niter, objfn = obj[1:niter])
  }
  
  best <- which.min(div)
  W <- Wlist[[best]]
  H <- Hlist[[best]]
  
  # Permutation of the columns
  if(arrange == TRUE){
    idx <- order(rowSums(H), decreasing = TRUE)
    W <- W[, idx]
    H <- H[idx, ]
  }
  
  Output <- list()
  Output$W <- W
  Output$H <- H
  Output$gkl <- div[best]
  Output$iter <- reslist[[best]]$iter
  Output$objfn <- reslist[[best]]$objfn
  
  return(Output)
}


# Min-Vol KL-NMF algorithm as described in Section 2.4...













# Eksempel 1

A = matrix(c(1, 1, 0,
             2, 2, 0,
             0, 1, 2,
             0, 2, 3),
           nrow = 4, ncol = 3)

A = matrix(c(1, 1, 0,
             2, 2, 0,
             0, 1, 2,
             0, 2, 3),
           nrow = 4, ncol = 3, byrow = TRUE)

fit2 <- KLNMF(A, k = 2, seed = 1:10, tol = 1e-8, maxiter = 1000)
fit2$iter

fit2$W %*% fit2$H

# Eksempel 2

V = matrix(c(20, 3, 24, 19,  2, 15,
             9, 14, 25, 30, 15, 10,
             30, 6, 12, 10, 11,  7,
             9, 27, 5, 11, 19, 15),
           nrow = 4, ncol = 6)

fit <- KLNMF(V, k = 5, seed = 1:10, tol = 1e-8)

t0 <- system.time(fit <- KLNMF(V, k = 5, seed = 1:10, tol = 1e-8))
t0["elapsed"]

fit$W
fit$H

plot(fit$objfn,  type = "l", xlab = "Iteration", ylab = "KL-divergence")
fit$gkl
