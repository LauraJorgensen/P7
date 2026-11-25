# Function for calculating the generalised KL-divergence
gkl.dev <- function(x, y){
  r <- y 
  p <- which(x > 0)   # 0 log(0) = 0
  r[p] <- (x * (log(x) - log(y)) - x + y)[p]   
  return(sum(r))
}

# Numerically stable log-determinant
logdet <- function(A) as.numeric(determinant(A, logarithm = TRUE)$modulus)

# Function for calculating the regularized generalised KL-divergence
gkl.dev.reg <- function(V, W, H, lambda, k, delta = 1){
  div <- gkl.dev(as.vector(V), as.vector(W %*% H))
  X <- t(W) %*% W + delta * diag(nrow = k)
  return(div + lambda * logdet(X))
}

# Function for normalizing W and H
normalize <- function(W, H){
  s <- colSums(W)
  list(W = t(t(W) / s), H = s * H)
}

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

# Testing different choices of k and finding BIC
k.choice <- function(trials = 3:10, V){
  Wlist <- list()
  Hlist <- list()
  biclist <- rep(NA_real_, length(trials))
  for (j in seq_along(trials)){
    k <- trials[j]
    fit <- KLNMF(V, k = k, seed = 1:10, tol = 1e-8, maxiter = 1000)
    Wlist[[j]] <- fit$W
    Hlist[[j]] <- fit$H
    biclist[j] <-gkl.bic(k, V, fit$W, fit$H)
  }
  best <- which.min(biclist)
  W <- Wlist[[best]]
  H <- Hlist[[best]]
  
  Output <- list()
  Output$W <- W
  Output$H <- H
  Output$bic <- biclist[best]
  Output$k <- trials[best]
  Output$bicall <- setNames(biclist, trials)
  return(Output)
}

# KL-NMF algorithm as described in Section 3.2 and Section 3.4
KLNMF = function(V, k, seed = 42, arrange = TRUE, tol = 1e-5, maxiter = 1000){
  
  m <- dim(V)[1]  
  n <- dim(V)[2]  
  eps <- 1e-10   # to avoid underflow
  
  div <- rep(0, length(seed)) 
  Wlist <- list()
  Hlist <- list()
  reslist <- list() # for each initialization, we store both the number of iterations and the KL-divergence for each iteration. 
  
  for(i in 1:length(seed)){
    set.seed(seed[i])
    
    # Initialize W and H
    W <- matrix(stats::runif(m*k), nrow = m, ncol = k)  
    H <- matrix(stats::runif(k*n), nrow = k, ncol = n)  
    
    ones <- matrix(1, nrow = m, ncol = n) 
    
    obj <- rep(0, maxiter) # to store 
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
    nm = normalize(W, H)
    W <- nm$W
    H <- nm$H
    
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


# Min-Vol KL-NMF algorithm as described in Section 3.3 and Section 3.4
MINVOLNMF = function(V, k, seed = 42, arrange = TRUE, tol = 1e-5, maxiter = 1000, lambda_tilde = 0.01){
  m <- dim(V)[1]
  n <- dim(V)[2]
  eps <- 1e-10   # to avoid underflow
  
  lambdas <- rep(0, length(seed))
  div <- rep(0, length(seed)) 
  objval <- rep(0, length(seed)) 
  Wlist <- list()
  Hlist <- list()
  reslist <- list() 
  
  for(i in 1:length(seed)){
    set.seed(seed[i])
    gamma <- 1
    
    # Initialize W and H
    W <- matrix(stats::runif(m*k), nrow = m, ncol = k)  
    H <- matrix(stats::runif(k*n), nrow = k, ncol = n)  
    
    WH0 <- W %*% H
    WH0[WH0 <= 0] <- eps
    lambda <- lambda_tilde * gkl.dev(as.vector(V), as.vector(WH0)) / max(abs(logdet(t(W) %*% W + diag(nrow = k))), eps)
    lambdas[i] <- lambda
    
    nm <- normalize(W, H)
    W <- nm$W
    H <- nm$H
    
    ones <- matrix(1, nrow = m, ncol = n) 
    
    obj <- rep(0, maxiter) 
    niter <- maxiter
    
    for(t in 1:maxiter){
      # Update of H
      WH <- W %*% H
      WH[WH <= 0] <- eps
      H <- H * (t(W) %*% (V / WH)) / (t(W) %*% ones)
      H[H <= 0] <- eps
      
      Y <- solve(t(W) %*% W + diag(nrow = k))
      
      Y_plus <- pmax(Y,0)
      Y_minus <- pmax(-Y,0)
      
      # Update of W 
      WH <- W %*% H
      WH[WH <= 0] <- eps
      
      Wup <- W * (((ones %*% t(H) - 4 * lambda * (W %*% Y_minus))^2 + (8 * lambda * W %*% (Y_plus + Y_minus)) * ((V / (WH)) %*% t(H)))^(1/2) - (ones %*% t(H) - 4 * lambda * W %*% Y_minus)) / (4 * lambda * W %*% (Y_plus + Y_minus) + eps)
      Wup[Wup <= 0] <- eps
      
      # Normalize
      nm <- normalize(Wup, H)
      Wnom <- nm$W
      Hnom <- nm$H
      
      # Line-search
      while (gkl.dev.reg(V, Wnom, Hnom, lambda, k) > gkl.dev.reg(V, W, H, lambda, k) &&
             gamma > eps) {
        gamma <- gamma * 0.8
        Wmix <- (1 - gamma) * W + gamma * Wup
        nm <- normalize(Wmix, H)
        Wnom <- nm$W
        Hnom <- nm$H
      }
      W <- Wnom
      H <- Hnom
      
      gamma <- min(1, gamma * 1.2)
      
      # Stopping Criterion
      obj[t] <- gkl.dev.reg(V, W, H, lambda, k)
      if(t > 1 && abs(obj[t-1] - obj[t]) / (abs(obj[t-1]) + eps) < tol){
        niter <- t; 
        break}
    }
    
    # Normalization of W
    H <- colSums(W) * H
    W <- t(t(W) * 1/colSums(W))
    
    Wlist[[i]] <- W
    Hlist[[i]] <- H
    div[i] <- gkl.dev(as.vector(V), as.vector(W %*% H))
    objval[i] <- gkl.dev.reg(V, W, H, lambda, k)
    reslist[[i]] <- list(iter = niter, objfn = obj[1:niter])
  }
  
  best <- which.min(objval)
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
  Output$obj <- objval[best] 
  Output$iter <- reslist[[best]]$iter
  Output$objfn <- reslist[[best]]$objfn
  Output$lambda <- lambdas[best]
  Output$lambda_tilde <- lambda_tilde
  
  return(Output)
}

# Function for tuning lambda tilde
tuning_lambda = function(V, k, seed = 42, grid = c(0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.2), rel_tol = 0.05, tol = 1e-5, maxiter = 1000){
  grid <- sort(grid)
  gkl <- rep(0, length(grid))
  fit <- vector("list", length(grid))
  
  for(j in seq_along(grid)){
    fit[[j]] <- MINVOLNMF(V, k, seed = seed, tol = tol, maxiter = maxiter, lambda_tilde = grid[j])
    gkl[j] <- fit[[j]]$gkl
  }
  
  ref <- gkl[1]
  ok <- is.finite(gkl) & gkl <= (1 + rel_tol) * ref
  best <- max(which(ok))
  
  list(lambda_tilde = grid[best],
       fit = fit[[best]],
       table = data.frame(lambda_tilde = grid, gkl = gkl, accepted = ok))
}
