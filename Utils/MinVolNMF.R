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