signature.match <- function(X,W, tol = 0.7){
  matched <- vector("list", ncol(W))
  best.cosine <- vector("list", ncol(W))
  all.best <- vector("list", ncol(W))
  
  for (i in seq_len(ncol(W))){
    cosine.list <- rep(NA_real_, ncol(X))
    over.tol <- list()
    w <- W[,i]
    
    for (j in seq_len(ncol(X))){
      x <- X[,j]
      c <- cosine(w,x)
      cosine.list[[j]] <- c
      
      if (!is.na(c) && c > tol) {
        over.tol[[length(over.tol)+1]] <- data.frame(cosmic_index = j, 
                                                     cosmic_name = colnames(X)[j],
                                                     cosine_sim = c)
      }
    }
    best <- which.max(cosine.list)
    
    best.cosine [[i]] <- cosine.list[[best]]
    matched[[i]] <- colnames(X)[best]
    
    all.best[[i]] <- if (length(over.tol) > 0) {
      do.call(rbind, over.tol)
    } else {
      NULL
    }
  }
  return(list(
    best_match = matched,
    best_cosine = best.cosine,
    above_tol = all.best
  ))
}
