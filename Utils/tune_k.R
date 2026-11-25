# Testing different choices of k and finding BIC
k.choice <- function(trials = 1:10, V, seed = 42, tol = 1e-8, maxiter = 1000){
  Wlist <- list()
  Hlist <- list()
  biclist <- rep(NA_real_, length(trials))
  for (j in seq_along(trials)){
    k <- trials[j]
    fit <- KLNMF(V, k = k, seed = seed, tol = tol, maxiter = maxiter)
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