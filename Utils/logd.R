# Numerically stable log-determinant
logdet <- function(A) as.numeric(determinant(A, logarithm = TRUE)$modulus)

