gen_eigm <- function(d) {
  # d is a vector of eigenvalues
  # returns a matrix with those eigenvalues
  sigma = diag(d) # diagonal matrix of eigenvalues
  nm = dim(sigma)
  V = matrix(rnorm(prod(nm)), ncol=nm[2])   # random matrix
  A = solve(V)%*%sigma%*%V        # random A with given eigenvalues
  return(A)
}


gen_eigm(c(1,1))
evals <- c(2,0.3)
A <- gen_eigm(evals)
SVDA <- svd(A)
# t(SVDA$u) %*% diag(SVDA$d) %*% SVDA$v
EIGA <- eigen(A)

SVDA$d;EIGA$values

norm(as.matrix(SVDA$d));norm(as.matrix(EIGA$values))

n <- 50
evals <- c(2,0.3)
SVDA <- matrix(NA, ncol=2, nrow=n)
for (i in 1:n) {
  A <- gen_eigm(evals)
  SVDA[i,] <- svd(A)$d
}
plot(rbind(evals,SVDA))
plot(SVDA[,1]*SVDA[,2])
     