gen_eigm <- function(d) {
  # d is a vector of eigenvalues
  # returns a matrix with those eigenvalues
  sigma = diag(d) # diagonal matrix of eigenvalues
  nm = dim(sigma)
  V = matrix(runif(prod(nm), -5, 5), ncol=nm[2])   # random matrix
  A = solve(V)%*%sigma%*%V        # random A with given eigenvalues
  return(A)
}

evals <- c(2,-2)
A <- gen_eigm(evals)
SVDA <- svd(A)
# SVDA$u %*% diag(SVDA$d) %*% t(SVDA$v)
EIGA <- eigen(A)

SVDA$d;EIGA$values

prod(SVDA$d);prod(EIGA$values)

n <- 50
evals <- c(2,-0.3)
SVDA <- matrix(NA, ncol=length(evals), nrow=n)
for (i in 1:n) {
  A <- gen_eigm(evals)
  SVDA[i,] <- svd(A)$d
}
plot(rbind(evals,SVDA))
plot(SVDA[,1]*SVDA[,2])
     
n <- 50
evals <- c(2,0.5,-0.3)
SVDA <- matrix(NA, ncol=length(evals), nrow=n)
for (i in 1:n) {
  A <- gen_eigm(evals)
  SVDA[i,] <- svd(A)$d
}
#plot(rbind(evals,SVDA))
plot(SVDA[,1]*SVDA[,2]*SVDA[,3])
