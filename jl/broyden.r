broyden <- function (Fx, Jx, x0, eps=1e-8, max.iter=25) {
  D0 <- J(x0)  # initial Jacobian
  Dinvk <- solve(D0) # initial Jacobian inverse
  d0 <-  - Dinvk %*% Fx(x0) # first step
  x <- x0 + d0 # current solution

  i <- 1
  while (norm(as.matrix(Fx(x)), "F") > eps && i <= max.iter) {
    uk <- Dinvk %*% Fx(x) # next step
    ck <- (t(d0) %*% (d0 + uk ))[1,1] # normalizing constant
    Dinv <- Dinvk - (t(d0 %*% t(uk))%*%Dinvk)/ck
    d0 <- - Dinv %*% Fx(x)
    alpha <- 1
   alpha <- backtrack(Fx, J(x), as.vector(d0), x)
  #  print(paste("alpha = ", alpha))
    x <- x + alpha*d0
    i <- i + 1
  #  print(i)
  #  print(norm(as.matrix(Fx(x)), "F"))
    
  }
  return(list(solution=x, iter=i))
}

x0 <- c(0,0)
broyden(Fx, J, x0, max.iter=100)

x0 <- c(0,1)
broyden(Fx, J, x0, max.iter=100)

x0 <- c(.5,1.5)
broyden(Fx, J, x0, max.iter=1000)
