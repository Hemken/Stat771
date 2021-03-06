newton <- function (Fx, Jx, x0, eps=1e-8, max.iter=25) {
  # Fx = function
  # Jx = Jacobian of Fx
  # x0 = initial guess
  # eps = stopping threshold
  # max.iter = maximum steps
  i <- 1
  while (norm(as.matrix(Fx(x0)), "F") > eps && i <= max.iter) {
    print(solve(Jx(x0)))
    x <- solve(Jx(x0),-Fx(x0))+x0
    print(x)
    x0 <- x
    i <- i+1
  }
  return(list(solution=x0, iter=i))
}

x0 <- c(0,0)
newton(Fx, J, x0)
