F <- function (x) {3 - 2*x^2} # an quadratic function with two roots
curve(F, -1.5, 1.5)
f <- function (x) {-4*x} # the derivative function of F
x <- 0.1 # note we cannot start at x==0
while (abs(F(x))> 1e-8) {
  x <- x - F(x)/f(x)  # Newton's Method
  print(c(x, F(x)))
}

x <- -0.1 # note we cannot start at x==0
while (abs(F(x))> 1e-8) {
  x <- x - F(x)/f(x)  # Newton's Method
  print(x)
}

Fx <- function (x) {
  c(x[1]*x[2],
    x[1]+x[2])
  }
J <- function (x) {
  x1 <- c(x[2],1)
  x2 <- c(x[1],1)
  cbind(x1,x2)
}

x0 <- c(2,1) 
solve(J(x0),-Fx(x0))+x0
x <- solve(J(x0),-Fx(x0))+x0
x0 <- x
while (norm(as.matrix(Fx(x0)), "F") > 1e-8) {
  x <- solve(J(x0),-Fx(x0))+x0
  x0 <- x
}
x

Fx <- function (x) {
  c(x[1]^2+x[2]^2+x[3]^2-3,
    x[1]^2+x[2]^2-x[3]-1,
    x[1]+x[2]+x[3]-3)
}
J <- function (x) {
  x1 <- c(2*x[1], 2*x[1], 1)
  x2 <- c(2*x[2], 2*x[2], 1)
  x3 <- c(2*x[3], -1, 1)
  cbind(x1,x2,x3)
}
x0 <- c(1,0,1) 
solve(J(x0),-Fx(x0))+x0

while (norm(as.matrix(Fx(x0)), "F") > 1e-8) {
  x <- solve(J(x0),-Fx(x0)) + x0
  x0 <- x
  print(x)
}
x