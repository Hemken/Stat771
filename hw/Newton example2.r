g <- function(v) {
  x <- v[1]
  y <- v[2]
  z <- v[3]
  g1 <- x^2 + y^2 + z^2 - 3
  g2 <- x^2 - y^2 - z   - 1
  g3 <- x   + y   + z   - 3
  return(matrix(c(g1,g2,g3),ncol=1))
}

g0 <- g(c(1,0,1))

Jg <- function(v) {
  x <- v[1]
  y <- v[2]
  z <- v[3]
  J11 <- 2*x
  J21 <- 2*x
  J31 <- 1
  J12 <- 2*y
  J22 <- 2*y
  J32 <- 1
  J13 <- 2*z
  J23 <- -1
  J33 <- 1
  return(matrix(c(J11,J21,J31,J12,J22,J32,J13,J23,J33),ncol=3))
}

Jg(c(1,0,1))

solve(Jg(c(1,0,1)), -g(c(1,0,1)))

g1 <- c(1,0,1) - solve(Jg(c(1,0,1)), -g(c(1,0,1)))
g0 <- g1
g1 <-  g0- solve(Jg(g0), g(g0))
g0 <- g1

# g0 <- c(0.5,1,0.5)
for (i in 1:10) {
  g1 <- solve(Jg(g0), -g(g0)) + g0
  g0 <- g1
  print(g1)
}
