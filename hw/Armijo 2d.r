Fx <- function (x) {
  x1 <- x[1] + x[2]
  x2 <- x[1] * x[2]
  return(c(x1,x2))
}

Fx(c(1,1))

fx <- function (x) {
  x11 <- 1
  x21 <- 1
  x12 <- x[2]
  x22 <- x[1]
  J <- matrix(c(x11,x12,x21,x22), ncol=2)
  return(J)
}

fx(c(1,1))

x0 <- c(4,1)
Fx(x0); fx(x0)

-fx(x0)%*%Fx(x0)

iter     = 1
alpha    = 5 
x0 <- c(4,1)
while (norm(as.matrix(Fx(x0)), "F")> 1e-8 && iter<=35) {
  p <- solve(fx(x0),Fx(x0))
  while (norm(as.matrix(Fx(x0)-alpha*p), "F") < norm(as.matrix(Fx(x0-alpha*p)),"F")) {
    alpha <- 0.5*alpha
    print("alpha reset")
  }
  print(alpha)
  x0 <- x0 - alpha*p  # Newton's Method
  # print(iter)
  iter<-iter+1
}
x0; (iter-1)
