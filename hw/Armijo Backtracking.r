Fx <- function (x) {3 - 2*x^2} # an quadratic function with two roots
f <- function (x) {-4*x} # the derivative function of F

curve(Fx, -1.5, 1.5)

x <- 0.1 # note we cannot start at x==0

alpha = 5    # blows up
alpha = 2.5  # blows up slowly
alpha = 2    # oscillates
alpha = 1.25 # converges in 17
alpha = 1.1  # converges in 12
alpha = 1    # converges in 7
alpha = 0.9  # converges in 13
alpha = 0.5  # converges in 34
iter = 1
x <- 0.1 # note we cannot start at x==0
while (abs(Fx(x))> 1e-8 && iter<=35) {
  x <- x - alpha*Fx(x)/f(x)  # Newton's Method
  print(c(iter, x, Fx(x)))
  iter<-iter+1
}

x <- 0.1 # note we cannot start at x==0
alpha <- 1
Fx0 <- Fx(x)
p <- - Fx(x)/f(x) # search direction and step size at current point
relax <- p

Fx1 <- Fx(x-p)
0.5*Fx1^2
0.5*Fx0^2+alpha*relax

Fx(x); norm(as.matrix(Fx(x)),"F")^2
Fx(x-p); norm(as.matrix(Fx(x-p)),"F")^2
Fx(x-p) <= Fx(x)

