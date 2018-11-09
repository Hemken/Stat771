Fx <- function (x) {3 - 2*x^2} # an quadratic function with two roots
f <- function (x) {-4*x} # the derivative function of F

curve(Fx, -1.5, 1.5)
abline(c(0,0))

x <- 0.1 # note we cannot start at x==0
Fx(x)  # response at starting x
f(x)   # slope at starting x
-sign(f(x)) # direction
abs(Fx(x)/f(x)) # step size

# without backtracking
alpha = 5    # blows up
alpha = 2.5  # blows up slowly
alpha = 2    # oscillates
alpha = 1.25 # converges in 17
alpha = 1.1  # converges in 12
alpha = 1    # converges in 7
alpha = 0.9  # converges in 13
alpha = 0.5  # converges in 34

max.iter <- 35
iter     = 1
alpha    = 1    # converges in 7
x        <- 0.1 # note we cannot start at x==0

iterlog <- data.frame(iter=rep(NA,max.iter),x=rep(NA,max.iter),
                      Fx=rep(NA,max.iter), f=rep(NA,max.iter))

while (abs(Fx(x))> 1e-8 && iter<=35) {
  iterlog[iter,1] <- iter
  iterlog[iter,2] <- x
  iterlog[iter,3] <- Fx(x)
  iterlog[iter,4] <- f(x)
  
  x <- x - alpha*Fx(x)/f(x)  # Newton's Method
  iter<-iter+1
}
print(iterlog)

# With backtracking
iter   <- 1
alpha0 <- 1
x      <- -0.1 # note we cannot start at x==0

iterlog <- data.frame(iter=rep(NA,10),x=rep(NA,10),
                      Fx=rep(NA,10), f=rep(NA,10), alpha=rep(NA,10))

while (abs(Fx(x))> 1e-8 && iter<=10) {
  alpha <- alpha0
  p <- Fx(x)/f(x)
  while (abs(Fx(x)) + alpha*p < abs(Fx(x-alpha*p))) {
    alpha <- 0.5*alpha
    print(alpha)
  }
  iterlog[iter,1] <- iter
  iterlog[iter,2] <- x
  iterlog[iter,3] <- Fx(x)
  iterlog[iter,4] <- f(x)
  iterlog[iter,5] <- alpha
  
  x <- x - alpha*Fx(x)/f(x)  # Newton's Method
  alpha <- alpha0
  iter<-iter+1
}
print(iterlog)
