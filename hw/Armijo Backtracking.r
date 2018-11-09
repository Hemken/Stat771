Fx <- function (x) {3 - 2*x^2} # an quadratic function with two roots
f <- function (x) {-4*x} # the derivative function of F

curve(Fx, -1.5, 1.5)
abline(c(0,0))

x <- 0.1 # note we cannot start at x==0
Fx(x)  # response at starting x
f(x)   # slope at starting x
-sign(f(x)) # direcction
abs(Fx(x)/f(x)) # step size
Fx(x-Fx(x)/f(x)) 
Fx(x)-Fx(x)/f(x)
for (a in seq(0.9, 1, 0.001)) {
  print(c(a, Fx(x-a*Fx(x)/f(x)), 
  (Fx(x)-a*Fx(x)/f(x))^2))
}

alpha = 5    # blows up
alpha = 2.5  # blows up slowly
alpha = 2    # oscillates
alpha = 1.25 # converges in 17
alpha = 1.1  # converges in 12
alpha = 1    # converges in 7
alpha = 0.9  # converges in 13
alpha = 0.5  # converges in 34

iter = 1
alpha = 2    # converges in 7
x <- 0.1 # note we cannot start at x==0
iterlog <- data.frame(iter=rep(NA,10),x=rep(NA,10),Fx=rep(NA,10), f=rep(NA,10))
while (abs(Fx(x))> 1e-8 && iter<=10) {
  iterlog[iter,1] <- iter
  iterlog[iter,2] <- x
  iterlog[iter,3] <- Fx(x)
  iterlog[iter,4] <- f(x)
  
  x <- x - alpha*Fx(x)/f(x)  # Newton's Method
  iter<-iter+1
}
print(iterlog)

# at alpha <- 2
alpha = 2
x <- 0.1 # arbitrary
Fx(x)/f(x) + Fx(x-alpha*Fx(x)/f(x))/f(x-alpha*Fx(x)/f(x)) == 0
# We need Fx(x)/f(x) + Fx(x-alpha*Fx(x)/f(x))/f(x-alpha*Fx(x)/f(x)) <= 0
# in order to converge
alpha = 1.9
Fx(x)/f(x) + Fx(x-alpha*Fx(x)/f(x))/f(x-alpha*Fx(x)/f(x)) < 0
alpha = 2.1
Fx(x)/f(x) + Fx(x-alpha*Fx(x)/f(x))/f(x-alpha*Fx(x)/f(x)) < 0

for (alpha in seq(1.5,2.5,0.1)) {
  print(Fx(x)/f(x) + Fx(x-alpha*Fx(x)/f(x))/f(x-alpha*Fx(x)/f(x)) >= 0)
}

# x <- 0.1 # note we cannot start at x==0
# alpha <- 1
# Fx0 <- Fx(x)
# p <- - Fx(x)/f(x) # search direction and step size at current point
# relax <- p
# 
# Fx1 <- Fx(x-p)
# 0.5*Fx1^2
# 0.5*Fx0^2+alpha*relax
# 
# Fx(x); norm(as.matrix(Fx(x)),"F")^2
# Fx(x-p); norm(as.matrix(Fx(x-p)),"F")^2
# Fx(x-p) <= Fx(x)
# 
iter <- 1
alpha0 <- 2    # converges in 7
x <- 0.1 # note we cannot start at x==0
iterlog <- data.frame(iter=rep(NA,10),x=rep(NA,10),
                      Fx=rep(NA,10), f=rep(NA,10), alpha=rep(NA,10))
while (abs(Fx(x))> 1e-8 && iter<=10) {
  alpha <- alpha0
  while (Fx(x)/f(x) + Fx(x-alpha*Fx(x)/f(x))/f(x-alpha*Fx(x)/f(x)) >= 0) {
    alpha <- 0.5*alpha
    print(alpha)
    print(Fx(x)/f(x) + Fx(x-alpha*Fx(x)/f(x))/f(x-alpha*Fx(x)/f(x)) >= 0)
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
