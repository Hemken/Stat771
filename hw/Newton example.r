x <- seq(-4.5, 4.5, length.out=91)
Fx <- function(x) (x^3-2*(x^2)+3*x-2)
fx <- function(x) (3*x^2-4*x+3)

plot(Fx(x)~x)

xiter <- rep(NA, 10)
xiter[1] <- x <- -5
for (i in 2:10){
  xiter[i] <- x - Fx(x)/fx(x)
  x <- xiter[i]
}

print(xiter)
print(x)
Fx(x)
