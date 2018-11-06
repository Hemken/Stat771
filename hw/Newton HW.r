x <- seq(-4.5, 4.5, length.out=91)
y <- seq(-4.5, 4.5, length.out=91)
support <- outer(x,y)

# F1 <- (1.5-x+x*y)*(y-1)+(2.25-x+x*y^2)*(y^2-1)+(2.625-x+x*y^3)*(y^3-1)
# F2 <- (1.5-x+x*y)*(x)+(2.25-x+x*y^2)*(2*x*y)+(2.625-x+x*y^3)*(3*y^2*x)

Fone <- function(x,y)(1.5-x+x*y)*(y-1)+(2.25-x+x*y^2)*(y^2-1)+(2.625-x+x*y^3)*(y^3-1)
Ftwo <- function(x,y)(1.5-x+x*y)*(x)+(2.25-x+x*y^2)*(2*x*y)+(2.625-x+x*y^3)*(3*y^2*x)

F1 <- outer(x,y,Fone)
F2 <- outer(x,y,Ftwo)

contour(x,y,F1)
contour(x,y,F2, add=TRUE)

f1y <- function(x,y) 1.5-6*x-2*x*y-6*x*y^2+4*x*y^3+6*x*y^5+3*2.625*y^2+4.5*y
f1x <- function(x,y) y^6+y^4-2*y^3-y^2-2*y+3

xiter <- rep(NA, 100)
xiter[1] <- x <- -5
for (i in 2:100){
  xiter[i] <- x - Fone(x,0)/f1x(x,0)
  x <- xiter[i]
}
print(xiter)
print(x)
Fone(x,0)
