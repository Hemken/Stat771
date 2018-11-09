backtrack <- function (Fx,f,p,x0,alpha0 = 1.0, rho = 0.5, c = 0.2) {
  Fx0 <- Fx(x0)
  fx0 <- f(x0)
  relaxSlope <- c*t(p)%*%(t(fx0)%*%Fx0)
    
  alpha <- alpha0
  #Fxtry <- Fx0
  Fxtry <- Fx(x0 + alpha*p)

  while (0.5*(Fxtry %*% Fxtry) > 0.5*(Fx0%*%Fx0) + alpha*relaxSlope[1]) {
    alpha = rho*alpha
    Fxtry = Fx(x0 + alpha*p)
  }
return(alpha)
}
