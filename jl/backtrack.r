backtrack <- function (Fx,J,p,x0,alpha0 = 1.0, rho = 0.5, c = 0.2) {
  # Fx = function
  # J  = Jacobian of Fx0
  # p  = search direction
  # x0 = current guess
  # alpha0 = initial alpha
  # rho= proportion to backtrack
  # c  = relaxation parameter
  Fx0 <- Fx(x0)
#  print("x0"); print(x0)
#  print("Fx0"); print(Fx0)
  # print("Jacobian") ;print(J)
  # print(p)
  relaxSlope <- c*t(p)%*%(t(J)%*%Fx0)
  
  alpha <- alpha0
  Fxtry <- Fx(x0 + alpha*p)
#  print(Fx0)
#  print(Fx0%*%Fx0)
  condition <- 0.5*(Fx0%*%Fx0) + alpha*relaxSlope[1]
  
  while (0.5*(Fxtry %*% Fxtry) > condition) {
    alpha = rho*alpha
    Fxtry = Fx(x0 + alpha*p)
    condition <- 0.5*(Fx0%*%Fx0) + alpha*relaxSlope[1]
  }
  return(alpha)
}
