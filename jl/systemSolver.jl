"""
Name: systemSolver \n
Author: Vivak Patel \n
Date: June 17, 2016 \n \n

Description: module implements newton's method with line search for solving a system of equation. \n
Exports: newtonMethod \n
For method specific information use the help feature.
"""
module systemSolver

export newtonMethod

"""
Name: newtonMethod \n
Description: implements newton's method for solving a system of equations with a line search (backtracking) \n
INPUTS: \n
1. S :: Function, evaluates the system
2. ∇S :: Function, evaluates the Jacobian of the system
3. X :: Array{Float64,1}, starting point
4. ϵ :: Float64, tolerance. Defaults to 10e-8 \n

KEYWORD INPUTS: \n
1. maxIter :: Int64, maximum number of iterations. Defaults to 25\n

OUTPUTS: \n
1. X :: Array{Float}, solution to the system of equations within tolerance
2. funcEvals :: Int64, number of function evaluations
3. jacobEvals :: Int64, number of jacobian evaluations
"""
function newtonMethod(S,∇S, X, ϵ=10e-8; maxIter = 25)
  funcEvals = 0
  jacobEvals = 0

  S₀ = S(X)
  funcEvals += 1
  while dot(S₀,S₀) > ϵ && jacobEvals <= maxIter
    #Search Direction
    ∇S₀ = ∇S(X)
    jacobEvals += 1
    p = -(∇S₀\S₀)

    #Step Length
    α, evals = backtrack(S,∇S₀,p,X)
    funcEvals += evals

    #Update Parameter
    X += α*p
    S₀ = S(X)
    funcEvals += 1
  end

  return X, funcEvals, jacobEvals
end

"""
Name: backtrack \n
Description: implements backtracking line search on 0.5‖S‖² \n
INPUTS: \n
1. S (Function), evaluates the system
2. ∇S₀ (Function), Jacobian evaluated at a point X₀
3. p (Array{Float64,1}), search direction
4. X₀ (Array{Float64,1}), starting point
5. α₀ (Float64), starting search length defaults to 1
6. ρ (Float64), decay of search length, defaults to 0.5
7. c (Float64), slope relaxation value, defaults to 0.2 \n

OUTPUTS: \n
1. α (Float64), step length
2. funcEvals (Int64), number of function evaluations
"""
function backtrack(S,∇S₀,p,X₀,α₀ = 1.0, ρ = 0.5, c = 0.2)
  S₀ = S(X₀)
  relaxSlope = c*p'*(∇S₀'*S₀)
  condition(step) = 0.5*dot(S₀,S₀) + step*relaxSlope[1]

  α = deepcopy(α₀)
  Sₐ = S₀ #Ensures appropriate scope of Sₐ
  try
    Sₐ = S(X₀ + α*p)
  end
  funcEvals = 2
  while 0.5*dot(Sₐ,Sₐ) > condition(α)
    α = ρ*α
    try
      Sₐ = S(X₀ + α*p)
    catch
      continue
    end
    funcEvals += 1
  end

  return α, funcEvals
end

end #END OF MODULE
