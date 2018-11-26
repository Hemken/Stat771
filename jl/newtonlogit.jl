function newtonlogit(S,∇S, B, alpha; 
        Y=Y, X=X[:,1:7], counts=X[:,8], ϵ=1e-8, maxiter = 25)
  i = 0
  S₀ = S(B, Y, X, counts)
  ∇S₀ = ∇S(B, Y, X, counts)
  while norm(∇S₀) > ϵ && i <= maxiter
    i += 1
    #Search Direction
    ∇S₀ = ∇S(B, Y, X, counts)
    p = ∇S₀'
    #Step Length
    α = alpha #, evals = backtrack(S,∇S₀,p,X)
    #Update Parameter
    B += α*p
    #Update Gradient
    S₀ = S(B, Y, X, counts)
  end

  return B, i
end