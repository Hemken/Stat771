function backtrack(S,∇S₀,p,X₀,Y, X, counts; α₀ = 7.0, ρ = 0.5, c = 0.2)
  S₀ = S(X₀, Y, X, counts)
  #  println(length(c*p'*(∇S₀'*S₀))) # Vivak claims this is a scalar?
  relaxSlope = norm(c*p'*(∇S₀'*S₀))
  condition(step) = 0.5*dot(S₀,S₀) + step*relaxSlope
  # Vivak Patel's condition, normed
  α = deepcopy(α₀)
  Sₐ = S₀ #Ensures appropriate scope of Sₐ
  try
    Sₐ = S(X₀ + α*p, Y, X, counts)
  catch
  end
  while 0.5*dot(Sₐ,Sₐ) > condition(α)
    α = ρ*α
    try
      Sₐ = S(X₀ + α*p, Y, X, counts)
    catch
      continue
    end
  end

  return α
end