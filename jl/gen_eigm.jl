function gen_eigm(d)
    # d is a vector of eigenvalues
    # returns a matrix with those eigenvalues
    Σ = diagm(0 => d) # diagonal matrix of eigenvalues
    n,m = size(Σ)
    V = rand(n,m)     # random matrix
    A = inv(V)*Σ*V        # random A with given eigenvalues
    return A
end