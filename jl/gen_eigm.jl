function gen_eigm(d)
    # d is a vector of eigenvalues
    # returns a matrix with those eigenvalues
    Σ = diagm(0 => d) # diagonal matrix of eigenvalues
    n,m = size(Σ)
    V = rand(n,m)     # random matrix
    Q,R = qr(V)       # orthonormalized
    A = Q*Σ*Q'        # random A with given eigenvalues
    return A
end