function givens(a,i,j)
    """
    Computes the Givens Rotation for a
    vector a at indices i and j, where
    the index at j is set to zero.
    """
    d = length(a)
    (i > d || j > d) && error("Index out of range.")
    l = sqrt(a[i]^2 + a[j]^2)
    λ = a[i]/l
    σ = a[j]/l
    G = ones(d)
    G[i] = λ
    G[j] = λ
    G = diagm(0 => G)
    G[i,j] = σ
    G[j,i] = -σ
    return G
end