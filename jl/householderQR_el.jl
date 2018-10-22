function householderQR(A)
    # from Erika Lee
    n, m = size(A)
    a = A[1:end, 1]
    H1 = householder(a)
    Q = H1
    Ai = H1*A
    for i in 2:m
        Ai = Ai[2:end,2:end]
        a = Ai[1:end, 1]
        Hi = householder(a)
        Ai = Hi*Ai
        hi = hcat(Matrix{Float64}(I,i-1,i-1),zeros(i-1,n-i+1))
        hi_ = hcat(zeros(n-i+1,i-1),Hi)
        Q = Q*vcat(hi, hi_)
    end
    R = Q'A
    return Q, R
end