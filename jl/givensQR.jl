function givensQR(A)
    n,m = size(A)
    R = copy(A)
    G = Matrix{Float64}(I,n,n)
    for j in 1:m
        for i = n:-1:(j+1)
            a = R[:,j]
            Gi = givens(a,j,i)
            G = G*Gi'
            R = Gi*R
        end
    end
    return R, G
end