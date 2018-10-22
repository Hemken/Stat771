# QR Solution
function consistentLS(A,b)
    """
    Solves a consistent linear system given
    the coefficient matrix A and the constant
    vector b. Assumes A is consistent.
    """
    n, m = size(A)
    F = qrfact(A,Val{true})
    d = F[:Q]'*b
    c = F[:R]\d[1:m]
    return F[:P]*c
end

function underLSQR(A,b; ϵ = 1e-14)
    """
    Solves an underdetermined linear system given
    the coefficient matrix A and the constant
    vector b. Returns the least norm solution.
    Uses the QR decomposition.
    """
    n, m = size(A)
    s = min(n,m)
    F = qrfact(A, Val{true})

    #Compute rank approximation r
    Rtrm = F[:R][1:s,1:s]
    r = maximum(find(abs.(diag(Rtrm)) .>= ϵ))
    l = m - r

    #Generate R and S
    R, S = F[:R][1:r,1:r], F[:R][1:r,r+1:end]
    d, P = R\(F[:Q]'*b)[1:r], R\S
    z2 = consistentLS(P'*P + eye(l), P'*d)
    z1 = d - P*z2
    return F[:P]*vcat(z1,z2)
end

# SVD Solution
function underLSSVD(A,b; ϵ = 1e-14)
    """
    Solves a (possibly) underdetermined linear
    systems given the coefficient matrix A and
    constant vector b using SVD of A. Returns
    the least norm solution.
    """
    n, m = size(A)
    U, Σ, V = svd(A, thin=false)
    c = U'*b

    #Determine c₁ and solution
    z = n < m ? vcat(c./Σ,zeros(m-n)) : c./Σ #Σ is vector
    return V*z
end;

# Example of Solving Linear System

n, m = 10, 20
print_with_color(:red,
"EXAMPLE 1: SVD vs QR\n")

println(
"Computing minimum-norm LS\n")

A = rand(n,m)
x = rand(m)
b = A*x

x₁ = underLSSVD(A,b)
x₂ = underLSQR(A,b)

println(
"Residual for SVD:
    $(norm(A*x₁ - b))")
println(
"Residual for QR:
    $(norm(A*x₂ - b))")
println(
"Error SVD and QR:
    $(norm(x₁-x₂))")

# Example Norm Calculation
n, m = 10, 20
print_with_color(:red,
"EXAMPLE 2: 2-Norm\n")

println(
"Computing Matrix 2-Norm\n")

A = rand(n,m)

println("Julia's 2-Norm of A:
    $(norm(A,2))\n")

println("Max Singular Value of A:
    $(maximum(svdvals(A)))")


#Approximate Rank K Calculation

using Distributions

function inexactRankK(A,k,c,p)
    """
    Inexact rank k approximation to a matrix A.
    c is an integer betweek k and m.
    p is an m-tuple probabilities summing to 1.
    """
    dist = Categorical(p)
    cols_ind = rand(dist,c)

    #Divide each col by sqrt(c*probability)
    C = A[:,cols_ind]./sqrt.(c*p[cols_ind]')

    #Compute eigenvalue decomposition
    D, V = eig(C'*C)

    #Compute top k singular values
    S = sqrt.(D[end-k+1:end])

    #Take top k right singular vec of C
    Z = V[:,end-k+1:end]

    #Compute left singular vectors
    H = C*(Z./S')
    return H, S
end

#
# Examples
#
#Generate Random Matrices

function generateLowRank(n,m)
    m <= n &&
    error("First argument ($n) must be smaller than second ($m)")
    r = rand(1:n)
    U, _ = qr(rand(n,n))
    V, _ = qr(rand(m,m))
    Σ = hcat(diagm(vcat(100*rand(r),zeros(n-r))),zeros(n,m-n))
    return U*Σ*V', r
end

#Example 3
n, m = 30,500
A, r = generateLowRank(n,m)
print_with_color(:red,
"EXAMPLE 3: $n by $m matrix with rank $r\n\n")

k, c, p = max(r-10,1), 30, ones(m)/m
H, S = inexactRankK(A,k,c,p)
err = vecnorm(A - H*H'*A,2)^2
println("k=$k and c=$c, error: $err\n")

k, c, p = max(r-10,1), 60, ones(m)/m
H, S = inexactRankK(A,k,c,p)
err = vecnorm(A - H*H'*A,2)^2
println("k=$k and c=$c, error: $err\n")

k, c, p = max(r-10,1), 120, ones(m)/m
H, S = inexactRankK(A,k,c,p)
err = vecnorm(A - H*H'*A,2)^2
println("k=$k and c=$c, error: $err\n")

k, c, p = r, 30, ones(m)/m
H, S = inexactRankK(A,k,c,p)
err = vecnorm(A - H*H'*A,2)^2
println("k=$k and c=$c, error: $err\n")

k, c, p = r, 60, ones(m)/m
H, S = inexactRankK(A,k,c,p)
err = vecnorm(A - H*H'*A,2)^2
println("k=$k and c=$c, error: $err\n")

k, c, p = r, 120, ones(m)/m
H, S = inexactRankK(A,k,c,p)
err = vecnorm(A - H*H'*A,2)^2
println("k=$k and c=$c, error: $err\n")
