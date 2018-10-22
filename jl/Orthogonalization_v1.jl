
#Solving Linear System
using LinearAlgebra
function consistentLS(A,b)
    """
    Solves a consistent linear system given
    the coefficient matrix A and the constant
    vector b. Assumes A is consistent.
    """
    n, m = size(A)
    F = qr(A,Val(true))
    d = F.Q'*b
    c = F.R\d[1:m]
    return F.P*c
end

#Example 1
printstyled(
"EXAMPLE 1: Randomly Generated System\n", color=:red)
n, m = 10, 4
println("Dimension of A: $n by $m")
A = rand(10,4)
x = rand(4)
b = A*x
println("Error between consistentLS and  'truth':
    $(norm(consistentLS(A,b) - x))")

# Under determined Least Squares
function underLS(A,b; ϵ = 1e-14)
    """
    Solves an underdetermined linear system given
    the coefficient matrix A and the constant
    vector b. Returns the least norm solution.
    """
    n, m = size(A)
    s = min(n,m)
    F = qr(A, Val(true))
    
    #Compute rank approximation r
    Rtrm = F.R[1:s,1:s]
    r = maximum(findall(abs.(diag(Rtrm)) .>= ϵ))
    l = m - r
    
    #Generate R and S
    R, S = F.R[1:r,1:r], F.R[1:r,r+1:end]
    d, P = R\(F.Q'*b)[1:r], R\S
    z2 = consistentLS(P'*P + Matrix(I,l,l), P'*d)
    z1 = d - P*z2
    return F.P*vcat(z1,z2)
end

# Example 2

printstyled(
"EXAMPLE 2: Fat Matrix\n", color=:red)
n, m = 4, 10
A = rand(n,m)
b = rand(n)
println("A is an $n by $m matrix")
println("Error between underLS and 'truth':
    $(norm(underLS(A,b) - A\b))")

# Example 3

printstyled(
"EXAMPLE 3: Repeating Tall Matrix\n", color=:red)
n, m = 5,10
A = rand(n,m)
b = rand(n)
A = vcat(A,A,A,A)
b = vcat(b,b,b,b)
println("A is an $(size(A,1)) by $(size(A,2)) matrix.")
println("Error between underLS and 'truth':
    $(norm(underLS(A,b) - A\b))")


# Householder Reflection
function householder(a)
    """
    Computes the householder reflection 
    given a nonzero vector a.
    """
    nrm_a = norm(a,2)
    nrm_a == 0 && error("Input vector is zero.")
    
    d = length(a)
    v = copy(a)
    v[1] = v[1] - nrm_a
    H = Matrix(I,d,d) - (2/dot(v,v))*v*v'
    return H
end

#Example of Householder Reflection

printstyled(
"EXAMPLE 4: Householder Reflection\n", color=:red)
println("Householder reflection for random 10
dimensional vector applied to the vector.")
a = rand(10)
H = householder(a)
println("Norm excluding first element:
    $(norm((H*a)[2:end]))")


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

#Example of Givens Rotations

printstyled(
"EXAMPLE 5: Givens Rotation\n", color=:red)
println("Givens Rotations for random 10
dimensional vector applied to the vector.")
n = 10
a = rand(n)
for i = 2:n
    G = givens(a,1,i)
    a = G*a
end
println("Norm excluding first element:
    $(norm(a[2:end]))")


# Example of Large Scale Linear Regression

# Generate Generic Linear Regression
m = 10
n = 100000
x₀ = randn(m)

# Generate Data
using Random
Random.seed!(3940)

function data(x)
    a = randn(m)
    b = dot(a,x) + randn()
    return hcat(a',b)
end

# Incremental QR

RC = zeros(m+1,m+1)
RC[1,:] = data(x₀)
for i = 2:n
    ind = min(i,m+1)
    RC[ind,:] = data(x₀)
    for j = 1:ind
        G = givens(RC[:,j],j,ind)
        RC = G*RC
    end
end

x = RC[1:m,1:m]\RC[1:m,end]

println("Error between estimate and truth:
$(norm(x - x₀))")


