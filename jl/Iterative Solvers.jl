#Implementation of Jacobi
function jacobiIteration(A, b, x; iterMax=100)
    """Implements the Jacobi iteration.
    A is the coefficient matrix.
    b is the constant matrix.
    x is an initial guess."""

    D = Diagonal(A)
    f = D\b
    G = D\(D - A)
    for j=1:iterMax
        x = G*x + f
    end
    return x, maximum(abs.(eigvals(G)))
end

#Implementation of Gauss-Seidel

function gsIteration(A, b, x; iterMax = 100)
    """Implements the GS iteration.
    A is the coefficient matrix.
    b is the constant matrix.
    x is an initial guess."""
    DE = tril(A)
    F = -triu(A,1)
    G = DE\F
    f = DE\b
    for j=1:iterMax
        x = G*x + f
    end
    return x, maximum(abs.(eigvals(G)))
end


function generateProblem(n)
    A = randn(n,n)  + sqrt(n)*eye(n)
    x = randn(n)
    b = A*x
    return A, b, x
end

print_with_color(:red,
"EXAMPLE 1: n=10\n")
A, b, x = generateProblem(10)
z_j, ρ_j = jacobiIteration(A,b,x)
z_gs, ρ_gs = gsIteration(A,b,x)
println("Jacobi Iteration")
println("\tAbsolute Error: $(norm(x-z_j))")
println("\tSpect. Radius: $ρ_j")
println("Gauss Seidel Iteration")
println("\tAbsolute Error: $(norm(x-z_gs))")
println("\tSpect. Radius: $ρ_gs\n\n")

print_with_color(:red,
"EXAMPLE 2: n=50\n")
A, b, x = generateProblem(50)
z_j, ρ_j = jacobiIteration(A,b,x)
z_gs, ρ_gs = gsIteration(A,b,x)
println("Jacobi Iteration")
println("\tAbsolute Error: $(norm(x-z_j))")
println("\tSpect. Radius: $ρ_j")
println("Gauss Seidel Iteration")
println("\tAbsolute Error: $(norm(x-z_gs))")
println("\tSpect. Radius: $ρ_gs\n\n")

print_with_color(:red,
"EXAMPLE 3: n=100\n")
A, b, x = generateProblem(100)
z_j, ρ_j = jacobiIteration(A,b,x)
z_gs, ρ_gs = gsIteration(A,b,x)
println("Jacobi Iteration")
println("\tAbsolute Error: $(norm(x-z_j))")
println("\tSpect. Radius: $ρ_j")
println("Gauss Seidel Iteration")
println("\tAbsolute Error: $(norm(x-z_gs))")
println("\tSpect. Radius: $ρ_gs")

#Kaczmarz Update

function kaczmarzUpdate(b,a,x)
    """
    The Kaczmarz update, where b is a scalar,
    a is vector, and x is the current iterate.
    Returns the next iterate.
    """
    return x + ((b - dot(a,x))/dot(a,a))*a
end


#Cyclical Kaczmarz
function cyclicalKaczmarz(b, A, x;maxIter=5)
    """
    Estimates solution to Ax=b using Kaczmarz
    update by cycling through the rows of A.
    b is constant vector.
    A is the coefficient matrix.
    x is an initial iterate.
    maxIter is the number of cycles over A.
    """
    n = length(b)
    for i = 0:(maxIter*n)-1
        k = rem(i,n)+1
        x = kaczmarzUpdate(b[k],A[k,:],x)
    end
    return x
end


#Randomized Kaczmarz
using Distributions
function randomizedKaczmarz(b, A, x; maxIter=1000)
    """
    Estimates solution to Ax=b using randomized
    Kaczmarz, where probability of sampling a row
    depends on the sum of squares of the row.
    Samples are taken with replacement.
    b is a constant vector.
    A is the coefficient matrix.
    x is an initial iterate.
    maxIter is number of rows sampled.
    """
    n = length(b)
    p = map(λ -> sum(A[λ,:].^2), 1:n)/sum(A.^2)
    dist = Categorical(p)
    indices = rand(dist,maxIter)
    for k in indices
        x = kaczmarzUpdate(b[k],A[k,:],x)
    end
    return x
end

#Examples
print_with_color(:red,
"EXAMPLE 4: n=100\n")
n = 100
d = n
A = randn(n,d)
x₊ = randn(d)
b = A*x₊

x₀ = randn(d)
x_cyc = cyclicalKaczmarz(b, A, x₀, maxIter=5)
x_ran = randomizedKaczmarz(b, A, x₀, maxIter=5*n)

println("Cyclical Kaczmarz")
println("\tNumber of Iterations: 5")
println("\tResidual Norm: $(norm(A*x_cyc-b))")
println("Randomized Kaczmarz")
println("\tNumber of Iterations: 5000")
println("\tResidual Norm: $(norm(A*x_ran-b))\n")

x_cyc = cyclicalKaczmarz(b, A, x₀, maxIter=10)
x_ran = randomizedKaczmarz(b, A, x₀, maxIter=10*n)

println("Cyclical Kaczmarz")
println("\tNumber of Iterations: 10")
println("\tResidual Norm: $(norm(A*x_cyc-b))")
println("Randomized Kaczmarz")
println("\tNumber of Iterations: 10000")
println("\tResidual Norm: $(norm(A*x_ran-b))\n")

println("Convergence Factor:")
println("\t$(1 - minimum(svdvals(A))^2/sum(A.^2))")

##Conjugated Solver

#Gram-Schmidt Conjugation
function conjGS(A,V; ϵ = 1e-15)
    """Implements classical Gram-Schmidt
    A-conjugation over a set of vectors
    in the vector V. Returns A-conjugated
    vectors using set in V"""

    normA(x) = sqrt(dot(x,A*x))

    Q = Vector(length(V))

    counter = 1
    for j = 1:length(V)
        if counter == 1
            #Check if zero vector
            normA(V[j]) <= ϵ && continue

            #Otherwise add to Conjugated Set
            Q[1] = V[j]/normA(V[j])

            counter+= 1
        else
            #Project
            proj = map(λ -> dot(Q[λ],A*V[j])*Q[λ] ,1:counter-1)
            a = V[j] - sum(proj)

            #Check if zero vector
            normA(a) <= ϵ && continue

            #Otherwise add to Conjugated Set
            Q[counter] = a/normA(a)

            counter +=1
        end
    end
    return Q
end


print_with_color(:red,
    "EXAMPLE 5: Sanity Check\n\n")
A = eye(Float64,10)
V = map(j -> A[:,j],1:5)

Q = conjGS(A,V)

println("Difference between V and Q: $(norm.(Q-V)) \n\n")

print_with_color(:red,
    "EXAMPLE 6: Conjugate of Standard Basis\n\n")

srand(10)
A = randn(10,10); A = A'*A;
V = map(j -> eye(10)[:,j], 1:10)

Q = conjGS(A,V)
Q = hcat(Q...) #Make Columns of Matrix
println("Inner Products: $((Q'*A*Q)[:,1])\n\n")

print_with_color(:red,
    "EXAMPLE 7: Failed Gram-Schmidt\n\n")
A = eye(3)
δ = 1e-8
V = Vector{Float64}[
    Float64[1, δ, δ ],
    Float64[1, δ, 0],
    Float64[1, 0, δ]]

Q = conjGS(A,V)
Q = hcat(Q...)
println("Inner Products: $((Q'*A*Q)) \n\n")

#Gram-Schmidt Conjugate 2
function conjGS2(A,V; ϵ = 1e-15)
    """Implements Classical GS conjugation
    twice. (Giraud et al. 2005)
    """

    normA(x) = sqrt(dot(x,A*x))

    Q = Vector(length(V))

    counter = 1
    for j = 1:length(V)
        if counter == 1
            #Check if zero vector
            normA(V[j]) <= ϵ && continue

            #Otherwise add to Conjugated Set
            Q[1] = V[j]/normA(V[j])

            counter+= 1
        else
            #Project
            proj = map(λ -> dot(Q[λ],A*V[j])*Q[λ] ,1:counter-1)
            a = V[j] - sum(proj)

            #Check if zero vector
            normA(a) <= ϵ && continue

            #Otherwise Repeat Projection
            proj = map(λ -> dot(Q[λ],A*a)*Q[λ], 1:counter-1)
            a = a - sum(proj)

            #Add Entry to Conjugated Set
            Q[counter] = a/normA(a)

            counter +=1
        end
    end
    return Q
end

print_with_color(:red,
    "EXAMPLE 8: Conjugate of Standard Basis\n\n")

srand(10)
A = randn(10,10); A = A'*A;
V = map(j -> eye(10)[:,j], 1:10)

Q = conjGS2(A,V)
Q = hcat(Q...) #Make Columns of Matrix
println("Inner Products: $((Q'*A*Q)[:,1]) \n\n")

print_with_color(:red,
    "EXAMPLE 9: Failed Gram-Schmidt\n\n")
A = eye(3)
δ = 1e-8
V = Vector{Float64}[
    Float64[1, δ, δ ],
    Float64[1, δ, 0],
    Float64[1, 0, δ]]

Q = conjGS2(A,V)
Q = hcat(Q...)
println("Inner Products: $((Q'*A*Q)) \n\n")

#Solver

function conjSolv(A,b,x₀, V; ϵ = 1e-15)
    """Solves Ax=b using A-conjugated directions.
    A is symmetric, positive definite matrix.
    b is an arbitrary coefficient vector.
    x₀ is an initial guess.
    V is a vector of linearly ind. vectors.
    """

    Q = conjGS2(A,V)
    r₀ = A*x₀ - b

    α = r₀'*hcat(Q...)
    x = x₀
    for j = 1:length(α)
        x = x - α[j]*Q[j]
    end
    return x
end

print_with_color(:red,
    "EXAMPLE 10: Conjugated Direction Solver\n\n")

A = randn(10); A = A'*A;
x_true = randn(10)
b = A*x_true
V = map(j -> eye(10)[:,j],1:10)
x₀ = randn(10)

y = conjSolv(A,b,x₀,V)

println("Difference between truth and computed:
    \t $(norm(y - x_true)) ")

# Conjguated Gradient Solver
function conjGraSolv(A,b,x₀; ϵ = 1e-15)
    """Implementation of CG method.
    A is symmetric positive definite
    b is a constant vector.
    x0 is an initial guess."""

    #Update Iterate & Residual
    x = x₀
    r = b - A*x;
    p = r;
    rsold = r' * r;

    for i = 1:length(b)
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-15
            println("Breaking at iteration $i")
              break;
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
    return x
end

print_with_color(:red,
    "EXAMPLE 11: Conjugated Gradient Solver\n\n")

A = randn(100); A = A'*A + Diagonal(10*rand(100));
x_true = randn(100)
b = A*x_true
x₀ = randn(100)

y = conjGraSolv(A,b,x₀)

println("Difference between truth and computed:
    \t $(norm(y - x_true)) ")
