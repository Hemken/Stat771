function gen_Axb(n)
    a = rand(n)
    A = [ones(n) a]
    e = randn(Float64, n)
#    display(e)
    x = [1,2]
    b = A*x + e
    return A, x, b
end