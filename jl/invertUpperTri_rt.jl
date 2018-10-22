function invertUpperTri(A)
    ## Get dimensions of A
    n, m = size(A)

    ## Setup empty array to hold result
    B = zeros(n,m)

    ## Fill out diagonal with inverse diagonal from A
    for k = 1:n
        B[k, k] = 1/A[k,k]
    end

    ## Starting in the lower right corner, fill out the rest of the matrix.
    for i = (n-1):-1:1
        for j = (i+1):n
            B[i,j] = -sum(A[i,(i+1):j].*B[(i+1):j,j])/A[i,i]
        end
    end

    return(B)
end