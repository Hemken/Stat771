function backsub(U,y)
    # U is an upper triangular matrix
    # y is a column vector
    cols = size(U,2)
    b = fill(0.0, cols)
    Uf = convert(Array{Float64}, U)
    yf = convert(Array{Float64}, y)
    for i in range(cols, stop=1, step=-1)
        if i==cols
            b[i]=yf[i]/Uf[i,i]
        else
            yf[i] = yf[i] - transpose(Uf[i,:])*b
            b[i] = yf[i]/Uf[i, i]
        end
    end
    return b
end