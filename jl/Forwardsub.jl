function fwdsub(L,y)
    # L is a lower triangular matrix
    # y is a column vector
    cols = size(L,2)
    b = fill(0.0, cols)
    Lf = convert(Array{Float64}, L)
    yf = convert(Array{Float64}, y)
    for i in range(1, stop=cols)
        if i==1
            b[i]=yf[i]/Lf[i,i]
        else
            yf[i] = yf[i] - transpose(Lf[i,:])*b
            b[i] = yf[i]/Lf[i, i]
        end
    end
    return b
end