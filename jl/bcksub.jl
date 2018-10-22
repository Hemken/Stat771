function bcksub(U,y)
    # U is an upper triangular matrix
    # y is a matrix
    cols = size(U,2)
    ycols = size(y,2)
    b = fill(0.0, (cols, ycols))
    Uf = convert(Array{Float64}, U)
    yf = convert(Array{Float64}, y)
    for j in 1:ycols
        yj = yf[:,j]
        for i in range(cols, stop=1, step=-1)
            if i==cols
                b[i,j]=yj[i]/Uf[i,i]
            else
                yj[i] = yj[i] - transpose(Uf[i,:])*b[:,j]
                b[i,j] = yj[i]/Uf[i, i]
            end
        end
    end
    return b
end