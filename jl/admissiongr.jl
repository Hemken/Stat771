function admitgr(B,Y,X,counts)
    Yf = convert(Array{Float64,1}, Y)
    gr = zeros(24,1)
    for i in 1:24
        gr[i,1] = (Yf[i].-1/(1+exp(-dot(X[i,:],B)))).*counts[i]
    end
    return((X[:,:]'*gr)./sum(counts))
end