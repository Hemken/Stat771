function admithessian(B,Y,X,counts)
    Yf = convert(Array{Float64,1}, Y)
    H = zeros(24,24)
    for i in 1:24
        H[i,i] = (1/(1+exp(1).^(-dot(X[i,:],B))))*(1-(1/(1+exp(-dot(X[i,:],B))))).*counts[i]
    end
    return(X'*H*X./sum(counts))
end