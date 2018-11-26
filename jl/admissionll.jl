function admitll(B, Y, X, counts)
    ll = zeros(24,1)
    for i = 1:24
        #println(dot(X[i,:],B))
        ll[i] = counts[i]*(Y[i]*dot(X[i,:],B)-log(1+exp(dot(X[i,:],B))))
    end
    return(sum(ll))
end