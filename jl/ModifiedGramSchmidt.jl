using LinearAlgebra
function mgs(X)
    n = size(X,1)
    p = size(X,2)
    Q = zeros(n, p)
    R = zeros(p,p)
    for i = 1:p # columns
        Q[:,i] = X[:,i] # copy next vector
        for j = 1:(i-1) # rows, use previous vectors one at a time
            R[j,i] = (Q[:,j])'*Q[:,i] # build next R value
            Q[:,i] = Q[:,i] - R[j,i]*Q[:,j]
#            display(Q)
        end
        R[i,i] = norm(Q[:,i]) # normalizing constant
        Q[:,i] = Q[:,i]./R[i,i] # normalize this vector
    end
    return Q, R
end