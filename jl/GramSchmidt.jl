using LinearAlgebra
function gs(X)
    n = size(X,1)
    p = size(X,2)
    Q = zeros(size(X))
    R = zeros(p,p)
    for i = 1:p
        Q[:,i] = X[:,i] # copy next vector
        if i>1
            R[1:(i-1),i] = (Q[:,1:(i-1)])'*Q[:,i] # coefficients to remove prev vecs
            Q[:,i] = Q[:,i] - Q[:,1:(i-1)]*R[1:(i-1),i] # new vec
        end
        R[i,i] = norm(Q[:,i])   # normalizing constant for this vector
        Q[:,i] = Q[:,i]./R[i,i] # normalize this vector
    end
    return Q, R
end