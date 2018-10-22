function Hvec(u,x)
    return x - u*(u'*x)
end

function house_qr(A)
    """ Householder reflections for QR decomposition.
    % R, the upper triangular factor, and
    % U, the reflector generators for use by house_apply.
    from Moler 
    https://blogs.mathworks.com/cleve/2016/10/03/householder-reflections-and-the-qr-decomposition
    """

    m,n = size(A)
    U = zeros(m,n)
    R = copy(A)
    for j = 1:min(m,n)
        u = house_gen(R[j:m,j])
        U[j:m,j] = u
        R[j:m,j:n] = Hvec(u,R[j:m,j:n])
        R[j+1:m,j] = repeat([0], m-j)
    end
    return R, U
end

function house_apply(U,X)
    """ Apply Householder reflections.
     without actually computing Q.  continuing from Moler
     Calculates Q*X, from U returned by house_qr()
    """
    Z = copy(X)
    m,n = size(U)
    for j = n:-1:1
        Z = Hvec(U[:,j],Z)
    end
    return Z
end