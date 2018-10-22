function consistentLS(A,b)
    """
    Solves a consistent linear system given
    the coefficient matrix A and the constant
    vector b. Assumes A is consistent.
    """
    n, m = size(A)
    F = qr(A,Val(true))
    d = F.Q'*b
    c = backsub(F.R, d[1:m])
    return F.P*c
end