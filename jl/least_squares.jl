function least_squares_dh(A,b)
    """
    Solves the linear squares regression problem given
    the coefficient matrix A and the constant vector b
    """

    n,m = size(A)

    # Get QR decomposition
    F = qr(A, Val(false))

    # Rank = m

    # Calculate c
    c = (F.Q'*b)
    c1 = c[1:m,:]
    c2 = c[m+1:n,:]

    # Calculate x
    x = backsub(F.R,c1)

    return x, norm(c2)
end