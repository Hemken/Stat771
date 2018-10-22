function getDecimal(base :: Int64, d :: Vector{Int64}, e :: Int64)
    """
    A function that takes a base, a vector of digits and an exponent,
    then returns a decimal represetntation of the number.
    """
    num = 0
    base = float(base)
    for i = 0:length(d)-1
        num += d[i+1]*base^(-i)
    end
    return num*base^e
end