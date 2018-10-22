function getDigits(decimalNum, base :: Int64, digits :: Int64)
    """
    A function that takes a decimal representation of a number and
    returns its representation in a specific base up to a certain
    number of digits.
    """
    base = float(base)
    e = floor(Int64,log(base,decimalNum))
    d = zeros(Int64,digits)
    num = decimalNum/(base^e)
    for j = 1:digits
        d[j] = floor(Int64,num)
        num = (num - d[j])*base
    end

    return d, e
end