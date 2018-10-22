#Algorithm to Generate Decimal Representation

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

#Algorithm to Generate a Base Representation

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

decimal = 16.625

#Example 1
base, digits = 2, 4
d, e = getDigits(decimal, base, digits)
printstyled("EXAMPLE 1:\n", color=:red)
println("$decimal has a expansion in base $base with precision $digits \nof $d with exponent $e\n")
approx = getDecimal(base,d,e)
println("The decimal representation of $base-ary representation is $approx.\n")

#Example 2
base, digits = 2, 8
d, e = getDigits(decimal, base, digits)
printstyled("EXAMPLE 2:\n", color=:red)
println("$decimal has a expansion in base $base with precision $digits \nof $d with exponent $e\n")
approx = getDecimal(base,d,e)
println("The decimal representation of $base-ary representation is $approx.\n")


#Example 3
base, digits = 3, 8
d, e = getDigits(decimal, base, digits)
printstyled("EXAMPLE 3:\n", color=:red)
println("$decimal has a expansion in base $base with precision $digits \nof $d with exponent $e\n")
approx = getDecimal(base,d,e)
println("The decimal representation of $base-ary representation is $approx.\n")

#Half Precision Values

printstyled("Nonnormalized Smallest Positive 16 Bit Float\n", color=:red)
num1 = nextfloat(Float16(0.0))
num2 = 2.0^(-24)
println("$num1 == 2.0^(-24) => $(num1 == num2)\n")

printstyled("Largest Positive 16 Bit Float\n", color=:red)
num1 = prevfloat(Inf16)
num2 = getDecimal(2, ones(Int64,11),15)
println("$num1 == getDecimal(2, ones(Int64,11),15) => $(num1 == num2)\n")

printstyled("Example representation\n", color=:red)
println("Float16(1.391) has representation\n$(bitstring(Float16(1.391)))\n")

printstyled("Note: Binary Representation\n", color=:red)
d,e = getDigits(1.391, 2, 11)
println("Digits: $d")
println("Exponent: $e")

#Single Precision Values

printstyled("Nonnormalized Smallest Positive 32 Bit Float\n", color=:red)
num1 = nextfloat(Float32(0.0))
num2 = 2.0^(-149)
println("$num1 == 2.0^(-149) => $(num1 == num2)\n")

printstyled("Largest Positive 32 Bit Float\n", color=:red)
num1 = prevfloat(Inf32)
num2 = getDecimal(2, ones(Int64,24),127)
println("$num1 == getDecimal(2, ones(Int64,24),127) => $(num1 == num2)\n")

#Floating Point Approximation Errors

#Example 1
printstyled("EXAMPLE 1:\n", color=:red)
val1 = 1f8
pert1 = 4f0
println("32 Bit: $val1 == $val1 + $pert1?
    $(val1+pert1 == val1)\n")

#Example 2
printstyled("EXAMPLE 2:\n", color=:red)
val1 = 1f3
pert1 = 4f0
pert2 = 1f-5
println("32 Bit: $val1 == $val1 + $pert1? $(val1+pert1 == val1)\n")
println("32 Bit: $val1 == $val1 + $pert2? $(val1+pert2 == val1)\n")

# Values corresponding to 1 ULP

#Example 3
printstyled("EXAMPLE 3:\n", color=:red)
val1 = 1.438e3
eps1 = eps(val1)
println("The value corresponding to 1 ULP for $val1 is $eps1 \n")

# Machine Precision

#Example 4
printstyled("EXAMPLE 4:\n", color=:red)
val1 = 1.0
eps1 = eps(val1)
println("The value corresponding to 64 bit machine epsilon is $eps1 \n")
val1 = 1f0
eps1 = eps(val1)
println("The value corresponding to 32 bit machine epsilon is $eps1 \n")
val1 = Float16(1.0)
eps1 = eps(val1)
println("The value corresponding to 16 bit machine epsilon is $eps1 \n")

# Addition Relative Error

#Example 4
printstyled("EXAMPLE 4:\n", color=:red)
val1 = 1.1
val2 = 0.1
val3 = 1.2
println("$val1 + $val2 = $val3, but in floating point arithmetic, $val1 + $val2 = $(val1 + val2), which is a relative error of $(((val1 + val2)- val3)/val3).\n")
println("Note, twice 64 bit machine epsilon is $(2*eps(1.0))")
