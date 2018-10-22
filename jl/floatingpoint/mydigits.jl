function mydigits(dec::Float64, base::Int64, prec::Int=64)
    order = log(dec)/log(base)
#    println("order=", order)
    exp = Int(floor(abs(order))*sign(order))
#    println("exp=", exp)
    d = zeros(Int,1,max(exp+1,prec))
    power = Int(exp)
    println(dec, " in base ", base, " to ", prec, " digits:")
    for i = 1:prec
#        println("i=", i, " power=",power)
        if dec/(float(base)^power) >= 1
            d[i] = floor(dec/(float(base)^power))
            dec = dec - d[i]*(float(base)^power)
        end
        power -= 1
    end
    println(d,"exp", exp)
    if dec < 1e-15
        println("exactly")
    else
        println("approx.")
    end
    println("")
end
