function armijo(F,F₀,d,x₀,Y,X,counts; α₀ = 7.0, ρ = 0.5, z = 1e-4)
    x₊ = x₀ + α₀*d
    F₊ = F(x₊, Y, X, counts)
    i = 0
    while abs(F₊) > (1-z)*α₀*abs(F₀+ α₀*norm(d)) && i<5
        i +=1
        α₀ = ρ*α₀
        x₊ = x₀ + α₀*d
        F₊ = F(x₊, Y, X, counts)
    end
    return α₀
end