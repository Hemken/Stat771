# Implemented in Julia 1.0
# Vivak Patel

#Import Data
using RDatasets
UCBAdmit = RDatasets.dataset("datasets","UCBAdmissions")

#Generate Observed Variables
Y = map(y -> y == "Admitted" ? 1 : 0, UCBAdmit[1])

#Generate Explanatory Variables
using LinearAlgebra
X0 = ones(nrow(UCBAdmit))
X1 = map(x -> x == "Female" ? 1 : 0, UCBAdmit[2])
function proc_dept(val)
    E = Matrix{Float64}(I,5,5)
    if val == "A"; return zeros(5); end
    if val == "B"; return E[:,1]; end
    if val == "C"; return E[:,2]; end
    if val == "D"; return E[:,3]; end
    if val == "E"; return E[:,4]; end
    if val == "F"; return E[:,5]; end
end

X2 = vcat(map(x -> proc_dept(x)', UCBAdmit[3])...)
X = hcat(X0,X1,X2,Float64[UCBAdmit[4]...])

#Output Comments
println("Variable `Y` has the observed variable.")
println("Variable `X` has the explanatory variables.")
