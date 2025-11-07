using LinearAlgebra
using NLsolve
using Pkg
using Polynomials
using Plots
using LaTeXStrings

#First we want to set up a function which constructs the Suryanarayana polynomial from the free variables.
#That is
#$$
#p'(x) = c\prod_{k=1}^{L} (x-l_k) \prod_{k=1}^{R} (x-r_k)
#$$
# K+2 free variables:
# r₁, r₂, ... rₖ, c₁ (equation c), c₂ (integration constant)

function p_ansatz(roots, c₁, c₂)
    x = Polynomial([0, 1])
    p_prim = c₁ * reduce(*, [x-rk for rk in roots]) * reduce(*, [x-(1-rk) for rk in roots])
    p = integrate(p_prim, c₂)
    return p, p_prim
end


# the purpose of the following function is to return the residual of the conditions for
# given free variables

function F!(conditions, freevar, lumo = 0.4)
    p, p_prim = p_ansatz(freevar[1:end-2], freevar[end-1], freevar[end])
    conditions[1] = p(freevar[1]) #p_prim(freevar[1])
    conditions[2] = p(lumo) - p(freevar[2])
    conditions[3] = p(freevar[3])
    conditions[4] = p(lumo) - p(0)
    conditions[5] = p(0.5) - 0.5
end

# function F2!(conditions, freevar, lumo = 0.40)
#     p, p_prim = p_ansatz(freevar[1:end-2], freevar[end-1], freevar[end])
#     conditions[1] = p_prim(freevar[1])
#     conditions[2] = p_prim(freevar[2])
#     conditions[3] = p_prim(freevar[3])
#     conditions[4] = p(lumo) - p(0) #corresponds to p(r_k+1) = p(lumo)
#     conditions[5] = p(0.5) - 0.5
# end

freevar0 = [collect(1/8:1/8:0.5-1/8); [-1000.0, 0.3]]
freevar0[end] += -p(0.5)+0.5;  # Adjust integration constant to match one condition exactly
lumo = 0.4
Ftmp!(x,y) = F!(x,y, lumo)
sol=nlsolve(Ftmp!,copy(freevar0))


# Continuation strategy:
#for i=2:10
#    global lumo += 9*10.0^(-i)
#    Ftmp!(x,y) = F!(x,y, lumo)
#    global sol = nlsolve(Ftmp!,sol.zero)
#end

p, p_prim = p_ansatz(sol.zero[1:3], sol.zero[4], sol.zero[5])
XX = [tmp for tmp in 0.0:0.01:1]   
plot(XX, p.(XX), label="Suryanarayana nlsolve", legend=:right)
plot!([-0.2, 1.0], [p(lumo), p(lumo)], linestyle=:dash, label=L"p(\lambda_{\textrm{lumo}})")
plot!([-0.2, 1.0], [p(1-lumo), p(1-lumo)], linestyle=:dash, label=L"p(\lambda_{\textrm{homo}})")
#plot!([sol.zero[1]-0.1, sol.zero[1]+ 0.1], [p(sol.zero[1]), p(sol.zero[1])], legend=:topleft, label="")
#plot!([sol.zero[2]-0.1, sol.zero[2]+ 0.1], [p(sol.zero[2]), p(sol.zero[2])], label="")
#plot!([sol.zero[3]-0.1, sol.zero[3]+ 0.1], [p(sol.zero[3]), p(sol.zero[3])], label="")

#scatter!([sol.zero[1:3]], [p.(sol.zero[1:3])], label="roots")

savefig("Suryanarayana_deg7.pdf")