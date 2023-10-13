using Parameters, LinearAlgebra, Random, Distributions, DataFrames, CSV

function realβ(s,n)
    z = rand(Normal(),n)
    W = exp.(s.*z)
    w = normalize(W,1)
    Nₒ = rand(Multinomial(n,w))
    N̂ₒ = normalize(Nₒ,1)
    β = cov(Nₒ,z)/var(z)
    β̂ = cov(N̂ₒ,z)/var(z)
    return β, β̂
end

realβ(1,10000000)