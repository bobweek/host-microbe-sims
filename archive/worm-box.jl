include("worm-box-defs.jl")

n = 100 # number of microbe species

λ = 10
ψ₀ˢ = 1.0
Cₑˢ = 1.0-1e-5/n
Cₐˢ = 1.0-2e-5/n
Nˢ = fill(1000,n)

E = Environment(λ=fill(λ,n),ψ₀=fill(ψ₀ˢ,n),Cₑ=fill(Cₑˢ,n),Cₐ=fill(Cₐˢ,n),N=Nˢ)

using Distributions

sd = 1e-8 # stddev in microbial effects on worm fitness
nml = Normal(0,sd)

f = rand(nml,n)
ε = 0.01
ψ₀ʰ = 1.0
Cₑʰ = 1.0-1e-4/n
Cₐʰ = 1.0-2e-4/n
Nʰ = zeros(n)

h = Host(f=f,ε=fill(ε,n),ψ₀=fill(ψ₀ʰ,n),Cₑ=fill(Cₑʰ,n),Cₐ=fill(Cₐʰ,n),N=Nʰ)

H = LinkedList{Host}()
K = 100
for k in 1:K
	push!(H,deepcopy(h))
end

sys = System(E=deepcopy(E),H=deepcopy(H),K=K,K₀=K,n=n)


model = Model(sys=deepcopy(sys))

hist = run_model(model)

hist[999].S.N

@save "hist.jld2" hist

@load "hist.jld2" hist

hist[999].S.N