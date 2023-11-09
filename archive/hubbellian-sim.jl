include("hubbellian-defs.jl")

# model parameters
P = ModelParameters()

# build model and initial conditions
sim = make_sim(P)

# run simulation
hist = run_sim(sim)

# save simulated data
@save "hist.jld2" hist

@load "hist.jld2" hist