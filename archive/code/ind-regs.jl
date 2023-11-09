include("hubbellian-defs.jl");

using Gadfly, DataFrames, GLM, Colors, Cairo, Fontconfig, Compose

nd = 5      # num dfe draw reps
ns = 20     # num sel exprmt reps per sim rep

ℓ = 1
ε = 1

# model parameters
P = ModelParameters(sd=1.0,L=5,t=ℓ,ω=ε);

# for loop across reps per par pair here

S = Vector{Float64}()
R = Vector{Float64}()

# rep dfe draws
for ndi in 1:nd

    # build model and initial conditions (baby parents)
    sim = make_sim(P);

    # rep sel exprmts
    for nsi in 1:ns

        # perform fitness experiment
        S₀, R₀ = fit_diffs(sim)
        push!(S,S₀)
        push!(R,R₀)

    end

end

# get slope (a²) and intercept (b)
df = DataFrame(S=S,R=R)
ols = lm(@formula(R ~ S), df)
b, m² = coef(ols);

m²ᵣ = round(m², digits=3)
m²ₛ = "m² = $m²ᵣ"

bᵣ = round(b, digits=3)
bₛ = "b = $bᵣ"

p = plot(df,x=:S,y=:R,Geom.point,
    intercept=[b], slope=[m²], Geom.abline(color="blue", style=:dash),
    Guide.xlabel("Selection Differential (S)"), Guide.ylabel("Response (R)"),
    Guide.annotation(compose(context(), Compose.text(-0.0025,0.0175,m²ₛ), fill("blue"))),
    Guide.annotation(compose(context(), Compose.text(-0.0025,0.0160,bₛ), fill("blue"))),
    Coord.cartesian(xmin=-0.005, xmax=0.02, ymin=-0.005, ymax=0.02), Guide.title("ℓ=$ℓ, ε=$ε"))

p |> PNG("hubbellian/evowibo/ℓ$ℓ-ε$ε.png")

# p = plot(df,x=:S,y=:R,Geom.point,
#     intercept=[b], slope=[m²], Geom.abline(color="blue", style=:dash),
#     Guide.xlabel("Selection Differential (S)"), Guide.ylabel("Response (R)"),
#     Guide.annotation(compose(context(), Compose.text(-0.0025,0.075,m²ₛ), fill("blue"))),
#     Guide.annotation(compose(context(), Compose.text(-0.0025,0.060,bₛ), fill("blue"))),
#     Coord.cartesian(xmin=-0.005, xmax=0.02, ymin=-0.1, ymax=0.1), Guide.title("ℓ=$ℓ, ε=$ε"))

# p |> PNG("hubbellian/evowibo/ℓ0-ε0.png")
