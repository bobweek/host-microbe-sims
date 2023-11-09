include("hubbellian-defs.jl");

using Gadfly, DataFrames, GLM, Colors, Cairo, Fontconfig

nd = 5      # num dfe draw reps
ns = 200    # num sel exprmt reps per sim rep
np = 5      # num reps per par pair

vs = 0.0:0.2:1.0
ωs = 0.0:0.2:1.0

dat = DataFrame()
dat[!,"a²"] = zeros(0)
dat[!,"b"] = zeros(0)
dat[!,"v"] = zeros(0)
dat[!,"ω"] = zeros(0)

# for loops across par pars here
for v₀ in vs
    for ω₀ in ωs

        a² = Vector{Float64}()
        b  = Vector{Float64}()

        # model parameters
        P = ModelParameters(sd=1.0,L=5,t=v₀,ω=ω₀);

        # for loop across reps per par pair here

        for npi in 1:np

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
            b₀, a²₀ = coef(ols);

            push!(a²,a²₀)
            push!(b,b₀)

        end

        n = length(a²)
        thing = DataFrame(a²=a²,b=b,v=fill(v₀,n),ω=fill(ω₀,n))
        append!(dat,thing)

    end

end

CSV.write("sel_exprmnts.csv",dat);
