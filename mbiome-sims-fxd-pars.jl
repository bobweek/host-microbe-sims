#########################################################################
#                                                                       #
# use mbiome expl sims to do the following:                             #
#                                                                       #
#                                                                       #
# (1) look at evolution of ḡ, m̄, G, and M                               #
# 	as fct of S ∈ {10, 100, 1000} with ℓ = 1, κ = 0                     #
# 	as fct of Nₑ ∈ {10, 100, 1000} with ℓ = 1, κ = 0                    #
# 	as fct of κ ∈ {0, 0.25, 0.5, 0.75, 1} with ℓ + κ = 1 and Nₑ = 100   #
#	as fct of ℓ,κ ∈ {0, 0.25, 0.5, 0.75, 1} with S = 1000               #
#		similar to analytical host trait evo across trans modes         #
#                                                                       #
# (2) using same data, look at Corr(g,m)                                #
#                                                                       #
# (3) using same data, "fit" analytical model                           #
#	in particular, set analytical model pars (e.g., G and M)            #
#	to simulated data, compute Δz̄ from simulated data                   #
#	compare with analytically predicted Δz̄                              #
#                                                                       #
#########################################################################

include("defs.jl");

print("\nPackages Loaded\n")

function Nₑ(Nₑs)

    # if no dir make dir
    if !isdir("Nₑ")
        mkdir("Nₑ")
    end

    for Nₑ in Nₑs

        fn = findall(x->x==Nₑ,Nₑs)[1]

        p = parameters(
            ℓ = 1,
            κ = 0,
            Nₑ = Nₑ,
            fname = lpad(fn,3,"0"),
            fldr = "Nₑ/",
            T = 200
            )
        
        ens(p)

    end

end

function S(Ss)

    # if no dir make dir
    if !isdir("S")
        mkdir("S")
    end

    for S in Ss

        fn = findall(x->x==S,Ss)[1]

        p = parameters(
            ℓ = 1,
            κ = 0,
            S = S,
            fname = lpad(fn,3,"0"),
            fldr = "S/",
            T = 200
            )
        
        ens(p)

    end

end

function κ(κs)

    # if no dir make dir
    if !isdir("κ")
        mkdir("κ")
    end

    for κ in κs

        fn = findall(x->x==κ,κs)[1]

        p = parameters(
            ℓ = 1-κ,
            κ = κ,
            fname = lpad(fn,3,"0"),
            fldr = "κ/",
            T = 200
            )
        
        ens(p)

    end

end

Nₑs = trunc.(Int64,10 .^ (1.0:0.25:3.0))
Ss  = trunc.(Int64,10 .^ (1.0:0.25:3.0))
κs  = 0:0.125:1

Nₑ(Nₑs)
S(Ss)
κ(κs)

convert("Nₑ/",1:length(Nₑs))
convert("S/",1:length(Ss))
convert("κ/",1:length(κs))
