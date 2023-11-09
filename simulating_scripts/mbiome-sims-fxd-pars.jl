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
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/Nₑ")
        mkdir("dat/Nₑ")
    end

    for Nₑ in Nₑs

        fn = findall(x->x==Nₑ,Nₑs)[1]

        p = parameters(
            ℓ = 1,
            κ = 0,
            Nₑ = Nₑ,
            fname = lpad(fn,3,"0"),
            fldr = "dat/Nₑ/",
            T = 200
            )
        
        ens(p)

    end

    convert("dat/Nₑ/",1:length(Nₑs))

    cmmd1 = `echo "Has finished simulating across Nₑ"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))

end

function S(Ss)

    # if no dir make dir
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/S")
        mkdir("dat/S")
    end    

    for S in Ss

        fn = findall(x->x==S,Ss)[1]

        p = parameters(
            ℓ = 1,
            κ = 0,
            S = S,
            fname = lpad(fn,3,"0"),
            fldr = "dat/S/",
            T = 200
            )
        
        ens(p)

    end

    convert("dat/S/",1:length(Ss))

    cmmd1 = `echo "Has finished simulating across S"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))

end

function κ(κs)

    # if no dir make dir
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/κ")
        mkdir("dat/κ")
    end    

    for κ in κs

        fn = findall(x->x==κ,κs)[1]

        p = parameters(
            ℓ = 1-κ,
            κ = κ,
            fname = lpad(fn,3,"0"),
            fldr = "dat/κ/",
            T = 200
            )
        
        ens(p)

    end

    convert("dat/κ/",1:length(κs))

    cmmd1 = `echo "Has finished simulating across κ"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))

end

function Jₑ(Jₑs)

    # if no dir make dir
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/Jₑ")
        mkdir("dat/Jₑ")
    end

    for Jₑ in Jₑs

        fn = findall(x->x==Jₑ,Jₑs)[1]

        p = parameters(
            ℓ = 1,
            κ = 0,
            Jₑᴹ = Jₑ,
            Jₑᴱ = Jₑ,
            fname = lpad(fn,3,"0"),
            fldr = "dat/Jₑ/",
            T = 200
            )
        
        ens(p)

    end

    convert("dat/Jₑ/",1:length(Jₑs))

    cmmd1 = `echo "Has finished simulating across Jₑ"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))

end

Nₑs = trunc.(Int64,10 .^ (1.0:1.0:3.0))
κs  = 0:0.5:1

Nₑ(Nₑs)
κ(κs)

# Ss  = trunc.(Int64,10 .^ (1.0:1.0:3.0))
# S(Ss)

# Jₑs = trunc.(Int64,10 .^ (2.0:2.0:6.0))
# Jₑ(Jₑs)

cmmd1 = `echo "Has finished all simulations"`
cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
run(pipeline(cmmd1,cmmd2))
