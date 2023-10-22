include("defs.jl");

print("\nPackages Loaded\n")

function SJₑ(n)

    # if no dir make dir
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/SJₑ")
        mkdir("dat/SJₑ")
    end

    Js = trunc.(Int64, 10 .^ rand(Uniform(1,5),n))
    Ss = trunc.(Int64, 10 .^ rand(Uniform(1,3),n))

    Π = Vector{parameters}()
    for r in 1:n
        π = parameters(
            ℓ = 1,
            κ = 0,
            Jₑᴹ= Js[r],
            Jₑᴱ= Js[r],
            S = Ss[r],
            K = 1,
            k = r,
            μ = 0,
            fldr="dat/SJₑ/"
            )
        push!(Π,π)
    end

    rp(Π)

    convert("dat/SJₑ/",n)

    cmmd1 = `echo "Has finished simulating across S & Jₑ"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))

end

n = 10^3
convert("dat/SJₑ/",n)
SJₑ(n)