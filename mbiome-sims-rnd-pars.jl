include("defs.jl");

print("\nPackages Loaded\n")

function Nₑ(n)

    # if no dir make dir
    if !isdir("Nₑ")
        mkdir("Nₑ")
    end

    Nₑs = trunc.(Int64, 10 .^ rand(Uniform(1,3),n))

    Π = Vector{parameters}()
    for r in 1:n
        π = parameters(
            ℓ = 1,
            κ = 0,
            Nₑ=Nₑs[r],
            K=1,
            k=r,
            fldr="Nₑ/"
            )
        push!(Π,π)
    end

    rp(Π)

    convert("Nₑ/",n)

    print("\nNₑ done\n")

    cmmd1 = `echo "Has finished simulating across Nₑ"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))

end

function S(n)

    # if no dir make dir
    if !isdir("S")
        mkdir("S")
    end

    Ss = trunc.(Int64, 10 .^ rand(Uniform(1,3),n))

    Π = Vector{parameters}()
    for r in 1:n
        π = parameters(
            ℓ = 1,
            κ = 0,
            S=Ss[r],
            K=1,
            k=r,
            fldr="S/"
            )
        push!(Π,π)
    end

    rp(Π)

    convert("S/",n)

    print("\nS done\n")

    cmmd1 = `echo "Has finished simulating across S"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))
    
end

function κ(n)

    # if no dir make dir
    if !isdir("κ")
        mkdir("κ")
    end
    
    κ = rand(Uniform(),n)

    Π = Vector{parameters}()
    for r in 1:n
        π = parameters(
            ℓ=1-κ[r],
            κ=κ[r],
            K=1,
            k=r,
            fldr="κ/")
        push!(Π,π)
    end

    rp(Π)

    convert("κ/",n)

    print("\nκ done\n")

    cmmd1 = `echo "Has finished simulating across S"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))
        
end

function ℓκ(n)

    # if no dir make dir
    if !isdir("ℓκ")
        mkdir("ℓκ")
    end
    
    ℓκ = rand(Uniform(),(2,n))

    Π = Vector{parameters}()
    for r in 1:n
        π = parameters(
            ℓ=ℓκ[1,r]*ℓκ[2,r],
            κ=ℓκ[1,r]*(1-ℓκ[2,r]),
            K=1,
            k=r,
            fldr="ℓκ/")
        push!(Π,π)
    end

    rp(Π)

    convert("ℓκ/",n)

    print("\nℓκ done\n")

    cmmd1 = `echo "Has finished simulating across S"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))
    
end

n = 10^3

Nₑ(n)
S(n)
κ(n)
ℓκ(n)

cmmd1 = `echo "Has finished your simulations"`
cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
run(pipeline(cmmd1,cmmd2))