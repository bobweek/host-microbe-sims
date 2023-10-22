include("defs.jl");

print("\nPackages Loaded\n")

function Nₑ(n)

    # if no dir make dir
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/negNₑ")
        mkdir("dat/negNₑ")
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
            s=-0.01,
            μ=0,
            fldr="dat/negNₑ/"
            )
        push!(Π,π)
    end

    rp(Π)

    convert("dat/negNₑ/",n)

    print("\nNₑ done\n")

    cmmd1 = `echo "Has finished simulating across Nₑ"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))

end

function S(n)

    # if no dir make dir
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/negS")
        mkdir("dat/negS")
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
            s=-0.01,
            μ=0,
            fldr="dat/negS/"
            )
        push!(Π,π)
    end

    rp(Π)

    convert("dat/negS/",n)

    print("\nS done\n")

    cmmd1 = `echo "Has finished simulating across S"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))
    
end

function κ(n)

    # if no dir make dir
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/negκ")
        mkdir("dat/negκ")
    end
    
    κ = rand(Uniform(),n)

    Π = Vector{parameters}()
    for r in 1:n
        π = parameters(
            ℓ=1-κ[r],
            κ=κ[r],
            K=1,
            k=r,
            s=-0.01,
            μ=0,
            fldr="dat/negκ/")
        push!(Π,π)
    end

    rp(Π)

    convert("dat/negκ/",n)

    print("\nκ done\n")

    cmmd1 = `echo "Has finished simulating across κ"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))
        
end

function κℓ(n)

    # if no dir make dir
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/negκℓ")
        mkdir("dat/negκℓ")
    end
    
    κℓ = rand(Uniform(),(2,n))

    Π = Vector{parameters}()
    for r in 1:n
        π = parameters(
            ℓ=κℓ[1,r]*κℓ[2,r],
            κ=κℓ[1,r]*(1-κℓ[2,r]),
            K=1,
            k=r,
            s=-0.01,
            μ=0,
            fldr="dat/negκℓ/")
        push!(Π,π)
    end

    rp(Π)

    convert("dat/negκℓ/",n)

    print("\nκℓ done\n")

    cmmd1 = `echo "Has finished simulating across κℓ"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))
    
end

n = 10^3

Nₑ(n)
S(n)
κ(n)
κℓ(n)

cmmd1 = `echo "Has finished negative selection simulations"`
cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
run(pipeline(cmmd1,cmmd2))