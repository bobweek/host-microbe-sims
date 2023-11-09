include("defs.jl");

print("\nPackages Loaded\n")

function Nₑ(n)

    # if no dir make dir
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/posNₑ")
        mkdir("dat/posNₑ")
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
            fldr="dat/posNₑ/"
            )
        push!(Π,π)
    end

    rp(Π)

    convert("dat/posNₑ/",n)

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
    if !isdir("dat/posS")
        mkdir("dat/posS")
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
            fldr="dat/posS/"
            )
        push!(Π,π)
    end

    rp(Π)

    convert("dat/posS/",n)

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
    if !isdir("dat/posκ")
        mkdir("dat/posκ")
    end
    
    κ = rand(Uniform(),n)

    Π = Vector{parameters}()
    for r in 1:n
        π = parameters(
            ℓ=1-κ[r],
            κ=κ[r],
            K=1,
            k=r,
            fldr="dat/posκ/")
        push!(Π,π)
    end

    rp(Π)

    convert("dat/posκ/",n)

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
    if !isdir("dat/posκℓ")
        mkdir("dat/posκℓ")
    end
    
    κℓ = rand(Uniform(),(2,n))

    Π = Vector{parameters}()
    for r in 1:n
        π = parameters(
            ℓ=κℓ[1,r]*κℓ[2,r],
            κ=κℓ[1,r]*(1-κℓ[2,r]),
            K=1,
            k=r,
            fldr="dat/posκℓ/")
        push!(Π,π)
    end

    rp(Π)

    convert("dat/posκℓ/",n)

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

cmmd1 = `echo "Has finished positive selection simulations"`
cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
run(pipeline(cmmd1,cmmd2))