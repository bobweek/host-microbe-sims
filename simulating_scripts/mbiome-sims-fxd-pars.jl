#
#
# simulate ensemble of replicates across fixed background parameters
#
#

include("defs.jl");

print("\nPackages Loaded\n")

Nₑs = trunc.(Int64,10 .^ (1.0:1.0:3.0))
ensNₑ(Nₑs)

κs  = 0:0.5:1
ensκ(κs)

Ss  = trunc.(Int64,10 .^ (1.0:1.0:3.0))
ensS(Ss)

Jₑs = trunc.(Int64,10 .^ (2.0:2.0:6.0))
ensJₑ(Jₑs)

cmmd1 = `echo "Has finished fixed parameter simulations"`
cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
run(pipeline(cmmd1,cmmd2))