#
#
# simulate across random background parameters
#
#

include("defs.jl");

print("\nPackages Loaded\n")

n = 10^3

rpNₑ(n)

rpS(n)

rpκ(n)

rpκℓ(n)

cmmd1 = `echo "Has finished random parameter simulations"`
cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
run(pipeline(cmmd1,cmmd2))