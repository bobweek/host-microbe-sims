#
# mds in r
#

@rimport stats
@rimport glue
@rimport ggplot2

pcoa = stats.cmdscale(dst,var"eig"=true)

# lists column one then columns two
pcoa[:points]

size(dat)

names(pcoa)

@rput pcoa

R"pcoa$points"

# can make own dist matrix to plug into cmdscale
# so might make kl-div fct for this

# maybe use canberra?? just to see??