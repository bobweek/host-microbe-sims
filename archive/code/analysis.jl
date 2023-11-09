include("analysis-defs.jl");

dat = CSV.read("dat.csv",DataFrame);
met = CSV.read("met.csv",DataFrame);

#
# pooled hosts
#

pdat, pmet = host_pool(dat,met);

plot(DataFrame(x=values(pdat[190,:])),x="x",Geom.histogram)

# bray

pY = mds(pdat,"bray",2);

cmbpmet = deepcopy(pmet);
for i in 1:2
    colnm = "Coord$i"
    cmbpmet[!,colnm] = pY[:,i]
end

s = pmet.Loc .< 1
br = pmet.Age/maximum(pmet.Age)
clrs = HSV.(255*br,1,1);

p = plot(cmbpmet,x=:Coord1,y=:Coord2,Geom.point,color=clrs,shape=s);
p |> PNG("pcoa.bray.png")

# robust.aitchison

pY = mds(pdat,"robust.aitchison",2);

cmbpmet = deepcopy(pmet);
for i in 1:2
    colnm = "Coord$i"
    cmbpmet[!,colnm] = pY[:,i]
end

s = pmet.Loc .< 1
br = pmet.Age/maximum(pmet.Age)
clrs = HSV.(255*br,1,1);

p = plot(cmbpmet,x=:Coord1,y=:Coord2,Geom.point,color=clrs,shape=s);
p |> PNG("pcoa.raitchison.png")

# jaccard

pY = mds(pdat,"jaccard",2);

cmbpmet = deepcopy(pmet);
for i in 1:2
    colnm = "Coord$i"
    cmbpmet[!,colnm] = pY[:,i]
end

s = pmet.Loc .< 1
br = pmet.Age/maximum(pmet.Age)
clrs = HSV.(255*br,1,1);

p = plot(cmbpmet,x=:Coord1,y=:Coord2,Geom.point,color=clrs,shape=s);
p |> PNG("pcoa.jaccard.png")

# Dst = dst_mat(pdat,"bray");
# pcoa = fit(MDS,Dst; distances=true)
# evlspcoa = eigvals(pcoa::MDS)
# sum(evlspcoa[1:2])/sum(evlspcoa)

thing = fit(MDS,Matrix(pdat)'; distances=false)
pca = predict(thing)

# evlspca = eigvals(pca::MDS);
# sum(evlspca[1:2])/sum(evlspca)

cmbpmet = deepcopy(pmet);
for i in 1:2
    colnm = "Coord$i"
    cmbpmet[!,colnm] = pca[i,:]
end

s = pmet.Loc .< 1
br = pmet.Age/maximum(pmet.Age)
clrs = HSV.(255*br,1,1);

p = plot(cmbpmet,x=:Coord1,y=:Coord2,Geom.point,color=clrs,shape=s);
p |> PNG("pca.png")

#
# all the hosts all the time (plot averaged host locs in ord with oval representing variance [one stdev in ea dir])
#

# Y = mds(dat,"bray",2);

# cmbmet = deepcopy(met);

# for i in 1:2
#     colnm = "Coord$i"
#     cmbmet[!,colnm] = Y[:,i]
# end

# s = met.l .< 1
# br = met.a/maximum(met.a)
# clrs = HSV.(255*br,1,1);

# Gadfly.plot(cmbmet,x=:Coord1,y=:Coord2,Geom.point,color=clrs,shape=s)

