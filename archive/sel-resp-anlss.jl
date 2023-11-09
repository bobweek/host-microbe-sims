using Gadfly, DataFrames, GLM, Colors, Cairo, Fontconfig, CSV, Compose

dat = CSV.read("sel_exprmnts.csv",DataFrame);

dfε0 = dat[dat.ε .==0.0,:]
dfℓ0 = dat[dat.ℓ .==0.0,:]

ols = lm(@formula(m² ~ ℓ), dfε0)
c₀, s₀ = coef(ols);
sᵣ = round(s₀, digits=3)
s = "slope = $sᵣ"
cᵣ = round(c₀, digits=3)
c = "y-int = $cᵣ"
p1 = plot(dfε0,x=:v,y=:a²,Geom.point,
    intercept=[c₀], slope=[s₀], Geom.abline(color="blue", style=:dash),
    Guide.xlabel("Vertical transmission (ℓ)"), Guide.ylabel("Microbial heritability (m²)"),
    Guide.title("Without Environmental Transmission, ε=0"),
    Guide.annotation(compose(context(), Compose.text(0.1,1.20,s), fill("blue"))),
    Guide.annotation(compose(context(), Compose.text(0.1,1.05,c), fill("blue"))))
p1 |> PNG("hubbellian/evowibo/ε0.png")

ols = lm(@formula(b ~ ℓ), dfε0)
c₀, s₀ = coef(ols);
sᵣ = round(s₀, digits=3)
s = "slope = $sᵣ"
cᵣ = round(c₀, digits=3)
c = "y-int = $cᵣ"
p1b = plot(dfε0,x=:ℓ,y=:b,Geom.point,
    intercept=[c₀], slope=[s₀], Geom.abline(color="blue", style=:dash),
    Guide.xlabel("Vertical transmission (ℓ)"), Guide.ylabel("Microbial bias (b)"),
    Guide.title("Without Environmental Transmission, ε=0"),
    Guide.annotation(compose(context(), Compose.text(0.1,-0.0020,s), fill("blue"))),
    Guide.annotation(compose(context(), Compose.text(0.1,-0.0025,c), fill("blue"))))
p1b |> PNG("hubbellian/evowibo/ε0b.png")

ols = lm(@formula(m² ~ ε), dfℓ0)
c₀, s₀ = coef(ols);
sᵣ = round(s₀, digits=3)
s = "slope = $sᵣ"
cᵣ = round(c₀, digits=3)
c = "y-int = $cᵣ"
p2 = plot(dfℓ0,x=:ε,y=:a²,Geom.point,
    intercept=[c₀], slope=[m²₀], Geom.abline(color="blue", style=:dash),
    Guide.xlabel("Environmental Transmission (ε)"), Guide.ylabel("Microbial heritability (m²)"),
    Guide.title("Without Lineal Transmission, ℓ=0"),
    Guide.annotation(compose(context(), Compose.text(0.1,0.80,s), fill("blue"))),
    Guide.annotation(compose(context(), Compose.text(0.1,0.70,c), fill("blue"))))
p2 |> PNG("hubbellian/evowibo/ℓ0.png")

ols = lm(@formula(b ~ ε), dfℓ0)
c₀, s₀ = coef(ols);    
sᵣ = round(s₀, digits=3)
s = "slope = $sᵣ"
cᵣ = round(c₀, digits=3)
c = "y-int = $cᵣ"
p2b = plot(dfℓ0,x=:ω,y=:b,Geom.point,
    intercept=[c₀], slope=[s₀], Geom.abline(color="blue", style=:dash),
    Guide.xlabel("Environmental transmission (ε)"), Guide.ylabel("Microbial bias (b)"),
    Guide.title("Without Lineal Transmission, ℓ=0"),
    Guide.annotation(compose(context(), Compose.text(0.1,-0.0020,s), fill("blue"))),
    Guide.annotation(compose(context(), Compose.text(0.1,-0.0025,c), fill("blue"))))
p2b |> PNG("hubbellian/evowibo/ℓ0b.png")

# interaction

dfℓε = dat[dat.ε .== dat.ℓ,:]

ols = lm(@formula(m² ~ ℓ + ℓ^2), dfℓε)
c₀, s₀, s²₀ = coef(ols);
s²ᵣ = round(s²₀, digits=3)
s² = "quadratic = $s²ᵣ"
sᵣ = round(s₀, digits=3)
s = "slope = $sᵣ"
cᵣ = round(c₀, digits=3)
c = "y-int = $cᵣ"
t = 0:0.01:1
p1 = Gadfly.plot(dfℓε,layer(x=:v,y=:a²,Geom.point),
    layer(x = t, y=c₀ .+ s₀*t .+ s²₀*t.^2, Geom.path),
    Guide.xlabel("Vertical transmission (ℓ)"), Guide.ylabel("Microbial heritability (m²)"),
    Guide.title("Equal Lineal and Environmental Transmission, ℓ=ε"),
    Guide.annotation(compose(context(), Compose.text(0.1,1.30,s²), fill("blue"))),
    Guide.annotation(compose(context(), Compose.text(0.1,1.10,s), fill("blue"))),
    Guide.annotation(compose(context(), Compose.text(0.1,0.9,c), fill("blue"))))
p1 |> PNG("hubbellian/evowibo/ℓε.png")

# full fit

ols = lm(@formula(m² ~ ε + ℓ + ε*ℓ + ε^2 + ℓ^2),dat)
d₀, εₛ, ℓₛ, ε²ₛ, ℓ²ₛ, εℓₛ = coef(ols)

fnrng, cb = 0.0:0.01:1.0, [colorant"black"]

D = DataFrame(x=repeat(fnrng, outer=length(fnrng)), y=repeat(fnrng, inner=length(fnrng)))

D.z = d₀ .+ εₛ*D.x .+ ℓₛ*D.y .+ ε²ₛ*D.x.^2 .+ ℓ²ₛ*D.y.^2 .+ εℓₛ*(D.x .* D.y)

plot(D, x=:x, y=:y,
  layer(z=:z, Geom.contour(levels=[-0.1:-0.1:-0.5;]), linestyle=[:dash], color=cb),
    layer(z=:z, Geom.contour(levels=[0.1:0.1:0.5;]),  color=cb))

plot(z=(x,y) -> d₀ .+ εₛ*x .+ ℓₛ*y .+ ε²ₛ*x.^2 .+ ℓ²ₛ*y.^2 .+ εℓₛ*(x .* y),
    xmin=[0], xmax=[1], ymin=[0], ymax=[1], Geom.contour,
    Guide.xlabel("Environmental transmission (ε)"), Guide.ylabel("Lineal Transmission (ℓ)"))

pW = plot(x -> exp(x), -2.0,2.0,
    Guide.xlabel("Phenotype (P)"), Guide.ylabel("Fitness (W)"), Guide.title("Fitness Function"))
pW |> PNG("hubbellian/evowibo/W.png")

using Plots; pyplot()
x=range(0,stop=1,length=100)
y=range(0,stop=1,length=100)
f(x,y) = d₀ .+ εₛ*x .+ ℓₛ*y .+ ε²ₛ*x.^2 .+ ℓ²ₛ*y.^2 .+ εℓₛ*(x .* y)
p = Plots.plot(x,y,f,st=:surface,camera=(-30,30),xlims=(0,1),ylims=(0,1),xlabel="ℓ", ylabel="ε",zlabel="m²",zguidefontrotation=90)
