#####################################
#									#
# definitions for					#
#	analytical model and			#
#	microbiome explicit simulations #
#									#
#####################################

using Parameters, LinearAlgebra, Random, Distributions, DataFrames, CSV, Documenter, DocumenterDiagrams, DocumenterTools

# 
# functions used in both analytical and simulation models
#

#region

"""
returns `ε` as a function of `ℓ` and `κ`,
	where `ε` is probability that a microbe 
	is environmentally acquired

ε = 1 - κ / (1 - ℓ)
"""
function ε₀(ℓ,κ)

	if ℓ < 1
		return 1 - κ / (1 - ℓ)
	else
		return 1
	end

end

"""
checks that inheritance parameters are sensible
- ℓ,κ ≥ 0
- ℓ+κ ≤ 1
"""
function inh_check(ℓ,κ)

	if (ℓ < 0) | (κ < 0)
		print("\nError: Negative Inheritance Parameters\n")
		exit()
	end
	if (ℓ+κ)>1
		print("\nError: Inheritance Parameters Sum > 1\n")
		exit()
	end

	return
end


#endregion

#
# analytical model
#

#region

"""
parameters and state variables for analytical model

model parameters
- `ℓ` = lineal inheritance
- `κ` = collective inheritance
- `ε` = probability of environmental ancestry
- `E` = variance of noise
- `G` = additive genetic variance
- `M` = additive microbial variance
- `P` = phenotypic variance
- `β` = selection gradient
	- equals selection strength `s` under assumptions in main text
- `T` = number of host gens to iterate
- `k` = parameter combination index

state variables
- `z̄` = mean trait value
- `ḡ` = mean additive genetic value
- `m̄` = mean additive microbial value
- `ξ` = environmental microbial value
- `t` = host generation
"""
@with_kw mutable struct asystem

	# model parameters
	ℓ::Float64	= (3-√5)/2	# lineal inheritance
	κ::Float64	= (3-√5)/2	# collective inheritance
	ε::Float64	= ε₀(ℓ,κ)	# prb of env ancestry
	E::Float64	= 0.0		# variance of noise
	G::Float64	= 1.0		# additive genetic variance
	M::Float64	= 1.0		# additive microbial variance
	P::Float64	= G+M+E		# phenotypic variance
	β::Float64	= 1.0		# selection gradient
	T::Int64	= 50		# number host gens to iterate
	k::Int64	= 0			# parameter combination

	# state variables
	ḡ::Float64 = 0.0		 # mean additive genetic value
	m̄::Float64 = 0.0		 # mean additive microbial value
	z̄::Float64 = ḡ+m̄		   # mean trait value
	ξ::Float64 = 0.0		# environmental microbial value
	t::Int64   = 0			# host generation

end

"""
iterates the analytical model following the equations
- Δḡ = Gβ
- Δm̄ = (κ+ℓ)Mβ + (1-κ-ℓ)(ξ-m̄)
- Δξ = (1-ε)(Mβ+m̄-ξ)
"""
function Δ(as::asystem)
	@unpack_asystem as

	# TODO: move note to md
	# note:
	# ξ' = εξ + (1-ε)m̄*
	# m̄*-m̄ = Mβ
	# ∴ ξ' = εξ + (1-ε)(Mβ+m̄)
	# ∴ Δξ = (1-ε)(Mβ + m̄ - ξ)
	#
	#		= κ (Mβ + m̄ - ξ) / (1 - ℓ)
	#
	#		so κ and ℓ cannot be combined into v...

	Δξm̄ = ξ - m̄

	ḡ += G*β
	m̄ += (κ+ℓ)*M*β + (1-κ-ℓ)*Δξm̄
	ξ += (1-ε)*(M*β-Δξm̄)

	z̄ = ḡ + m̄

	t += 1

	@pack_asystem! as
	return as
end

"""
writes analytical model parameters and state variables to csv
- if `k,t=0`, creates new csv
	- otherwise appends to csv
"""
function asave(as::asystem)
	@unpack_asystem as

	dat = DataFrame(
		z = z̄,
		g = ḡ,
		m = m̄,
		ξ = ξ,
		t = t,
		ℓ = ℓ,
		κ = κ,
		ε = ε,
		E = E,
		G = G,
		P = P,
		β = β,
		T = T,
		k = k
	)

	# if no dir make dir
	if !isdir("dat")
		mkdir("dat")
	end
	
	# if first entry, make new file
	if (k == 0) & (t == 0)
		CSV.write("dat/a_dat.csv", dat)
	else # append to current file
		CSV.write("dat/a_dat.csv", dat, append=true)
	end

	return
end

"""
runs the analytical model for `T` host generations
- `T` is specified in type `asystem`
- saves time series using function `asave`
"""
function arun(as::asystem)

	# make sure ℓ + κ ≤ 1 and ℓ,κ > 0
	inh_check(as.ℓ,as.κ)

	asave(as)
	for t in 1:as.T
		Δ(as)
		asave(as)
	end

	return
end

#endregion

#
# simulation model
#

#region

"""
variance of per capita additive effects of microbial species
- designes to make variance of microbiome dynamics comparable to variance of mutation
"""
function σ²(S,Jₑ,L)
	return S/(Jₑ^2*(S-1)*((1-(1-1/Jₑ)^L)))
end

"""
parameters for simulation model

system parameters
- `ℓ`  = lineal inheritance
- `κ`  = collective inheritance
- `ε`  = probability of environmental ancestry
- `L`  = microbiome dynamics per host generation
- `S`  = number of microbe species
- `Nₑ`= number of hosts
- `s`  = host selection strength
- `μ`  = variance of mutation
- `Jₑᴱ` = number of microbes in environment
- `Jₑᴹ` = number of microbes in a host
- `aqsz` = size of sample for acquired microbiomes
- `shsz` = size of sample for shedded microbiomes

simulation parameters
- `T` = number of host generations to run model
- `K` = number of replicates
- `σ` = std-dev of per capita microbial effects
- `B` = host gens burnin period
- `fldr` = folder to save data in
- `fname` = file name for saving data and parameters

simulation variable
- `k` = replicate or parameter combination index
"""
@with_kw mutable struct parameters

	# system parameters
	ℓ::Float64	= 1/3		# lineal inheritance
	κ::Float64	= 1/3		# collective inheritance
	ε::Float64	= ε₀(ℓ,κ)	# prb of env ancestry
	L::Int64	= 30		# microbiome dynamics per host generation
	S::Int64	= 100		# number of microbe species
	Nₑ::Int64	= 200		# number of hosts
	s::Float64	= 0.01		# host selection strength
	μ::Float64	= 1.0		# variance of mutation
	Jₑᴱ::Int64	= 10^6		# number of microbes in environment
	Jₑᴹ::Int64	= 10^6		# number of microbes in a host
	aqsz::Int64 = 100		# size of sample for acquired microbiomes
	shsz::Int64 = 100		# size of sample for shedded microbiomes

	# simulation parameters
	T::Int = 50					# number of host generations to run model
	K::Int = 100				# number of replicates
	σ::Float64 = √σ²(S,Jₑᴹ,L)	# std-dev of per capita microbial effects
	B::Int64 = 0				# host gens burnin period
	fname::String = "000"		# file name for saving data/parameters
	fldr::String = "tst/"		# folder to save data in
	
	# simulation variable
	k::Int = 1					# current replicate

end

"""
attributes of a host individual
- `α` = vector of per capita effect on host trait across microbe species
- `f` = vector of microbe relatlive abundances in host microbiome
- `p` = index of host parent in previous generation (zero implies first generation)
- `g` = additive genetic value
- `m` = additive microbial value
- `c` = additive microbial offset (s.t. m̄₀ = ḡ₀ = 0)
- `z` = trait value
- `w` = realized relative fitness
- `W` = realized absolute fitness
"""
@with_kw mutable struct host

	α::Vector{Float64}	# per capita microbe effect on host trait
	f::Vector{Float64}	# microbe rel abundances in host microbiome
	p::Int64   = 0		# index of host parent in previous generation (zero => first gen)
	g::Float64 = 0.0	# additive genetic value
	m::Float64 = 0.0	# additive microbial value
	c::Float64 = 0.0	# additive microbial offset (s.t. m̄₀ = ḡ₀ = 0)
	z::Float64 = 0.0	# trait value
	w::Float64 = 1.0	# realized relative fitness
	W::Int64   = 1		# realized absolute fitness

end

"""
attributes of simulation model
- `fᴱ`= microbe relatlive abundances in environment
- `c` = environmental microbial offset (s.t. ξ₀ = 0)
- `H` = hosts
- `t` = host generation
"""
@with_kw mutable struct system

	fᴱ::Vector{Float64}	# microbe rel abundances in environment
	c::Float64 = 0.0	# env microbial offset (s.t. ξ₀=0)
	H::Vector{host}		# hosts
	t::Int64 = 0		# host generation

end

"""
single iteration of stochastic microbiome dynamics,
	where `Jₑ` is microbiome size and `f` is initial vector of
	microbe relative abundances
"""
function hubbell(Jₑ,f)
	f = abs.(f)
	normalize!(f,1)
	f = normalize(rand(Multinomial(Jₑ,f)),1)
	return f
end

"""
computes expected relative fitnesses of host individuals
- expected absolute fitness given by exp(s*z)
"""
function fitness(sys::system,π::parameters)
	@unpack_system sys
	@unpack_parameters π

	W = Vector{Float64}()
	for h in H
		push!(W, exp(s*h.z))
	end
	w = normalize(W,1)

	return w
end

"""
performs reproduction of host individuals, which includes:
- genetic mutation
- microbiome transmission
- host parents stochastically shedding microbes
- host offspring stochastically sampling microbes
"""
function hrep(sys::system,π::parameters)
	@unpack_system sys
	@unpack_parameters π

	# form pooled microbiome of selected hosts
	fᴾ = zeros(S)
	for h in H
		fᴾ .+= h.W .* h.f
	end
	normalize!(fᴾ,1)
	
	# draw genetic mutations
	μˢ = rand(Normal(0,√μ),Nₑ)

	# form host offspring population
	Hₒ = Vector{host}()
	i = 1
	k = 1
	for h in H
		
		# offspring of host h
		for j in 1:h.W

			# sample from environmental microbiome
			fᴱᴾ = normalize(rand(Multinomial(aqsz,fᴱ)),1)
			
			# sample from pooled parental microbiome
			fᴾᴾ = normalize(rand(Multinomial(aqsz,fᴾ)),1)

			# initial offspring microbiome
			f = ((1-ℓ-κ) .* fᴱᴾ) .+ 
				(ℓ .* (h.f)) .+ (κ .* fᴾᴾ)

			push!(Hₒ, host(
				α = deepcopy(h.α),
				f = deepcopy(f),
				p = deepcopy(i),
				c = deepcopy(h.c),
				g = h.g + μˢ[k]
				))
			k += 1
		end

		i += 1

	end
	H = Hₒ

	# environment after hosts shed microbes
	fᴾᴾ = normalize(rand(Multinomial(shsz,fᴾ)),1)
	fᴱ = (ε .* fᴱ) .+ ((1-ε) .* fᴾᴾ)

	@pack_system! sys
	return sys
end

"""
performs one iteration of simulation model
- equals one host generation
- includes:
	- host reproduction
	- host development
	- computation of host traits
	- computation of realized fitnesses
"""
function hgen(sys::system,π::parameters)
	
	# host reproduction
	sys = hrep(sys,π)	

	@unpack_system sys
	@unpack_parameters π

	# host development
	for i in 1:L

		# environmental microbiome dynamic
		fᴱ = hubbell(Jₑᴱ,fᴱ)

		# host microbiome dynamic
		for h in H
			h.f = hubbell(Jₑᴹ,h.f)
		end
		
	end
	
	# compute host traits
	for h in H
		# scaling by Jₑᴹ makes α a per capita effect
		h.m = Jₑᴹ*dot(h.α,h.f) - h.c
		h.z = h.g + h.m
	end

	# compute host expected fitnesses
	w = fitness(sys,π)

	# compute host realized fitnesses
	Nₒ = rand(Multinomial(Nₑ,w))
	for i in 1:Nₑ

		H[i].W = Nₒ[i]
		H[i].w = Nₒ[i]/Nₑ

		# note:
		#	because sum(Nₒ)/Nₑ = 1
		#	Nₒ/Nₑ are the realized rel fitnesses of each adult host

	end

	@pack_system! sys
	return sys
end

"""
generates initial conditions
	using `B` iterations of `hgen`,
	where `B` is a property of `parameters`
"""
function burnin(π::parameters)
	@unpack_parameters π

	# per capita microbial effects
	α = rand(Normal(0,σ),S)

	# env mbiome
	fᴱ = fill(1/S,S)

	# host vector
	f = fill(1/S,S)
	c = Jₑᴹ*dot(α,f)
	h = host(α=α, f=f, c=c)
	H = Vector{host}()
	for k in 1:Nₑ
		push!(H,deepcopy(h))
	end

	# make system
	sys = system(fᴱ=fᴱ, c=c, H=H)

	# do burnin
	for i in 1:B
		hgen(sys,π)
	end

	@pack_parameters! π
	return sys, π
end

"""
creates dataframe where each row corresponds to a host individual, and columns correspond to the following:
- h    = host individual index
- p    = host parent index
- g    = additive genetic value
- m    = additive microbial value
- z    = trait value
- W    = realized absolute fitness
- t    = host generation
- k    = simulation replicate number
- type = indicator of whether environment or host
- f1   = realative abundance of microbe spp 1
- ⋮    = ⋮
- fS   = realative abundance of microbe spp S
- a1   = per capita additive effect of microbe spp 1
- ⋮    = ⋮
- aS   = per capita additive effect of microbe spp S
"""
function ind_df(sys::system,π::parameters)
	@unpack_system sys
	@unpack_parameters π

	hi    = Vector{Int64}()
	pi    = Vector{Int64}()
	gs    = Vector{Float64}()
	ms    = Vector{Float64}()
	zs    = Vector{Float64}()
	Ws    = Vector{Float64}()
	αs    = Array{Float64}(undef,0,S)
	fs    = Array{Float64}(undef,0,S)
	types = Vector{String}()

	# hi    = vector of host individual indices
	# pi    = vector of host parent indices
	# gs    = vector of additive genetic values
	# ms    = vector of additive microbial values
	# zs    = vector of trait values
	# Ws    = vector of realized absolute fitnesses
	# αs    = vector of per capita additive effects
	# fs    = vector of microbe realative abundances
	# types = vector of whether environment or host
	
	i = 0
	push!(hi,i)
	push!(pi,0)
	push!(gs,0)
	push!(ms,Jₑᴱ*dot(fᴱ,H[1].α)-c)
	push!(zs,0)
	push!(Ws,0)
	αs = [αs;H[1].α']
	fs = [fs;fᴱ']	
	push!(types,"env")

	for h in H
		i += 1
		push!(hi,i)
		push!(pi,h.p)
		push!(gs,h.g)
		push!(ms,h.m)
		push!(zs,h.z)
		push!(Ws,h.W)
		αs = [αs;h.α']
		fs = [fs;h.f']
		push!(types,"hst")
	end

	idat = DataFrame(
		h=hi,
		p=pi,
		g=gs,
		m=ms,
		z=zs,
		W=Ws,
		t=fill(t,Nₑ+1),
		k=fill(k,Nₑ+1),
		type=types
	)

	for i in 1:S
		colname = "f$i"
		idat[!,colname] = fs[:,i]
	end
	
	for i in 1:S
		colname = "a$i"
		idat[!,colname] = αs[:,i]
	end

	return idat

end

"""
creates dataframe with single row summarizing a population with columns corresponding to the following:
- z = mean trait
- g = mean additive genetic value
- m = mean additive microbial value
- ξ = environmental microbial value
- P = phenotypic variance
- G = additive genetic variance
- M = additive microbial variance
- corr = correlation of genetic and microbial values
- β = realized selection gradient
- t = host generation
- k = replicate or parameter combination index
"""
function pop_df(sys::system,π::parameters)

	idat = ind_df(sys,π)

	envdat = idat[idat.type .== "env", :]
	hstdat = idat[idat.type .== "hst", :]

	return DataFrame(
		z = mean(hstdat.z),
		g = mean(hstdat.g),
		m = mean(hstdat.m),
		ξ = envdat.m[1],
		P = var(hstdat.z),
		G = var(hstdat.g),
		M = var(hstdat.m),
		corr = cov(hstdat.g,hstdat.m)/sqrt(G*M),
		β = cov(hstdat.W,hstdat.z)/P,
		t = hstdat.t[1],
		k = hstdat.k[1]
	)

end

"""
old method to save parameters.
not currently being used.
"""
function savepar(π::parameters)
	@unpack_parameters π
	par = DataFrame(
		sig = σ,
		ell = ℓ,
		kap = κ,
		eps = ε,
		L   = L,
		S   = S,
		Ne  = Nₑ,
		s   = s,
		mu  = μ,
		JeE = Jₑᴱ,
		JeM = Jₑᴹ,
		K   = K,
		T   = T,
		B   = B
	)
	CSV.write(fldr*"par_"*fname*".csv", par)
end

"""
appends population data to a dataframe with columns defined in pop_df
"""
function popcat(sys::system,π::parameters,popdf::DataFrame)
	pdf = pop_df(sys,π)
	return vcat(popdf,pdf)
end

"""
runs a simulation

- includes burnin
- results in a population dataframe
	- with `T` rows
	- columns defined by function `pop_df`
"""
function sim(π::parameters)

	inh_check(π.ℓ,π.κ)

	sys, π = burnin(π)
	popdf = pop_df(sys,π)

	for t in 1:π.T
		hgen(sys,π)
		sys.t += 1
		popdf = popcat(sys,π,popdf)
	end

	return popdf
end

"""
converts `parameters` to a single-row dataframe
"""
function π2df(π::parameters)
	@unpack_parameters π
	df = DataFrame(
		sig = σ,
		ell = ℓ,
		kap = κ,
		eps = ε,
		L   = L,
		S   = S,
		Ne  = Nₑ,
		s   = s,
		mu  = μ,
		JeE = Jₑᴱ,
		JeM = Jₑᴹ,
		K   = K,
		T   = T,
		B   = B,
		k   = k
	)
	return df
end

"""
runs ensemble of replicate simulations
- designed for iterating across fixed sets of parameters
- writes parameters to csv
- concatenates simulated population dataframes into one large dataframe, and writes to csv
	- location and name of csv determined by `fldr` and `fname` properties of `parameters`
- uses multithreading to run replicates in parallel
"""
function ens(π::parameters)

	# build parameters dataframe across replicates
	rps = DataFrame()
	for k in 1:π.K
		rps = vcat(rps,π2df(π))
	end

	# write parameter data to csv
	CSV.write(π.fldr*"par_"*π.fname*".csv", rps)

	# build vector of population dataframes across replicates
	repdfs = Vector{DataFrame}(undef, π.K)
	Threads.@threads for k in 1:π.K
		local πₚ = deepcopy(π)
		πₚ.k = k
		repdfs[k] = sim(πₚ)
	end

	# convert vector of dataframes to dataframe
	pdat = repdfs[1]
	for k in 2:π.K
		pdat = vcat(pdat,repdfs[k])
	end

	# write population-level data to csv
	CSV.write(π.fldr*"pop_dat_"*π.fname*".csv", pdat)

end

"""
runs collection of simulations across randomly drawn parameters
- `rp` stands for _random parameters_
- takes vector of `parameters` type with randomly drawn entries
- concatenates simulated population dataframes into one large dataframe, and writes to csv
	- location of csv determined by `fldr` property of `parameters`
- writes `parameters` to csv
- uses multithreading to run simulations in parallel

# Example
```julia-repl
n = 10
Nₑs = trunc.(Int64, 10 .^ rand(Uniform(1,3),n))
Π = Vector{parameters}()
for r in 1:n
	π = parameters(
		Nₑ=Nₑs[r],
		K=1,
		k=r,
		fldr="dat/posNₑ/"
		)
	push!(Π,π)
end
rp(Π)
```
"""
function rp(Π::Vector{parameters})

	# build dataframe of parameters across reps
	rps = DataFrame()
	for π in Π
		rps = vcat(rps,π2df(π))
	end

	# write parameter data to csv
	CSV.write(Π[1].fldr*"rndpars.csv", rps)

	# build vector of population dataframes across replicates
	pdfs = Vector{DataFrame}(undef, length(Π))
	Threads.@threads for π in Π
		local k = π.k
		pdfs[k] = sim(π)
	end

	# convert vector of dataframes to dataframe
	pdat = DataFrame()
	for π in Π
		pdat = vcat(pdat,pdfs[π.k])
	end

	# write population-level data to csv
	CSV.write(Π[1].fldr*"rndpar_dat.csv", pdat)

end

#endregion

#
# qunatifying analytical models ability to predict simulation model
#

#region

"""
summarizes state of analytical model
- used in predicting state of simulation model

state variables
- ḡ = mean additive genetic value
- m̄ = mean additive microbial value
- ξ = environmental microbial value
- z̄ = mean trait value
"""
@with_kw mutable struct astate
	ḡ::Float64 = 0.0	 # mean additive genetic value
	m̄::Float64 = 0.0	 # mean additive microbial value
	ξ::Float64 = 0.0	# environmental microbial value
	z̄::Float64 = 0.0	 # mean trait value
end

"""
reads parameter and population data from csv's
- designed to read output of `ens`
"""
function read(fldr::String,fname::String)

	pdat = CSV.read(fldr*"pop_dat_"*fname*".csv", DataFrame)
	par = CSV.read(fldr*"par_"*fname*".csv", DataFrame)

	return pdat, par
	
end

"""
reads parameter and population data from csv's
- designed to read output of `rp`
"""
function read(fldr::String)

	pdat = CSV.read(fldr*"rndpar_dat.csv", DataFrame)
	par = CSV.read(fldr*"rndpars.csv", DataFrame)

	return pdat, par
	
end

"""
predicts simulated state using analytical model, compares using cohen's _d_, and outputs dataframe with relevant parameters and data
- used in the `convert` functions
- `pdat` and `par` are typically obtained using the `read` function
- `k` specifies replicate or parameter combination index
- `t` specifies the current host generation
	- must satisfy `t < T`
	- analytical prediction will be compared against simulated state at `t+1`

columns of outputted dataframe
- `κ`  = collective inheritance
- `κℓ`= total microbial inheritance (κ + ℓ)
- `k`  = replicate or parameter combination number
- `β` = realized selection gradient
- `P` = phenotypic variance
- `G` = additive genetic variance
- `M` = additive microbial variance
- `Dza` = analytically predicted change in mean trait
- `Dzs` = simulated change in mean trait
- `Dga` = analytically predicted change in mean genetic value
- `Dgs` = simulated change in mean genetic value
- `Dma` = analytically predicted change in mean microbial value
- `Dms` = simulated change in mean microbial value
- `dz` = cohen's _d_ for trait value: (z̄ₐ - z̄ₛ)/√P
	- z̄ₐ = analytical prediction
	- z̄ₛ = simulated value
- `dg` = cohen's _d_ for genetic value: (ḡₐ - ḡₛ)/√G
	- ḡₐ = analytical prediction
	- ḡₛ = simulated value
- `dm` = cohen's _d_ for microbial value: (m̄ₐ - m̄ₛ)/√M
	- m̄ₐ = analytical prediction
	- m̄ₛ = simulated value
- `corr` = corr(g,m)
"""
function cohen(pdat::DataFrame,par::DataFrame,t::Int64,k::Int64)

	subdat = pdat[pdat.k .== k, :]

	# current simulated state
	crtdat = subdat[subdat.t .== t, :]

	if length(crtdat.z) == 0
		print("Error: attempted to read at t > T")
		exit()
	end

	# next simulated state
	nxtdat = subdat[subdat.t .== (t+1), :]

	if length(nxtdat.z) == 0
		print("Error: attempted to read at t+1 > T")
		exit()
	end

	# current state
	cs = asystem(

		# model parameters
		ℓ = par.ell[k],		# lineal inheritance
		κ = par.kap[k],		# collective inheritance
		ε = par.eps[k],		# prb of env ancestry
		E = 0.0,			# variance of noise
		G = crtdat.G[1],	# additive genetic variance
		M = crtdat.M[1],	# additive microbial variance
		P = crtdat.P[1],	# phenotypic variance
		β = crtdat.β[1],	# selection gradient (= s under model assumptions in main text)

		# state variables
		ḡ = crtdat.g[1],	 # mean additive genetic value
		m̄ = crtdat.m[1],	 # mean additive microbial value
		z̄ = crtdat.z[1],	 # mean trait value
		ξ = crtdat.ξ[1]		# env additive microbial value

	)

	# simulated state
	ss = astate(
		ḡ = nxtdat.g[1],	 # mean additive genetic value
		m̄ = nxtdat.m[1],	 # mean additive microbial value
		z̄ = nxtdat.z[1],	 # mean trait value
		ξ = nxtdat.ξ[1]		# env additive microbial value
	)

	as  = Δ(cs)
	
	ρ = DataFrame(
		κ  = cs.κ,
		κℓ = cs.κ + cs.ℓ,
		k  = k,
		β = cs.β,
		P = cs.P,
		G = cs.G,
		M = cs.M,
		Dza = as.z̄ - cs.z̄,
		Dzs = ss.z̄ - cs.z̄,
		Dga = as.ḡ - cs.ḡ,
		Dgs = ss.ḡ - cs.ḡ,
		Dma = as.m̄ - cs.m̄,
		Dms = ss.m̄ - cs.m̄,
		dz = (as.z̄ - ss.z̄)/sqrt(cs.P),
		dg = (as.ḡ - ss.ḡ)/sqrt(cs.G),
		dm = (as.m̄ - ss.m̄)/sqrt(cs.M),
		corr = crtdat.corr[1]
	)

	return ρ

end

"""
uses output of `ens` to quantify analytical models ability to predict simulation model
- writes results in csv file with location determined by `fldr` property of `parameters`
- adds three additional columns to dataframe made by `cohen` (`S`,`Ne`,`Je`)

columns of outputted csv
- `S` = number of microbe species
- `Ne`= host population size
- `Je`= number of microbes in microbiome
- `κ`  = collective inheritance
- `κℓ`= total microbial inheritance (κ + ℓ)
- `k`  = replicate number
- `β` = realized selection gradient
- `P` = phenotypic variance
- `G` = additive genetic variance
- `M` = additive microbial variance
- `Dza` = analytically predicted change in mean trait
- `Dzs` = simulated change in mean trait
- `Dga` = analytically predicted change in mean genetic value
- `Dgs` = simulated change in mean genetic value
- `Dma` = analytically predicted change in mean microbial value
- `Dms` = simulated change in mean microbial value
- `dz` = cohen's _d_ for trait value: (z̄ₐ - z̄ₛ)/√P
	- z̄ₐ = analytical prediction
	- z̄ₛ = simulated value
- `dg` = cohen's _d_ for genetic value: (ḡₐ - ḡₛ)/√G
	- ḡₐ = analytical prediction
	- ḡₛ = simulated value
- `dm` = cohen's _d_ for microbial value: (m̄ₐ - m̄ₛ)/√M
	- m̄ₐ = analytical prediction
	- m̄ₛ = simulated value
- `corr` = corr(g,m)
"""
function convert(fldr::String,pcombs::UnitRange{Int64})

	fnldat = DataFrame()

	for i in pcombs

		fname = lpad(i,3,"0")

		pdat, par = read(fldr,fname)

		for t in 1:(par.T[1]-1)
			for k in 1:par.K[1]
				ρ = cohen(pdat,par,t,k)
				tmpdat = hcat(DataFrame(
					S=par.S[k],
					Ne=par.Ne[k],
					Je=par.JeM[k]),ρ)
				fnldat = vcat(fnldat,tmpdat)
			end
		end
	end

	CSV.write(fldr*"fnl_fxdpar_dat.csv", fnldat)

end

"""
uses output of `rp` to quantify analytical models ability to predict simulation model
- writes results in csv file with location determined by `fldr` property of `parameters`
- adds three additional columns to dataframe made by `cohen` (`S`,`Ne`,`Je`)

columns of outputted csv
- `S` = number of microbe species
- `Ne`= host population size
- `Je`= number of microbes in microbiome
- `κ`  = collective inheritance
- `κℓ`= total microbial inheritance (κ + ℓ)
- `k`  = parameter combination number
- `β` = realized selection gradient
- `P` = phenotypic variance
- `G` = additive genetic variance
- `M` = additive microbial variance
- `Dza` = analytically predicted change in mean trait
- `Dzs` = simulated change in mean trait
- `Dga` = analytically predicted change in mean genetic value
- `Dgs` = simulated change in mean genetic value
- `Dma` = analytically predicted change in mean microbial value
- `Dms` = simulated change in mean microbial value
- `dz` = cohen's _d_ for trait value: (z̄ₐ - z̄ₛ)/√P
	- z̄ₐ = analytical prediction
	- z̄ₛ = simulated value
- `dg` = cohen's _d_ for genetic value: (ḡₐ - ḡₛ)/√G
	- ḡₐ = analytical prediction
	- ḡₛ = simulated value
- `dm` = cohen's _d_ for microbial value: (m̄ₐ - m̄ₛ)/√M
	- m̄ₐ = analytical prediction
	- m̄ₛ = simulated value
- `corr` = corr(g,m)
"""
function convert(fldr::String,n::Int64)

	# n is number replicates

	fnldat = DataFrame()

	pdat, par = read(fldr)

	tmpdat = Vector{DataFrame}(undef, n)
	for t in 1:(par.T[1]-1)
		Threads.@threads for k in 1:n
			ρ = cohen(pdat,par,t,k)
			tmpdat[k] = hcat(DataFrame(
				S=par.S[k],
				Ne=par.Ne[k],
				Je=par.JeM[k]),ρ)
		end
		for k in 1:n
			fnldat = vcat(fnldat,tmpdat[k])
		end
	end

	CSV.write(fldr*"fnl_rndpar_dat.csv", fnldat)

end

#endregion

#
# functions that simulate and analyze results across randomly drawn parameters
#

#region

"""
simulate and analyze results across randomly drawn host population size `Nₑ`
- `n` is number of random draws
- if no local folder `dat`, `rpNₑ` creates this folder
	- if no subfolder `Nₑ`, `rpNₑ` creates this folder
- outputs results to csv file in `/dat/Nₑ`

columns of outputted csv
- `S` = number of microbe species
- `Ne`= host population size
- `Je`= number of microbes in microbiome
- `κ`  = collective inheritance
- `κℓ`= total microbial inheritance (κ + ℓ)
- `k`  = replicate number
- `β` = realized selection gradient
- `P` = phenotypic variance
- `G` = additive genetic variance
- `M` = additive microbial variance
- `Dza` = analytically predicted change in mean trait
- `Dzs` = simulated change in mean trait
- `Dga` = analytically predicted change in mean genetic value
- `Dgs` = simulated change in mean genetic value
- `Dma` = analytically predicted change in mean microbial value
- `Dms` = simulated change in mean microbial value
- `dz` = cohen's _d_ for trait value: (z̄ₐ - z̄ₛ)/√P
	- z̄ₐ = analytical prediction
	- z̄ₛ = simulated value
- `dg` = cohen's _d_ for genetic value: (ḡₐ - ḡₛ)/√G
	- ḡₐ = analytical prediction
	- ḡₛ = simulated value
- `dm` = cohen's _d_ for microbial value: (m̄ₐ - m̄ₛ)/√M
	- m̄ₐ = analytical prediction
	- m̄ₛ = simulated value
- `corr` = corr(g,m)
"""
function rpNₑ(n)

    # if no dir make dir
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/Nₑ")
        mkdir("dat/Nₑ")
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
            fldr="dat/Nₑ/"
            )
        push!(Π,π)
    end

    rp(Π)

    convert("dat/Nₑ/",n)

    print("\nNₑ done\n")

    cmmd1 = `echo "Has finished simulating across Nₑ"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))

end

"""
simulate and analyze results across randomly drawn number of microbe species `S`
- `n` is number of random draws
- if no local folder `dat`, `rpS` creates this folder
	- if no subfolder `S`, `rpS` creates this folder
- outputs results to csv file in `/dat/S`

columns of outputted csv
- `S` = number of microbe species
- `Ne`= host population size
- `Je`= number of microbes in microbiome
- `κ`  = collective inheritance
- `κℓ`= total microbial inheritance (κ + ℓ)
- `k`  = replicate number
- `β` = realized selection gradient
- `P` = phenotypic variance
- `G` = additive genetic variance
- `M` = additive microbial variance
- `Dza` = analytically predicted change in mean trait
- `Dzs` = simulated change in mean trait
- `Dga` = analytically predicted change in mean genetic value
- `Dgs` = simulated change in mean genetic value
- `Dma` = analytically predicted change in mean microbial value
- `Dms` = simulated change in mean microbial value
- `dz` = cohen's _d_ for trait value: (z̄ₐ - z̄ₛ)/√P
	- z̄ₐ = analytical prediction
	- z̄ₛ = simulated value
- `dg` = cohen's _d_ for genetic value: (ḡₐ - ḡₛ)/√G
	- ḡₐ = analytical prediction
	- ḡₛ = simulated value
- `dm` = cohen's _d_ for microbial value: (m̄ₐ - m̄ₛ)/√M
	- m̄ₐ = analytical prediction
	- m̄ₛ = simulated value
- `corr` = corr(g,m)
"""
function rpS(n)

    # if no dir make dir
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/S")
        mkdir("dat/S")
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
            fldr="dat/S/"
            )
        push!(Π,π)
    end

    rp(Π)

    convert("dat/S/",n)

    print("\nS done\n")

    cmmd1 = `echo "Has finished simulating across S"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))
    
end

"""
simulate and analyze results across randomly drawn collective inheritance values `κ`
- `n` is number of random draws
- if no local folder `dat`, `rpκ` creates this folder
	- if no subfolder `κ`, `rpκ` creates this folder
- outputs results to csv file in `/dat/κ`

columns of outputted csv
- `S` = number of microbe species
- `Ne`= host population size
- `Je`= number of microbes in microbiome
- `κ`  = collective inheritance
- `κℓ`= total microbial inheritance (κ + ℓ)
- `k`  = replicate number
- `β` = realized selection gradient
- `P` = phenotypic variance
- `G` = additive genetic variance
- `M` = additive microbial variance
- `Dza` = analytically predicted change in mean trait
- `Dzs` = simulated change in mean trait
- `Dga` = analytically predicted change in mean genetic value
- `Dgs` = simulated change in mean genetic value
- `Dma` = analytically predicted change in mean microbial value
- `Dms` = simulated change in mean microbial value
- `dz` = cohen's _d_ for trait value: (z̄ₐ - z̄ₛ)/√P
	- z̄ₐ = analytical prediction
	- z̄ₛ = simulated value
- `dg` = cohen's _d_ for genetic value: (ḡₐ - ḡₛ)/√G
	- ḡₐ = analytical prediction
	- ḡₛ = simulated value
- `dm` = cohen's _d_ for microbial value: (m̄ₐ - m̄ₛ)/√M
	- m̄ₐ = analytical prediction
	- m̄ₛ = simulated value
- `corr` = corr(g,m)
"""
function rpκ(n)

    # if no dir make dir
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/κ")
        mkdir("dat/κ")
    end
    
    κ = rand(Uniform(),n)

    Π = Vector{parameters}()
    for r in 1:n
        π = parameters(
            ℓ=1-κ[r],
            κ=κ[r],
            K=1,
            k=r,
            fldr="dat/κ/")
        push!(Π,π)
    end

    rp(Π)

    convert("dat/κ/",n)

    print("\nκ done\n")

    cmmd1 = `echo "Has finished simulating across κ"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))
        
end

"""
simulate and analyze results across randomly drawn total microbial inheritance `κ+ℓ`
- `n` is number of random draws
- if no local folder `dat`, `rpκℓ` creates this folder
	- if no subfolder `κℓ`, `rpκℓ` creates this folder
- outputs results to csv file in `/dat/κℓ`

columns of outputted csv
- `S` = number of microbe species
- `Ne`= host population size
- `Je`= number of microbes in microbiome
- `κ`  = collective inheritance
- `κℓ`= total microbial inheritance (κ + ℓ)
- `k`  = replicate number
- `β` = realized selection gradient
- `P` = phenotypic variance
- `G` = additive genetic variance
- `M` = additive microbial variance
- `Dza` = analytically predicted change in mean trait
- `Dzs` = simulated change in mean trait
- `Dga` = analytically predicted change in mean genetic value
- `Dgs` = simulated change in mean genetic value
- `Dma` = analytically predicted change in mean microbial value
- `Dms` = simulated change in mean microbial value
- `dz` = cohen's _d_ for trait value: (z̄ₐ - z̄ₛ)/√P
	- z̄ₐ = analytical prediction
	- z̄ₛ = simulated value
- `dg` = cohen's _d_ for genetic value: (ḡₐ - ḡₛ)/√G
	- ḡₐ = analytical prediction
	- ḡₛ = simulated value
- `dm` = cohen's _d_ for microbial value: (m̄ₐ - m̄ₛ)/√M
	- m̄ₐ = analytical prediction
	- m̄ₛ = simulated value
- `corr` = corr(g,m)
"""
function rpκℓ(n)

    # if no dir make dir
    if !isdir("dat")
        mkdir("dat")
    end
    if !isdir("dat/κℓ")
        mkdir("dat/κℓ")
    end
    
    κℓ = rand(Uniform(),(2,n))

    Π = Vector{parameters}()
    for r in 1:n
        π = parameters(
            ℓ=κℓ[1,r]*κℓ[2,r],
            κ=κℓ[1,r]*(1-κℓ[2,r]),
            K=1,
            k=r,
            fldr="dat/κℓ/")
        push!(Π,π)
    end

    rp(Π)

    convert("dat/κℓ/",n)

    print("\nκℓ done\n")

    cmmd1 = `echo "Has finished simulating across κℓ"`
    cmmd2 = `mail -s "Cresko Cluster" bweek@uoregon.edu`
    run(pipeline(cmmd1,cmmd2))
    
end

#endregion

#
# functions that simulate and analyze results across fixed sets of parameters
#

#region

"""
simulate and analyze results across fixed set of host population sizes
- `Nₑs` is vector of host population sizes
- if no local folder `dat`, `ensNₑ` creates this folder
	- if no subfolder `Nₑ`, `ensNₑ` creates this folder
- outputs results to csv file in `/dat/Nₑ`

columns of outputted csv
- `S` = number of microbe species
- `Ne`= host population size
- `Je`= number of microbes in microbiome
- `κ`  = collective inheritance
- `κℓ`= total microbial inheritance (κ + ℓ)
- `k`  = replicate number
- `β` = realized selection gradient
- `P` = phenotypic variance
- `G` = additive genetic variance
- `M` = additive microbial variance
- `Dza` = analytically predicted change in mean trait
- `Dzs` = simulated change in mean trait
- `Dga` = analytically predicted change in mean genetic value
- `Dgs` = simulated change in mean genetic value
- `Dma` = analytically predicted change in mean microbial value
- `Dms` = simulated change in mean microbial value
- `dz` = cohen's _d_ for trait value: (z̄ₐ - z̄ₛ)/√P
	- z̄ₐ = analytical prediction
	- z̄ₛ = simulated value
- `dg` = cohen's _d_ for genetic value: (ḡₐ - ḡₛ)/√G
	- ḡₐ = analytical prediction
	- ḡₛ = simulated value
- `dm` = cohen's _d_ for microbial value: (m̄ₐ - m̄ₛ)/√M
	- m̄ₐ = analytical prediction
	- m̄ₛ = simulated value
- `corr` = corr(g,m)
"""
function ensNₑ(Nₑs)

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

"""
simulate and analyze results across fixed set of microbiome richnesses
- `Ss` is vector of microbiome richnesses
- if no local folder `dat`, `ensS` creates this folder
	- if no subfolder `S`, `ensS` creates this folder
- outputs results to csv file in `/dat/S`

columns of outputted csv
- `S` = number of microbe species (microbiome richness)
- `Ne`= host population size
- `Je`= number of microbes in microbiome
- `κ`  = collective inheritance
- `κℓ`= total microbial inheritance (κ + ℓ)
- `k`  = replicate number
- `β` = realized selection gradient
- `P` = phenotypic variance
- `G` = additive genetic variance
- `M` = additive microbial variance
- `Dza` = analytically predicted change in mean trait
- `Dzs` = simulated change in mean trait
- `Dga` = analytically predicted change in mean genetic value
- `Dgs` = simulated change in mean genetic value
- `Dma` = analytically predicted change in mean microbial value
- `Dms` = simulated change in mean microbial value
- `dz` = cohen's _d_ for trait value: (z̄ₐ - z̄ₛ)/√P
	- z̄ₐ = analytical prediction
	- z̄ₛ = simulated value
- `dg` = cohen's _d_ for genetic value: (ḡₐ - ḡₛ)/√G
	- ḡₐ = analytical prediction
	- ḡₛ = simulated value
- `dm` = cohen's _d_ for microbial value: (m̄ₐ - m̄ₛ)/√M
	- m̄ₐ = analytical prediction
	- m̄ₛ = simulated value
- `corr` = corr(g,m)
"""
function ensS(Ss)

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

"""
simulate and analyze results across fixed set of collective inheritances
- `κs` is vector of collective inheritances
- if no local folder `dat`, `ensκ` creates this folder
	- if no subfolder `κ`, `ensκ` creates this folder
- outputs results to csv file in `/dat/κ`

columns of outputted csv
- `S` = number of microbe species
- `Ne`= host population size
- `Je`= number of microbes in microbiome
- `κ`  = collective inheritance
- `κℓ`= total microbial inheritance (κ + ℓ)
- `k`  = replicate number
- `β` = realized selection gradient
- `P` = phenotypic variance
- `G` = additive genetic variance
- `M` = additive microbial variance
- `Dza` = analytically predicted change in mean trait
- `Dzs` = simulated change in mean trait
- `Dga` = analytically predicted change in mean genetic value
- `Dgs` = simulated change in mean genetic value
- `Dma` = analytically predicted change in mean microbial value
- `Dms` = simulated change in mean microbial value
- `dz` = cohen's _d_ for trait value: (z̄ₐ - z̄ₛ)/√P
	- z̄ₐ = analytical prediction
	- z̄ₛ = simulated value
- `dg` = cohen's _d_ for genetic value: (ḡₐ - ḡₛ)/√G
	- ḡₐ = analytical prediction
	- ḡₛ = simulated value
- `dm` = cohen's _d_ for microbial value: (m̄ₐ - m̄ₛ)/√M
	- m̄ₐ = analytical prediction
	- m̄ₛ = simulated value
- `corr` = corr(g,m)
"""
function ensκ(κs)

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

"""
simulate and analyze results across fixed set of microbiome sizes
- `Jₑs` is vector of microbiome sizes
- if no local folder `dat`, `ensκ` creates this folder
	- if no subfolder `Jₑ`, `ensκ` creates this folder
- outputs results to csv file in `/dat/Jₑ`

columns of outputted csv
- `S` = number of microbe species
- `Ne`= host population size
- `Je`= number of microbes in microbiome (microbiome size)
- `κ`  = collective inheritance
- `κℓ`= total microbial inheritance (κ + ℓ)
- `k`  = replicate number
- `β` = realized selection gradient
- `P` = phenotypic variance
- `G` = additive genetic variance
- `M` = additive microbial variance
- `Dza` = analytically predicted change in mean trait
- `Dzs` = simulated change in mean trait
- `Dga` = analytically predicted change in mean genetic value
- `Dgs` = simulated change in mean genetic value
- `Dma` = analytically predicted change in mean microbial value
- `Dms` = simulated change in mean microbial value
- `dz` = cohen's _d_ for trait value: (z̄ₐ - z̄ₛ)/√P
	- z̄ₐ = analytical prediction
	- z̄ₛ = simulated value
- `dg` = cohen's _d_ for genetic value: (ḡₐ - ḡₛ)/√G
	- ḡₐ = analytical prediction
	- ḡₛ = simulated value
- `dm` = cohen's _d_ for microbial value: (m̄ₐ - m̄ₛ)/√M
	- m̄ₐ = analytical prediction
	- m̄ₛ = simulated value
- `corr` = corr(g,m)
"""
function ensJₑ(Jₑs)

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

#endregion