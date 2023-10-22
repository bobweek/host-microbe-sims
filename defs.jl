########################
#                      #
# mbiome explicit sims #
#                      #
########################

using Parameters, LinearAlgebra, Random, Distributions, DataFrames, CSV

#region

# computes ε as a function of ℓ and κ
# used in both analytical and simulation models
function ε₀(ℓ,κ)

	if ℓ < 1
		return 1 - κ / (1 - ℓ)
	else
		return 1
	end

end

# check inheritance parameters
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

@with_kw mutable struct asystem

	# model parameters
	ℓ::Float64	= (3-√5)/2	# lineal inheritance
	κ::Float64	= (3-√5)/2	# collective inheritance
	ε::Float64	= ε₀(ℓ,κ)	# prb of env ancestry
	E::Float64	= 0.0		# variance of noise
	G::Float64	= 1.0		# additive genetic variance
	M::Float64	= 1.0		# additive microbial variance
	P::Float64	= G+M+E		# phenotypic variance
	β::Float64	= 1.0		# selection gradient (= s under model assumptions in main text)
	T::Int64	= 50		# number host gens to iterate
	k::Int64	= 0			# parameter combination

	# state variables
	ḡ::Float64 = 0.0		 # mean additive genetic value
	m̄::Float64 = 0.0		 # mean additive microbial value
	z̄::Float64 = ḡ+m̄		   # mean trait value
	ξ::Float64 = 0.0		# env additive microbial value
	t::Int64   = 0			# host generation

end

# Δḡ = Gβ
# Δm̄ = qMβ+(1-q)(ξ-m̄)
# Δξ = (1-ε)(Mβ+m̄-ξ)
function Δ(as::asystem)
	@unpack_asystem as

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

# save system state and parameter values
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

	if (k == 0) & (t == 0)
		CSV.write("a_dat.csv", dat)
	else
		CSV.write("a_dat.csv", dat, append=true)
	end

	return
end

# run the analytical model
function arun(as::asystem)

	# make sure ℓ + κ ≤ 1
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

# variance of per-capita additive effects (see main text)
function σ²(S,Jₑ,L)
	return S/(Jₑ^2*(S-1)*((1-(1-1/Jₑ)^L)))
end

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
	Jₑᴱ::Int64	= 10000		# number of microbes in environment
	Jₑᴹ::Int64	= 10000		# number of microbes in a host
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

@with_kw mutable struct host

	α::Vector{Float64}	# per-capita microbe effect on host trait
    f::Vector{Float64}	# microbe rel abundances in host microbiome
	p::Int64   = 0		# index of host parent in previous generation (zero => first gen)
	g::Float64 = 0.0	# additive genetic value
	m::Float64 = 0.0	# additive microbial value
	c::Float64 = 0.0	# additive microbial offset (s.t. m̄₀=ḡ₀=0)
	z::Float64 = 0.0	# trait
	w::Float64 = 1.0	# realized relative fitness
	W::Int64   = 1		# realized absolute fitness

	# note: using realized instead of expected fitnesses
	#			because it should make analytical model more accurate
	#			when calculating β = cov(w,z)/P

end

@with_kw mutable struct system

	fᴱ::Vector{Float64}	# microbe rel abundances in environment
	c::Float64 = 0.0	# env microbial offset (s.t. ξ₀=0)
	H::Vector{host}		# hosts
	t::Int64 = 0		# host generation

end

# a microbiome dynamic
function hubbell(Jₑ,f)
	f = abs.(f)
	normalize!(f,1)
	f = normalize(rand(Multinomial(Jₑ,f)),1)
	return f
end

# host fitness
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

# host reproduction
function hrep(sys::system,π::parameters)
	@unpack_system sys
	@unpack_parameters π

	# form pooled mbiome of selected hosts
	fᴾ = zeros(S)
	for h in H
		fᴾ .+= h.W .* h.f
	end
	normalize!(fᴾ,1)
	
	# mutations
	μˢ = rand(Normal(0,√μ),Nₑ)

	# form host offspring
	Hₒ = Vector{host}()
	i = 1
	k = 1
	for h in H
		
		# offspring of host i
		for j in 1:h.W

			# sample from env
			fᴱᴾ = normalize(rand(Multinomial(aqsz,fᴱ)),1)
			
			# sample from pool
			fᴾᴾ = normalize(rand(Multinomial(aqsz,fᴾ)),1)

			# init offspring mbiome
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

	# env post shedding
	fᴾᴾ = normalize(rand(Multinomial(shsz,fᴾ)),1)
	fᴱ = (ε .* fᴱ) .+ ((1-ε) .* fᴾᴾ)

	@pack_system! sys
	return sys
end

# one host generation
function hgen(sys::system,π::parameters)
	
	# host reproduction
	sys = hrep(sys,π)	

	@unpack_system sys
	@unpack_parameters π

	# host development
	for i in 1:L

		# env mbiome dyn
		fᴱ = hubbell(Jₑᴱ,fᴱ)

		# hst mbiome dyn
		for h in H
			h.f = hubbell(Jₑᴹ,h.f)
		end
		
	end
	
	# host traits
	for h in H
		# Jₑᴹ makes α per capita effect
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
	end

	# note:
	# 	because sum(Nₒ)/Nₑ = 1
	#	Nₒ/Nₑ are the realized rel fit's of each adult host

	@pack_system! sys
	return sys
end

# generate initial condition
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

function ind_df(sys::system,π::parameters)
	@unpack_system sys
	@unpack_parameters π

	hi    = Vector{Int64}()
	pi    = Vector{Int64}()
	gs    = Vector{Float64}()
	ms    = Vector{Float64}()
	zs    = Vector{Float64}()
	Ws	  = Vector{Float64}()
	αs	  = Array{Float64}(undef,0,S)
	fs	  = Array{Float64}(undef,0,S)
	types = Vector{String}()
	
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

function pop_df(sys::system,π::parameters)

	idat = ind_df(sys,π)

	envdat = idat[idat.type .== "env", :]
	hstdat = idat[idat.type .== "hst", :]

	P = var(hstdat.z)
	G = var(hstdat.g)
	M = var(hstdat.m)

	return DataFrame(
		z = mean(hstdat.z),
		g = mean(hstdat.g),
		m = mean(hstdat.m),
		ξ = envdat.m[1],
		P = P,
		G = G,
		M = M,
		corr = cov(hstdat.g,hstdat.m)/sqrt(G*M),
		β = cov(hstdat.W,hstdat.z)/P,
		t = hstdat.t[1],
		k = hstdat.k[1]
	)

	# note: lande & arnold use rel fit to calculate β,
	#			which is what i do above as well

end

# function savepar(π::parameters)
# 	@unpack_parameters π
# 	par = DataFrame(
# 		sig = σ,
# 		ell = ℓ,
# 		kap = κ,
# 		eps = ε,
# 		L   = L,
# 		S   = S,
# 		Ne  = Nₑ,
# 		s   = s,
# 		mu  = μ,
# 		JeE = Jₑᴱ,
# 		JeM = Jₑᴹ,
# 		K   = K,
# 		T   = T,
# 		B   = B
# 	)
# 	CSV.write(fldr*"par_"*fname*".csv", par)
# end

function popcat(sys::system,π::parameters,popdf::DataFrame)
	pdf = pop_df(sys,π)
	return vcat(popdf,pdf)
end

# simulate system
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

# run ensemble
function ens(π::parameters)

	rps = DataFrame()
	for k in 1:π.K
		rps = vcat(rps,π2df(π))
	end
	CSV.write(π.fldr*"par_"*π.fname*".csv", rps)

	repdfs = Vector{DataFrame}(undef, π.K)
	Threads.@threads for k in 1:π.K
		local πₚ = deepcopy(π)
		πₚ.k = k
		repdfs[k] = sim(πₚ)
	end

	pdat = repdfs[1]
	for k in 2:π.K
		pdat = vcat(pdat,repdfs[k])
	end

	CSV.write(π.fldr*"pop_dat_"*π.fname*".csv", pdat)

end

function rp(Π::Vector{parameters})

	rps = DataFrame()
	for π in Π
		rps = vcat(rps,π2df(π))
	end

	CSV.write(Π[1].fldr*"rndpars.csv", rps)

	pdfs = Vector{DataFrame}(undef, length(Π))
	Threads.@threads for π in Π
		local k = π.k
		pdfs[k] = sim(π)
	end

	pdat = DataFrame()
	for π in Π
		pdat = vcat(pdat,pdfs[π.k])
	end

	CSV.write(Π[1].fldr*"rndpar_dat.csv", pdat)

end

#endregion

#
# simulation vs analytical model
#

#region

@with_kw mutable struct astate
	ḡ::Float64 = 0.0	 # mean additive genetic value
	m̄::Float64 = 0.0	 # mean additive microbial value
	z̄::Float64 = 0.0	 # mean trait value
	ξ::Float64 = 0.0	# env additive microbial value
end

function read(fldr::String,fname::String)

	pdat = CSV.read(fldr*"pop_dat_"*fname*".csv", DataFrame)
	par = CSV.read(fldr*"par_"*fname*".csv", DataFrame)

	return pdat, par
	
end

function read(fldr::String)

	pdat = CSV.read(fldr*"rndpar_dat.csv", DataFrame)
	par = CSV.read(fldr*"rndpars.csv", DataFrame)

	return pdat, par
	
end

function cohen(pdat::DataFrame,par::DataFrame,t::Int64,k::Int64)

	subdat = pdat[pdat.k .== k, :]

	# current
	crtdat = subdat[subdat.t .== t, :]

	if length(crtdat.z) == 0
		print("Error: attempted to read at t > T")
		exit()
	end

	# next
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

# make csv with columns equal to S, Nₑ, ℓ, Corr(m,g), z̄ᵈ, ḡᵈ, m̄ᵈ, k
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

function convert(fldr::String,n::Int64)

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