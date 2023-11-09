################################
#                              #
# definitions of functions and #
# data structs for worm-box.jl #
#                              #
################################

using Parameters, LinearAlgebra, Random, Distributions, JLD2

function flr(x::Float64)
    return floor(Int,x)
end

@with_kw mutable struct ModelParameters

	# simulation level
	T::Int = flr(1e3)        # number microbe generations to run model
	tres::Int = 50           # temporal resolution for saving data

	# system level
	t::Float64  = 0              # transmission is purely environmental
	S::Int  = 100           # number of microbe species
	L::Int  = 100            # host lifetime in microbe generations
	K::Int  = 100           # number of hosts
	O::Vector{Float64} = fill(1/S,S)    # relative microbe abundances outside the system
	Jₒ::Int = 100            # size of immigrating microbe pool into system

	# environmental microbiome
	EJ::Int = flr(1e4)  # size of environmental microbiome

	# host microbiome
	HJ::Int = flr(1e3)   # size of host microbiome
	Jₕ::Int = 5         # number microbes host expels
	Jₑ::Int = 10         # number microbes host acquires

	# stddev of microbial effects on host fitness (normally distr, mean zero)
	sd::Float64 = 1e-8

end

@with_kw mutable struct Host

	# parameters	
	f::Vector{Float64}	# per-capita microbe effect on host fitness
	J::Int64			# host microbiome size
	Jₕ::Int				# size of migrating pool from host to env
	Jₑ::Int				# size of migrating pool from env to host

	# variables
    p::Vector{Float64}	# relative abundances in host
    # age::Int64			# host age (for deterministic or geo branching times)

end

@with_kw mutable struct Environment

	# environmental microbiome size
	J::Int64

	# relative abundances in environment
	p::Vector{Float64}

end

@with_kw mutable struct System

	# system-level parameters
	t::Float64			# 0 => pure env, 1 => max vert
	S::Int				# number of microbe species	
	L::Float64			# expected host lifetime in microbe gens
	K::Int				# number of hosts
	O::Vector{Float64}	# outer microbiome	
	Jₒ::Int				# size of migrating pool from outside to env

	# system components
	E::Environment		# environmental microbiome
	H::Vector{Host}	# hosts
	age::Int = 0		# system age

end

@with_kw mutable struct Sim

	T::Int64 = 1000  # number microbe gens to run sim for
	tres::Int64 = 10 # temporal resolution for exporting data

	sys::System

end

# deterministic microbial movement
function deter_move(srcp,tgtp,tgtJ,Jₘ)

	# compute immigrant weight
	q = Jₘ/(tgtJ+Jₘ)

	# compute relative abundances of target community
	tgtp = q*srcp .+ (1-q)*tgtp

	return tgtp
end

function hubbell_step(J,p)

	mnl = Multinomial(J,p)
	p = normalize(rand(mnl),1)
	
	return p
end



function host_reproduction(H,E,t)	

	K = length(H)

	# compute host fitnesses
	W = Vector{Float64}()
	for h in H
		push!(W, exp(dot(h.f,h.p)))
	end
	normalize!(W,1)

	# draw number offspring per host
	mnl = Multinomial(K,W)
	Off = rand(mnl)

	# create offspring hosts
	OH = Vector{Host}()
	for i in 1:K
		for j in 1:Off[i]

			# offspring microbiome is weighted average
			# of parental microbiome and environmental
			Op = t*H[i].p .+ (1-t)*E.p
			push!(OH, Host(
				f  = deepcopy(H[i].f),
				J  = deepcopy(H[i].J),
				Jₕ = deepcopy(H[i].Jₕ),
				Jₑ = deepcopy(H[i].Jₑ),
				p  = deepcopy(Op)
				))

		end
	end

	return OH
end

# one microbe generation
function microbe_gen(sys::System)
	@unpack_System sys

	E.p = deter_move(O, E.p, E.J, Jₒ)
	E.p = hubbell_step(E.J,E.p)

	# to be used for moving microbes from hosts to env
	pooled_Jₕ = 0
	pooled_p = zeros(K)	

	for h in H		

		# microbes move from env into host
		h.p = deter_move(E.p, h.p, h.J, h.Jₑ)

		# hubbelian dynamics play out in host
		h.p = hubbell_step(h.J,h.p)

		# host ages by one microbe generation
		# h.age += 1

		# accumulate relative abundances of
		# microbes leaving hosts
		pooled_p += h.p .* h.Jₕ

		# accumulate number of
		# microbes leaving hosts
		pooled_Jₕ += h.Jₕ

	end

	# microbes move from hosts to env
	E.p = deter_move(normalize(pooled_p,1), E.p, E.J, pooled_Jₕ)
	
	# if system age is divisible by host lifespan
	if mod(age,L)==0
		H = host_reproduction(H,E,t)
	end

	# entire systems ages by one microbe generation
	age += 1
	
	@pack_System! sys
	return sys
end

function make_sim(P::ModelParameters)
	@unpack_ModelParameters P

	# draw one set of microbial fitness effects for all hosts
	nml = Normal(0,sd)
	f   = rand(nml,S)

	# make host vector
	h = Host(f=f,J=HJ,Jₕ=Jₕ,Jₑ=Jₑ,p=O)
	H = Vector{Host}()
	for k in 1:K
		push!(H,deepcopy(h))
	end

	# make environment
	E = Environment(J=EJ,p=O)

	# make system
	sys = System(t=t,S=S,L=L,K=K,O=O,Jₒ=Jₒ,E=E,H=H)

	# make simulation
	sim = Sim(sys=deepcopy(sys))

	return sim
end

function run_sim(sim::Sim)
	@unpack_Sim sim

	hist = Vector{System}()

	for time in 1:T

		# append history
		push!(hist,deepcopy(sys))

		# iterate model
		sys = microbe_gen(sys)

	end

	return hist
end

# unfinished
# save only relative abundances
# function hist2csv(hist)

# 	K = Vector{Float64}()
# 	Nˢ= Vector{Vector{Float64}}()

# 	# export things in csv's
# 	for sys::System in hist
# 		push!(K,sys.K)
# 		push!(Nˢ,sys.S.N)
# 	end

# end

# unfinished
# stochastic microbe movement
# ex usage: comm1, comm2 = stoch_move(comm1,comm2,100)
# function stoch_move(src::Communittee,tgt::Communittee,Jₘ::Float64)

# 	# take random sample from source community as immigrating pool
# 	imm_mnl = Multinomial(flr(Jₘ),src.p)
# 	imm_pl	= rand(imm_mnl)

# 	# compute updated relative abundances for target
# 	tgt.p = normalize(abs_abund(tgt)+imm_pl,1)

# 	# compute updated relative abundances for source
# 	abs_src = abs_abund(src)
# 	# first, scale imm_pl if needed
# 	# (otherwise abs_src .- imm_pl could have neg entries)
# 	mx, ind = findmax(imm_pl .- abs_src)
# 	if mx > 0
# 		imm_pl *= abs_src/imm_pl[ind]
# 	end
# 	src.p = normalize(abs_src .- imm_pl,1)

# 	return tgt, src
# end
