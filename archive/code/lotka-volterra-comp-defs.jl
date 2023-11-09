################################
#                              #
# definitions of functions and #
# data structs for worm-box.jl #
#                              #
################################

using Parameters, LinearAlgebra, Random, Distributions, LinkedLists, JLD2

@with_kw mutable struct Environment

	# parameters
	λ::Vector{Float64}  # expected number immigrating microbes
	ψ₀::Vector{Float64} # base probability microbe survival in environment
	Cₑ::Vector{Float64} # microbe interspecific competition in environment
	Cₐ::Vector{Float64} # microbe intraspecific competition in environment
	
	# variables
	N::Vector{Int64}	# microbe abundances
	
end

@with_kw mutable struct Host

	# parameters
	J::Int64 = 1e5		# size of microbe community in host
	L::Float64 = 10     # expected host lifetime in microbe gens
	f::Vector{Float64}  # per-capita microbe effect on host fitness
	ε::Vector{Float64}  # per-microbe shedding probability	
	ψ₀::Vector{Float64} # base probability microbe survival in host
	Cₑ::Vector{Float64} # microbe interspecific competition in host
	Cₐ::Vector{Float64} # microbe intraspecific competition in host

	# variables
	N::Vector{Int64}	# microbe abundances in host	
	age::Int64 = 0		# host age (for deterministic or geo branching times)

end

@with_kw mutable struct System

	# system-level parameters
	t::Float64 = 0			# 0 => pure env, 1 => max vert
	S::Int64 = 100			# number of microbe species	
	K::Int64  = 100			# number of hosts
	maxFrac::Float64 = 0.1	# maximum fraction of substrate microbes acquired

	# system components
	E::Substrate			# microbe abundances in environment
	H::LinkedList{Host}		# linked list of hosts
	age::Int64 = 0			# system age

end

@with_kw mutable struct Model

	T::Int64 = 1000  # number microbe gens to run sim for
	tres::Int64 = 10 # temporal resolution for exporting data

	sys::System

end

# computes post competition abundance in a given compartment for a given microbe species
# this includes microbe reproduction
function microbe_comp(h::Host)
	@unpack_Host h

	CommN = sum(N)
	ψ = ψ₀.*(Cₐ.^N).*(Cₑ.^(CommN.-N))
	# ψ[ψ.>1].=1
	bin = Binomial.(N,ψ)
	N = reduce(vcat,2*rand.(bin))
	
	@pack_Host! h
	return h
end

function microbe_comp(S::Substrate)
	@unpack_Substrate S

	CommN = sum(N)
	ψ = ψ₀.*(Cₐ.^N).*(Cₑ.^(CommN.-N))
	bnml = Binomial.(N,ψ)
	N = reduce(vcat,2*rand.(bnml))
	
	@pack_Substrate! S
	return S
end

function comp_loop(sys::System)
	@unpack_System sys

	for h::Host in H
		h = microbe_comp(h)
	end

	@pack_System! sys
	return sys
end


# returns post immigration microbe abundance in substrate
function microbe_imm(S::Substrate)
	@unpack_Substrate S

	# add incoming microbes to substrate abundances
	N += floor.(Int64,λ) # replace with poisson if want more noise
	
	@pack_Substrate! S
	return S
end

function substr_stuff(sys::System)
	@unpack_System sys

	S = microbe_imm(S)
	S = microbe_comp(S)

	@pack_System! sys
	return sys
end

# returns microbe abundances after hosts acquire microbes
function host_acquire(sys::System)
	@unpack_System sys
	
	# compute total number microbes acquired
	p = K/K₀ # proportion microbes acquired by hosts
	p = min(p,1)
	p *= maxFrac
	Na = floor.(Int64,p*S.N)

	# microbe abundances in substrate after acquisition by hosts
	S.N -= Na

	# divide acquired microbes among hosts
	mnd = Multinomial.(Na,fill(fill(1/K,K),n))
	Nacq = reduce(vcat,transpose.(rand.(mnd))) # (n×K)

	# host linked list with updated microbe abundances
	i = 0
	for h::Host in H
		i += 1
		h.N += Nacq[:,i]
	end

	# updated microbe abundances for whole system
	@pack_System! sys
	return sys
end

# computes microbe abundances after hosts shed microbes
function host_shed(sys::System)
	@unpack_System sys

	for h::Host in H
		bmd = Binomial.(h.N,h.ε)
		Nε = rand.(bmd)
		S.N += Nε
		h.N -= Nε
	end

	# updated microbe abundances for whole system
	@pack_System! sys
	return sys
end

function host_branching(sys::System)
	@unpack_System sys

	for index in keys(H)
		h = getindex(H,index)
		if h.age >= floor(Int64,h.L)
			# compute number of offspring
			offMean = exp(dot(h.N,h.f))
			offDistr = Poisson(offMean)
			W = rand(offDistr)
			W = min(W,K₀-K+1)
			if W > 0
				if t == 0 # environmental transmission

					# host parent microbes go into substrate
					S.N += h.N

					# host offspring begin w no microbes
					for i in 1:W
						pushfirst!(H, Host(L=h.L, f=h.f, ε=h.ε, ψ₀=h.ψ₀, Cₑ=h.Cₑ, Cₐ=h.Cₐ, N=zeros(n)))
					end					

				else # strict vertical transmission

					# host parent microbes distributed among offspring

				end
			else
				# host microbes go into substrate
				S.N += h.N
			end
			# remove parental host
			deleteat!(H,index)
		end
		K = length(H)
	end	

	# updated microbe abundances for whole system
	@pack_System! sys
	return sys
end

function host_reproduction(sys::System)
	@unpack_System sys


end


function host_age(sys::System)
	@unpack_System sys

	# age each host by one
	for h::Host in H
		h.age += 1
	end

	# updated microbe abundances for whole system
	@pack_System! sys
	return sys
end

function host_shuffle(sys::System)
	@unpack_System sys

	hvec = Vector{Host}(undef,K)
	i = 0
	for h::Host in H
		i += 1
		hvec[i] = deepcopy(h)
	end
	shuffle!(hvec)
	i = 0
	for h::Host in H
		i += 1
		h = deepcopy(hvec[i])
	end

	@pack_System! sys
	return sys
end

# one microbe generation
function microbe_gen(sys::System)

	# perform model components
	sys = substr_stuff(sys)
	if sys.K > 0
		sys = host_age(sys)	
		sys = host_acquire(sys)
		sys = comp_loop(sys)
		sys = host_shed(sys)
		sys = host_branching(sys)
		sys = host_shuffle(sys)
	end
	sys.age += 1
	
	return sys
end

function run_model(model::Model)
	@unpack_Model model

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
function hist2csv(hist)

	K = Vector{Float64}()
	Nˢ= Vector{Vector{Float64}}()

	# export things in csv's
	for sys::System in hist
		push!(K,sys.K)
		push!(Nˢ,sys.S.N)
	end

end