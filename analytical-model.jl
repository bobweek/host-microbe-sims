#############################
#							#
# plot host trait evolution #
# across transmission modes #
#							#
#############################

include("defs.jl");

using RCall

i = 0
ℓs = [0, 0.25, 0.5, 0.75, 1.0]

for ℓ in ℓs

	# restrict parameter values
	# such that ℓ + κ ≤ 1
	κs = [0, 0.25, 0.5, 0.75, 1.0]
	if ℓ == 0.25
		κs = [0.0, 0.25, 0.5, 0.75]
	elseif ℓ == 0.5
		κs = [0.0, 0.25, 0.5]
	elseif ℓ == 0.75
		κs = [0.0, 0.25]
	elseif ℓ == 1.0
		κs = [0.0]
	end

	for κ in κs

		if κ == 1.0
			ℓ = 0.0
		end		

		as = asystem(ℓ = ℓ, κ = κ, G = 0, T = 5, k = i)
		arun(as)
		i += 1

	end

end

R"source('analytical-trans-mode-figure.r')"