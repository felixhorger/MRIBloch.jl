
module MRIBloch
	
	using LinearAlgebra

	include("rotation_matrices.jl")

	# TODO: MAKE GENERATED?

	function simulate(mode::Val{:pulse_frame}, m0::NTuple{3, <: Real}, ω1::AbstractVector{<: Real}, Δω0::AbstractVector{<: Real}, dt::Real)
		# Memory for magnetisation
		m = Array{Float64}(undef, 3, 1+length(Δω0))
		m[:, 1] .= m0
		# Precompute flip angles
		α = @. dt * sqrt(Δω0^2 + ω1^2)
		# Run
		R = Matrix{Float64}(undef, 3, 3)
		for t = 1:length(α)
			# Compute rotation axis and matrix
			axis = (ω1[t], -Δω0[t])
			axis = axis ./ norm(axis)
			rotation_matrix(Val(:xz), axis, α[t]; out=R)
			# Apply
			@views m[:, t+1] = R * m[:, t]
		end
		return m
	end

	#function simulate(mode::Val{:B0_frame}, m0::NTuple{3, <: Real}, ω1::AbstractVector{<: Real}, ϕ::AbstractVector{<: Real}, dt::Real)
	#	# Memory for magnetisation
	#	m = Array{Float64}(undef, 3, 1+length(ϕ))
	#	m[:, 1] .= m0
	#	# Run
	#	R = Matrix{Float64}(undef, 3, 3)
	#	for t = 1:length(α)
	#		# Compute rotation axis and matrix
	#		α = dt * ω1[t]
	#		rotation_matrix(Val(:xy), ϕ[t], α; out=R)
	#		# Apply
	#		@views m[:, t+1] = R * m[:, t]
	#	end
	#	return m
	#end

	function simulate(mode::Val{:B0_frame}, m0::NTuple{3, <: Real}, ω1::AbstractVector{<: Complex}, dt::Real)
		# Memory for magnetisation
		m = Array{Float64}(undef, 3, 1+length(ω1))
		m[:, 1] .= m0
		# Run
		R = Matrix{Float64}(undef, 3, 3)
		for t = 1:length(ω1)
			# Compute rotation axis and matrix
			abs_ω1 = abs(ω1[t])
			α = dt * abs_ω1
			axis = ω1[t] / abs_ω1
			rotation_matrix(Val(:xy), reim(axis), α; out=R)
			# Apply
			@views m[:, t+1] = R * m[:, t]
		end
		return m
	end
end

