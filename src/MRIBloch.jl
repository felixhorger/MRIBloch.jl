
module MRIBloch
	
	using LinearAlgebra

	# TODO: Rotations.jl
	function rotation_matrix(α::Real, n::NTuple{3, <: Real}; out=Matrix{Float64}(undef, 3, 3))
		# Arbitrary axis
		# Apply vector from the right
		(sinα, cosα) = sincos(α)
		one_minus_cosα = (1 - cosα)
		# First column
		out[1, 1] = cosα + n[1]^2 * one_minus_cosα
		out[2, 1] = n[1] * n[2] * one_minus_cosα + n[3] * sinα
		out[3, 1] = n[1] * n[3] * one_minus_cosα - n[2] * sinα
		# Second column
		out[1, 2] = n[1] * n[2] * one_minus_cosα - n[3] * sinα
		out[2, 2] = cosα + n[2]^2 * one_minus_cosα
		out[3, 2] = n[2] * n[3] * one_minus_cosα + n[1] * sinα
		# Third column
		out[1, 3] = n[1] * n[3] * one_minus_cosα + n[2] * sinα
		out[2, 3] = n[2] * n[3] * one_minus_cosα - n[1] * sinα
		out[3, 3] = cosα + n[3]^2 * one_minus_cosα
		return out
	end

	function rotation_matrix(α::Real, θ::Real; out=Matrix{Float64}(undef, 3, 3))
		# This is for B1 in xy plane
		# Apply vector from the right
		(sinα, cosα) = sincos(α)
		(sinθ, cosθ) = sincos(θ)
		# First column
		out[1, 1] = sinθ^2 + cosα * cosθ^2
		out[2, 1] = sinθ * cosθ - cosα * sinθ * cosθ
		out[3, 1] = -sinα * cosθ
		# Second column
		out[1, 2] = sinθ * cosθ - cosα * sinθ * cosθ
		out[2, 2] = cosθ^2 + cosα * sinθ^2
		out[3, 2] = sinα * sinθ
		# Third column
		out[1, 3] = sinα * cosθ
		out[2, 3] = -sinα * sinθ
		out[3, 3] = cosα

		return out
	end


	function simulate(m0, Δω0, ω1, dt)
		# Memory for magnetisation
		m = Array{Float64}(undef, 3, 1+length(Δω0))
		m[:, 1] = m0
		# Precompute flip angles
		α = @. dt * sqrt(Δω0^2 + ω1^2)
		# Run
		R = Matrix{Float64}(undef, 3, 3)
		for t = 1:length(α)
			# Compute rotation axis and matrix
			n = (ω1[t], 0.0, -Δω0[t])
			n = n ./ norm(n)
			rotation_matrix(α[t], n; out=R)
			# Apply
			m[:, t+1] = R * @view m[:, t]
		end
		return m
	end

end

