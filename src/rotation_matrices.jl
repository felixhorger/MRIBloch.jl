
# TODO: Rotations.jl

# General axis
function rotation_matrix(axis::NTuple{3, <: Real}, α::Real; out=Matrix{Float64}(undef, 3, 3))
	# Arbitrary axis
	# Apply vector from the right
	(sinα, cosα) = sincos(α)
	one_minus_cosα = (1 - cosα)
	# First column
	out[1, 1] = cosα + axis[1]^2 * one_minus_cosα
	out[2, 1] = axis[1] * axis[2] * one_minus_cosα + axis[3] * sinα
	out[3, 1] = axis[1] * axis[3] * one_minus_cosα - axis[2] * sinα
	# Second column
	out[1, 2] = axis[1] * axis[2] * one_minus_cosα - axis[3] * sinα
	out[2, 2] = cosα + axis[2]^2 * one_minus_cosα
	out[3, 2] = axis[2] * axis[3] * one_minus_cosα + axis[1] * sinα
	# Third column
	out[1, 3] = axis[1] * axis[3] * one_minus_cosα + axis[2] * sinα
	out[2, 3] = axis[2] * axis[3] * one_minus_cosα - axis[1] * sinα
	out[3, 3] = cosα + axis[3]^2 * one_minus_cosα
	return out
end


# Around axis in plane alinged with coordinate system
function rotation_matrix(plane::Val{:xy}, axis::NTuple{2, <: Real}, α::Real; out=Matrix{Float64}(undef, 3, 3))
	# axis in xy plane, axis must be normalised!
	# Apply vector from the right
	(sinα, cosα) = sincos(α)
	one_minus_cosα = (1 - cosα)
	# First column
	out[1, 1] = cosα + axis[1]^2 * one_minus_cosα
	out[2, 1] = axis[1] * axis[2] * one_minus_cosα
	out[3, 1] = -axis[2] * sinα
	# Second column
	out[1, 2] = axis[1] * axis[2] * one_minus_cosα
	out[2, 2] = cosα + axis[2]^2 * one_minus_cosα
	out[3, 2] = axis[1] * sinα
	# Third column
	out[1, 3] = axis[2] * sinα
	out[2, 3] = -axis[1] * sinα
	out[3, 3] = cosα
	return out
end

# Around axis in plane alinged with coordinate system
function rotation_matrix(plane::Val{:xz}, axis::NTuple{2, <: Real}, α::Real; out=Matrix{Float64}(undef, 3, 3))
	# Arbitrary axis
	# Apply vector from the right
	(sinα, cosα) = sincos(α)
	one_minus_cosα = (1 - cosα)
	# First column
	out[1, 1] = cosα + axis[1]^2 * one_minus_cosα
	out[2, 1] = axis[2] * sinα
	out[3, 1] = axis[1] * axis[2] * one_minus_cosα
	# Second column
	out[1, 2] = -axis[2] * sinα
	out[2, 2] = cosα
	out[3, 2] = axis[1] * sinα
	# Third column
	out[1, 3] = axis[1] * axis[2] * one_minus_cosα
	out[2, 3] = -axis[1] * sinα
	out[3, 3] = cosα + axis[2]^2 * one_minus_cosα
	return out
end


# Axis defined by phase
function rotation_matrix(plane::Val{:xy}, θ::Real, α::Real; out=Matrix{Float64}(undef, 3, 3))
	# Rotate by α, around axis specified by θ
	# Right handed coordinate system, x,y,z
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

function rotation_matrix(plane::Val{:xz}, θ::Real, α::Real; out=Matrix{Float64}(undef, 3, 3))
	# Rotate by α, around axis specified by θ
	# Right handed coordinate system, x,y,z
	# This is for B1 in xz plane
	# Apply vector from the right
	(sinα, cosα) = sincos(α)
	(sinθ, cosθ) = sincos(θ)
	# First column
	out[1, 1] = sinθ^2 + cosα * cosθ^2
	out[2, 1] = -sinα * cosθ
	out[3, 1] = sinθ * cosθ - cosα * sinθ * cosθ
	# Second column
	out[1, 2] = sinα * cosθ
	out[2, 2] = cosα
	out[3, 2] = -sinα * sinθ
	# Third column
	out[1, 3] = sinθ * cosθ - cosα * sinθ * cosθ
	out[2, 3] = sinα * sinθ
	out[3, 3] = cosθ^2 + cosα * sinθ^2

	return out
end

