abstract type AbstractBasis end
Base.iterate(x::AbstractBasis) = (x, nothing)
Base.iterate(x::AbstractBasis, ::Any) = nothing
Base.isempty(x::AbstractBasis) = false
in(x::AbstractBasis, y::AbstractBasis) = x == y

"""Vibrational basis state for a given vibrational quantum number `vee` represented
by a `vector` on the `Hamiltonian` `grid`, with eigenvalue `energy`.
"""
struct Vib <: AbstractBasis
    vee    :: Int64
    energy :: Float64
    vector :: Vector{Float64}
end

"""Electronic basis state for a given `state` with electronic angular momentum 
projection `lambda`, spin `ess` and vibrational basis `vib`.
"""
struct Ele <: AbstractBasis
    state   :: Int64
    lambda  :: Float64
    ess     :: Float64
end

struct Vibronic <: AbstractBasis
    ele :: Ele
    vib :: Vib
end

struct Rot <: AbstractBasis
    jay   :: Float64
    omega :: Float64
end

struct Rovibronic <: AbstractBasis
    rot :: Rot
    vib :: Vib
    ele :: Ele
end