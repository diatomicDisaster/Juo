abstract type AbstractBasis end

struct Electronic{N<:Number, L<:Real} <: AbstractBasis
    nstate::N 
    lambda::L
end

struct Spin{T<:Real} <: AbstractBasis
    ess::T
    sigma::T
end

struct Vibrational{N<:Number, L<:Real} <: AbstractBasis
    nstate::N
    vee::Real
end

struct Rotational{T<:Real} <: AbstractBasis
    jay::T
    omega::T
end

struct Rovibrational <: AbstractBasis
    e :: Electronic
    s :: Spin
    v :: Vibrational
    r :: Rotational
end

abstract type BasisSet end
