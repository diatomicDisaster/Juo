"""
The Potentials module contains functional forms of the various potential energy curves 
available.
"""
module Potentials

"""Functional form of the simple Morse potential energy curve.
    Arguments
        r :: Float64
            Internuclear distance at which to calculate the potential energy.
        r_e :: Float64
            Equilibrium internuclear distance (r coordinate of potential minimum).
        a :: Float64
            Coefficient that determines the width of the potential well.
        D_e :: Float64
            Dissociation energy (potential well depth).
    Returns
        V :: Vector{Float}
            Potential energy at the given internuclear distance.
"""
function morse(r::Float64; r_e::Float64=1.0, a::Float64=1.0, D_e::Float64=1.0)
    return D_e*(1 - exp(-a*(r - r_e)))^2
end

"""Functional form of the quantum harmonic oscillator potential energy curve.
    Arguments
        r :: Float64
            Internuclear distance at which to calculate the potential energy.
        r_e :: Float64
            Equilibrium internuclear distance (r coordinate of potential minimum).
        w :: Float64
            Frequency of the harmonic oscillator (determines width of potential well).
        m :: Float64
            Reduced mass of the harmonic oscillator.
    Returns
        V :: Float64
            Potential energy at the given internuclear distance.
"""
function qho(r::Float64; r_e::Float64=0.0, w::Float64=1.0, m::Float64=1.0)
    return .5*m*(w*(r - r_e))^2
end

end