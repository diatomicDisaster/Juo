"""
The Potentials module contains functional forms of the various potential energy curves 
available.
"""
module Potentials

using Plots

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

""" Functional form of the simple Morse potential energy curve with polynomial expansion parameters
    Arguments
        r :: Float64
        Internuclear distance at which to calculate the potential energy.
        r_e :: Float64
        Equilibrium internuclear distance (r coordinate of potential minimum).
        a :: Float64
        Coefficient that determines the width of the potential well.
        D_e :: Float64
        Dissociation energy (potential well depth).
        V_0::Float64
        Energy  of the potential at equillibrium (r_e), in comparison to the overall model.
        expansion_parameters :: Vector{Float64} 
        The polynomial expansion parameters for the morse potential
    Returns
        V   :: Vector{Float}
            Vector of potential energy values at each internuclear distance.
"""
function morse_polynomial(r::Float64; r_e::Float64=1.0, a::Float64=1.0, D_e::Float64=1.0, V_0::Float64=0.0, expansion_parameters::AbstractVector{Float64}=Float64[])
    y=1 - exp(-a*(r - r_e))
    V = D_e*y^2
    if isempty(expansion_parameters) == false
        for( i, k) in enumerate(expansion_parameters)
            V=V+k*y^(i+3)
        end
    end
    return V
end

""" Functional form of the dampened Morse potential energy curve with polynomial expansion parameters
    Arguments
        r :: Float64
        Internuclear distance at which to calculate the potential energy.
        r_e :: Float64
        Equilibrium internuclear distance (r coordinate of potential minimum).
        a :: Float64
        Coefficient that determines the width of the potential well.
        damp::Float64
        The dampening constant
        D_e :: Float64
        Dissociation energy (potential well depth).
        V_0::Float64
        Energy  of the potential at equillibrium (r_e), in comparison to the overall model.
        expansion_parameters :: Vector{Float64} 
        The polynomial expansion parameters for the dampened morse potential
    Returns
        V   :: Vector{Float}
            Vector of potential energy values at each internuclear distance.
"""
function morse_damp_polynomial(r::Float64; r_e::Float64=1.0, a::Float64=1.0, damp::Float64=0, D_e::Float64=1.0, V_0::Float64=0.0, expansion_parameters::AbstractVector{Float64}=Float64[])
    y=1 - exp(-a*(r - r_e))
    V_long = D_e*y^2
    V_damp=exp(-damp*(r-r_e))
    z=(r-r_e)/(r+r_e)
    V=0
    if isempty(expansion_parameters) == false
        for( i, k) in enumerate(expansion_parameters)
            V=V+k*z^(i+2)
        end
    end
    V=V_0+V*V_damp+V_long
    return V
end

""" Functional form of the Modified Morse potential energy curve with polynomial expansion parameters
    Arguments
        r :: Float64
        Internuclear distance at which to calculate the potential energy.
        r_e :: Float64
        Equilibrium internuclear distance (r coordinate of potential minimum).
        D_e :: Float64
        Dissociation energy (potential well depth).
        V_0::Float64
        Energy  of the potential at equillibrium (r_e), in comparison to the overall model.
        expansion_parameters :: Vector{Float64} 
        The polynomial expansion parameters for the modified morse potential
    Returns
        V   :: Vector{Float}
            Vector of potential energy values at each internuclear distance.
"""
function morse_modified(r::Float64; r_e::Float64=1.0, D_e::Float64=1.0, V_0::Float64=1.0, expansion_parameters::AbstractVector{Float64}=Float64[])
    z=(r-r_e)/(r+r_e)
    Denominator=sum(expansion_parameters)
    Numerator=0
    if isempty(expansion_parameters) == false
        for( i, k) in enumerate(expansion_parameters)
            Numerator=Numerator+k*z^(i+1)
        end
    end
    V=V_0+D_e*(1-exp(-Numerator))^2/(1-exp(-Denominator))^2
    return V
end

end