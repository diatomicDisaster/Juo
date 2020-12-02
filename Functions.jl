"""
The Potentials module contains functional forms of the various potential energy curves 
available.
"""
module Potentials

export morse

"""Functional form of the simple Morse potential energy curve.
    Arguments
        r  :: Float64
            The vector of internuclear distances at which to calculate the potential energy.
        a   :: Float64
            Coefficient that determines the width of the potential well.
        r_e :: Float64
            The equilibrium internuclear distance (r coordinate of potential minimum).
        D_e :: Float64
            The dissociation energy (potential well depth).
    Returns
        V   :: Vector{Float}
            Vector of potential energy values at each internuclear distance.
"""
function morse(r::Float64; a::Float64=1.0, r_e::Float64=1.0, D_e::Float64=1.0, mu::Float64=1.0)
    return D_e*(1 - exp(-a*(r - r_e)))^2
end

end