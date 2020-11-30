"""
The Potentials module contains functional forms of the various potential energy curves 
available.
"""
module PotentialCurves

export mor_pot

"""Functional form of the simple Morse potential energy curve.
    Arguments
        rs  :: Vector{Float}
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
function mor_pot(rs::Vector{Float64}; a::Float64=1.0, r_e::Float64=1.0, D_e::Float64=1.0)
    V = Vector{Float64}(undef, size(rs, 1))
    for (i, r) in enumerate(rs)
        v[i] = D_e*(1 - exp(-a*(r - r_e)))^2
    end
    return v
end

end
