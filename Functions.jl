"""
The Potentials module contains functional forms of the various potential energy curves 
available.
"""
module PotentialCurves

export mor_pot

"""Functional form of the Morse potential.
"""
function mor_pot(r::Float64; a::Float64=1.0, r_e::Float64=1.0, D_e::Float64=1.0)
    return D_e*(1-exp(-a*(r - r_e)))^2
end

function mor_pot(rs::Vector{Float64}; a::Float64=1.0, r_e::Float64=1.0, D_e::Float64=1.0)
    v = Vector{Float64}(undef, size(rs, 1))
    for (i, r) in enumerate(rs)
        v[i] = D_e*(1-exp(-a*(r - r_e)))^2
    end
    return v
end

end
