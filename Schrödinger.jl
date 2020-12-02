module Vibrational

using LinearAlgebra

"""Kinetic energy operator T_ij for the sinc-DVR method on the interval [0, inf].
    See: D. Colbert and W. Miller, J. Chem. Phys., 96 (1992), doi:10.1063/1.462100
"""
function sinc_dvr_radial(N::Int64)
    T = Array{Float64}(undef, (N, N))
    dx = 1.0/N
    for j = 1:N
        for i = 1:j-1
            T[i, j] = (1 / (2 * m * dx^2)) * ((-1)^(i - j))*( (2 / (i - j)^2) - (2 / (i + j)^2) )
        end
        T[j, j] = (1 / (2 * m * dx^2)) * ((pi^2 / 3) - 1/(2*j^2))
    end
    return Symmetric(T)
end

"""Kinetic energy operator T_ij for the sinc-DVR method on the interval [-inf, inf].
    See: D. Colbert and W. Miller, J. Chem. Phys., 96 (1992), doi:10.1063/1.462100
"""
function sinc_dvr_kinetic(rs, interval::Float64, m::Float64)
    T = Array{Float64}(undef, (length(rs), length(rs)))
    for j = 1:length(rs)
        for i = 1:j-1
            T[j, i] = (1 / (2*m * interval^2)) * 2*( ((-1)^(i - j)) / (i - j)^2 )
            T[i, j] = T[j, i]
        end
        T[j, j] = (1 / (2*m * interval^2)) * (pi^2 / 3)
    end
    return T
end


end