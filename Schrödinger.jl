module Vibrational

"""Kinetic energy operator T_ij for the sinc-DVR method on the interval [0, inf].

See: D. Colbert and W. Miller, J. Chem. Phys., 96 (1992), doi:10.1063/1.462100
"""
function sinc_dvr(N::Int64)
    T = Array{Float64}(undef, (N, N))
    dx = 1.0/N
    for i = 1:N
        for j = 1:N
            if i == j
                T[i, j] = ((-1/2)^(i-j)) * ((pi^2/3) - 1/(2*i^2))
            else
                T[i, j] = ((-1/(2*dx^2))^(i-j)) * ((2/(i-j)^2) - (2/(i+j)^2))
            end
        end
    end
    return T
end
end