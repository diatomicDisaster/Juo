"""
The Vibrational module contains methods for solving the vibrational Schrödinger equation.
"""
module Vibrational

using LinearAlgebra

"""Kinetic energy operator T_ij for the sinc-DVR method on the interval [-inf, inf].
    See: D. Colbert and W. Miller, J. Chem. Phys., 96 (1992), doi:10.1063/1.462100
"""
function build_sincdvr_kinetic(npoints::Int, interval::Float64)
    T = Array{Float64}(undef, (npoints, npoints))
    for j = 1:npoints
        for i = 1:j-1
            T[i, j] = T[j, i] = (1 / (2* interval^2)) * 2*( ((-1)^(i - j)) / (i - j)^2 )
        end
        T[j, j] = (1 / (2* interval^2)) * (pi^2 / 3)
    end
    return T
end

"""Calculate eigenvectors and eigenvalues of the sinc-DVR Hamiltonian for a grid of 
potential energy values.
    Arguments
        potential_grid :: AbstractVector{Float64}
            An array containing the value of the potential energy at each grid point.
        interval :: Float64
            The spatial interval Δr for the vibrational grid.
        threshold :: Float64
            Energy threshold above which to discard grid points.
        mass :: Float64
            The reduced mass of the molecule, appearing the kinetic energy operator.
    Returns
        vibrational_eigen :: LinearAlgebra.Eigen
            Julia LinearAlgebra.Eigen object containing the eigenvectors and eigenvalues
            of the sinc-DVR Hamiltonian matrix.
        hamiltonian_matrix :: AbstractMatrix{Float64}
            Hamiltonian matrix in the sinc-DVR basis (i.e prior to diagonalisation).
""" 
function sincdvr(;
    potential_grid::AbstractVector{Float64}=nothing, 
    interval::Float64=nothing, 
    threshold::Float64=-1.0, 
    mass::Float64=1.0
    )
    if threshold < 0.
        potential_matrix = Diagonal(potential_grid)
    else
        usepoints = isless.(potential_grid, threshold)
        potential_matrix   = Diagonal(potential_grid[usepoints])
    end
    kinetic_matrix     = build_sincdvr_kinetic(size(potential_matrix)[1], interval)
    hamiltonian_matrix = potential_matrix + kinetic_matrix/mass
    vibrational_eigen  = eigen(hamiltonian_matrix)
    return (vibrational_eigen, hamiltonian_matrix)
end

end