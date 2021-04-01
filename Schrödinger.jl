"""
The Vibrational module contains methods for solving the vibrational Schrödinger equation.
"""
module Vibrational

using LinearAlgebra

"""Kinetic energy operator T_ij for the sinc-DVR method on the interval [-inf, inf].
    See: D. Colbert and W. Miller, J. Chem. Phys., 96 (1992), doi:10.1063/1.462100
"""
function kinetic_operator(npoints::Int, interval::Float64)
    operator = Array{Float64}(undef, (npoints, npoints))
    for j = 1:npoints
        for i = 1:j-1
            operator[i, j] = (1 / (2* interval^2)) * 2*( ((-1)^(i - j)) / (i - j)^2 )
        end
        operator[j, j] = (1 / (2* interval^2)) * (pi^2 / 3)
    end
    return Symmetric(operator)
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
function solve(;
    potential_grid :: Array{Float64, 1} = nothing, 
    interval       :: Float64 = nothing, 
    threshold      :: Float64 = -1.0, 
    mass           :: Float64 = 1.0
    )
    if threshold < 0.
        potential_operator = Diagonal(potential_grid)
    else
        usepoints = isless.(potential_grid, threshold)
        potential_operator = Diagonal(potential_grid[usepoints])
    end
    kinetic_operator = sincdvr_kinetic(size(potential_operator)[1], interval)
    hamiltonian = potential_operator + kinetic_operator/mass
    vibrational_eigen = eigen(hamiltonian)
    return (vibrational_eigen, hamiltonian)
end

end

module Rotational

using LinearAlgebra
using SparseArrays

function raising_operator(;
    Jmax :: Float64 = 0.0
    )
    operator = sparse([1], [1], [0.])
    for ang_mom = 1:Jmax
        block = spdiagm(1 => [(ang_mom*(ang_mom + 1) - ang_mom_proj*(ang_mom_proj + 1))^.5 for ang_mom_proj=-ang_mom:ang_mom-1])
        operator = blockdiag(operator, block)
    end
    return operator
end

function integrate(bra, ket, interval)
    integrand = sum(interval * (bra + ket)/2)
    return integrand
end
end