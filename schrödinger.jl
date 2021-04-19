"""
Contains operators that appear in diatomic Born-Oppenheimer Schrödinger equation.
"""
module Hamiltonian

include("data.jl")

using LinearAlgebra
using SparseArrays

"""Build the vibrational Hamiltonian for a uniformly spaced potential energy 
grid using the sinc-DVR kinetic energy operator.

For details of the sinc-DVR method, see: 
D. Colbert and W. Miller, J. Chem. Phys., 96 (1992), doi:10.1063/1.462100
"""
function sincdvr(potential_grid::Array{Float64, 1}, dr::Float64, mu::Float64)
    npoints = length(potential_grid)
    kinetic_operator = Array{Float64}(undef, (npoints, npoints))
    for j = 1:npoints
        for i = 1:j-1
            kinetic_operator[i, j] = (1 / (2*mu*dr^2)) * 2*( ((-1)^(i - j)) / (i - j)^2 )
        end
        kinetic_operator[j, j] = (1 / (2*mu*dr^2)) * (pi^2 / 3)
    end
    return Symmetric(kinetic_operator) + Diagonal(potential_grid)
end

"""
Integrate the vibrational matrix element of a given function over the grid.
    """
function integrate(left::Array{Rovibronic}, grid_values::Array{Float64}, right::Array{Rovibronic}; dr::Float64=1.0)
    return dr*sum(left .* grid_values .* right)
end

"""
Builds the raising operator matrix for a given angular momentum quantum number.
"""
function raising(J :: Float64)
    ndimensions = floor(Int64, 2*J + 1)
    rows = Int64[i for i=2:ndimensions]
    cols = Int64[i for i=1:ndimensions-1]
    vals = Array{Float64, 1}(undef, ndimensions-1)
    for (i, M) in enumerate(-J:J-1)
        vals[i] = (J*(J + 1) - M*(M + 1))^.5
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

"""
Builds the lowering operator matrix for a given angular momentum quantum number.
"""
function lowering(J :: Float64)
    ndimensions = floor(Int64, 2*J + 1)
    rows = Int64[i for i=2:ndimensions]
    cols = Int64[i for i=1:ndimensions-1]
    vals = Array{Float64, 1}(undef, ndimensions-1)
    for (i, M) in enumerate(-J+1:J)
        vals[i] = (J*(J + 1) - M*(M - 1))^.5
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

"""
Builds the total momentum operator matrix for a given angular momentum quantum number.
"""
function momentum_squared(J :: Float64)
    ndimensions = floor(Int64, 2*J + 1)
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions)
    for (i, M) in enumerate(-J:J)
        vals[i] = J*(J + 1)
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

"""
Builds the momentum projection operator matrix for a given angular momentum quantum number.
"""
function momentum_projection(J :: Float64)
    ndimensions = floor(Int64, 2*J + 1)
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions)
    for (i, M) in enumerate(-J:J)
        vals[i] = M
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

"""
Build the vibronic basis set for a list of potential energy curves.
"""
function vibronic_basis(
    veemax       :: Int64,
    potential    :: Vector{Float64},
    separation   :: Float64, 
    reduced_mass :: Float64
)
    nstates = length(potential_list)
    vibbasis = Vector{Vibronic}(undef, veemax + 1)
    vibhamiltonian = sincdvr(potential, separation, reduced_mass)
    vibeigen = eigen(vibhamiltonian)
    for (v, vee) in enumerate(0:veemax_list[p])
        i = sum(veemax_list[1:p - 1]) + v
        val = vibeigen.values[v]
        vec = vibeigen.vectors[:,v]
        vibbasis[i] = Vibrational(vee, val, vec)
    end
    return vibbasis
end

"""
Build the rovibronic basis for a given total angular momentum quantum number and list of electronic quantum numbers.
"""
function rovibronic_basis(
    jay         :: Vector{Float64},
    veemax_list :: Vector{Int64},
    lambda_list :: Vector{Float64},
    ess_list    :: Vector{Float64}
)
    length(veemax_list) == length(lambda_list) == length(ess_list) || throw(ArgumentError("Length of veemax_list (=$(length(veemax_list))), lambda_list (=$(length(lambda_list))), and ess_list (=$(length(ess_list))) should be equal."))
    vibronicdimen_list = Array{Int64}(undef, length(veemax_list))
    @. vibronicdimen_list = (veemax_list + 1) * floor(Int64, 2*ess_list + 1)
    basis = Vector{Rovibronic}(undef, sum(vibronicdimen_list))
    for state = 1:nstates
        lambda = lambda_list[state]
        ess    = ess_list[state]
        veemax = veemax_list[state]
        for (s, sigma) in enumerate(-ess:ess)
            omega = lambda + sigma
            for vee = 0:veemax
                i = (j - 1) * (sum(vibronicdimen_list)) + sum(vibronicdimen_list[1:(state-1)]) + (s-1)*(veemax+1) + (vee+1)
                basis[i] = Rovibronic(state, vee, lambda, ess, sigma, jay, omega)
            end
        end
    end
    return basis
end



end