"""
The Vibrational module contains methods for solving the vibrational Schrödinger equation.
"""
module Vibrational

using LinearAlgebra

"""Build the vibrational Hamiltonian for a uniformly spaced potential energy 
grid using the sinc-DVR kinetic energy operator.

For details of the sinc-DVR method, see: 
D. Colbert and W. Miller, J. Chem. Phys., 96 (1992), doi:10.1063/1.462100
"""
function sincdvr(
    potential_grid :: Array{Float64, 1}, 
    sep :: Float64,
    reduced_mass :: Float64
    )
    npoints = length(potential_grid)
    kinetic_operator = Array{Float64}(undef, (npoints, npoints))
    for j = 1:npoints
        for i = 1:j-1
            kinetic_operator[i, j] = (1 / (2*reduced_mass*sep^2)) * 2*( ((-1)^(i - j)) / (i - j)^2 )
        end
        kinetic_operator[j, j] = (1 / (2*reduced_mass*sep^2)) * (pi^2 / 3)
    end
    return Symmetric(kinetic_operator) + Diagonal(potential_grid)
end

end

module Rotational

using LinearAlgebra
using SparseArrays
using Arpack
using BenchmarkTools

"""
Builds the raising operator matrix for the basis of angular momentum quantum numbers provided.
    Arguments
        Jlist :: Array{Float64, 1}
            List of angular moment quantum numbers to build the matrix for, can be non-sequential.
        ndimensions :: Int64
            The dimension of the rotational basis. For a complete set of rotational quantum 
            numbers, this is equal to (J_max + 1)^2 if J is integer, or Jmax(Jmax + 3) + 2 
            if J is half-integer.
    Returns
        SparseArrays.AbstractSparseMatrixCSC{Float64, 2}
            CSC sparse matrix representing the raising operator for the rotational basis.
"""
function raising_operator(
    Jlist :: Array{Float64, 1},
    ndimensions :: Int64
    )
    rows = Int64[i for i=2:ndimensions]
    cols = Int64[i for i=1:ndimensions-1]
    vals = Array{Float64, 1}(undef, ndimensions-1)
    i = 1
    for J in Jlist
        for M = -J:J-1
            vals[i] = (J*(J + 1) - M*(M + 1))^.5
            i += 1
        end
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function raising_operator(
    J :: Float64,
    ndimensions :: Int64
    )
    rows = Int64[i for i=2:ndimensions]
    cols = Int64[i for i=1:ndimensions-1]
    vals = Array{Float64, 1}(undef, ndimensions-1)
    i = 1
    for M = -J:J-1
        vals[i] = (J*(J + 1) - M*(M + 1))^.5
        i += 1
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end


"""
Builds the lowering operator matrix for the basis of angular momentum quantum numbers provided.
    Arguments
        Jlist :: Array{Float64, 1}
            List of angular moment quantum numbers to build the matrix for, can be non-sequential.
        ndimensions :: Int64
            The dimension of the rotational basis. For a complete set of rotational quantum 
            numbers, this is equal to (J_max + 1)^2 if J is integer, or Jmax(Jmax + 3) + 2 
            if J is half-integer.
    Returns
        SparseArrays.AbstractSparseMatrixCSC{Float64, 2}
            CSC sparse matrix representing the lowering operator for the rotational basis.
"""
function lowering_operator(
    Jlist :: Array{Float64, 1},
    ndimensions :: Int64
    )
    rows = Int64[i for i=1:ndimensions-1]
    cols = Int64[i for i=2:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions-1)
    i = 0
    for J in Jlist
        for M = -J:J
            if i > 0
                vals[i] = (J*(J + 1) - M*(M - 1))^.5
                i += 1
            end
        end
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function lowering_operator(
    J :: Float64,
    ndimensions :: Int64
    )
    rows = Int64[i for i=2:ndimensions]
    cols = Int64[i for i=1:ndimensions-1]
    vals = Array{Float64, 1}(undef, ndimensions-1)
    i = 0
    for M = -J:J
        if i > 0
            vals[i] = (J*(J + 1) - M*(M - 1))^.5
            i += 1
        end
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end


"""
Builds the squared momentum operator matrix for the basis of angular momentum quantum numbers provided.
    Arguments
        Jlist :: Array{Float64, 1}
            List of angular moment quantum numbers to build the matrix for, can be non-sequential.
        ndimensions :: Int64
            The dimension of the rotational basis. For a complete set of rotational quantum 
            numbers, this is equal to (J_max + 1)^2 if J is integer, or Jmax(Jmax + 3) + 2 
            if J is half-integer.
    Returns
        SparseArrays.AbstractSparseMatrixCSC{Float64, 2}
            CSC sparse matrix representing the squared momentum operator for the rotational basis.
"""
function momentum_operator(
    Jlist :: Vector{Float64},
    ndimensions :: Int64
    )
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions)
    i = 1
    for J in Jlist
        for M = -J:J
            vals[i] = J*(J + 1)
            i += 1
        end
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function momentum_operator(
    J :: Float64,
    ndimensions :: Int64
    )
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions)
    i = 1
    for M = -J:J
        vals[i] = J*(J + 1)
        i += 1
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end


"""
Builds the momentum projection operator matrix for the basis of angular momentum quantum numbers provided.
    Arguments
        Jlist :: Array{Float64, 1}
            List of angular moment quantum numbers to build the matrix for, can be non-sequential.
        ndimensions :: Int64
            The dimension of the rotational basis. For a complete set of rotational quantum 
            numbers, this is equal to (J_max + 1)^2 if J is integer, or Jmax(Jmax + 3) + 2 
            if J is half-integer.
    Returns
        AbstractSparseMatrixCSC{Float64, 2}
            CSC sparse matrix representing the momentum projection operator for the rotational basis.
"""
function momentum_projection_operator(
    Jlist :: Array{Float64, 1},
    ndimensions :: Int64
    )
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions)
    i = 1
    for J in Jlist
        for M = -J:J
            vals[i] = M
            i += 1
        end
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function momentum_projection_operator(
    J :: Float64,
    ndimensions :: Int64
    )
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions)
    i = 1
    for M = -J:J
        vals[i] = M
        i += 1
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end


"""
Builds the rigid-rotor Hamiltonian for the specified list of angular momentum quantum numbers,
and given value of total spin angular momentum.
    Arguments
        Jlist :: Array{Float64, 1}
            List of angular moment quantum numbers to build the matrix for, can be non-sequential.
        spin :: Float64
            Total spin angular momentum for the system.
        mass :: Float64
            Reduced mass of the molecule.
        r :: Float64
            Rigid rotor length (internuclear distance).
    Returns
        SparseArrays.AbstractSparseMatrixCSC{Float64, 2}
            CSC sparse matrix representing the momentum projection operator for the rotational basis.
"""
function hamiltonian(
    jaylist :: Array{Float64, 1},
    spin  :: Float64,
    reduced_mass  :: Float64,
    r     :: Float64
    )

    ndimensions = _countdimensions(jaylist)
    ndimensions_spin = floor(Int64, 2*spin + 1)
    
    #Build the J, M operators in the tensor product basis
    identity = Matrix{Float64}(1.0I, (ndimensions_spin, ndimensions_spin))
    jaysq = momentum_operator(jaylist, ndimensions)
    em = momentum_projection_operator(jaylist, ndimensions)
    jaysq_spinid = kron(jaysq, identity)
    em_spinid = kron(em, identity)

    #Build the S, Σ operators in the tensor product basis
    identity = Matrix{Float64}(I, (ndimensions, ndimensions))
    spinsq = momentum_operator(spin, ndimensions_spin)
    sigma = momentum_projection_operator(spin, ndimensions_spin)
    jayid_spinsq = kron(identity, spinsq)
    jayid_sigma = kron(identity, sigma)

    #Build S-uncoupling operators in the tensor product basis
    jaypl = raising_operator(jaylist, ndimensions)
    jaymi = lowering_operator(jaylist, ndimensions)
    spinpl = raising_operator(spin, ndimensions_spin)
    spinmi = lowering_operator(spin, ndimensions_spin)
    jaypl_spinmi = kron(jaypl, spinmi)
    jaymi_spinpl = kron(jaypl, spinmi)

    #Return sum of Hamiltonian terms
    return ((jaysq_spinid - em_spinid.^2) + (jayid_spinsq - jayid_sigma.^2) + (jaypl_spinmi + jaymi_spinpl))/(2*reduced_mass*r^2)
end

end

module Vibronic

include("data.jl")
struct VibrationalGrid
    rmin :: Float64
    rmax :: Float64
    sep  :: Float64
    num  :: Int64
    vals :: Vector{Float64}
end

function hamiltonian(
    sincdvr_sep :: VibrationalGrid,
    potentials :: AbstractArray{Array{Float64, 1}, 1};
    lambda :: AbstractArray{Float64, 1} = [0.0],
    spin   :: AbstractArray{Float64, 1} = [0.0],
    mu     :: Float64 = 1.0
    )
    vib_hamils = Vibrational.sincdvr.(potentials, sincdvr_sep)
    vib_eigens = eigen.(vib_hamils)

end




end