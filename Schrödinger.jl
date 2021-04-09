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
function hamiltonian(
    potential_grid :: Array{Float64, 1}, 
    interval :: Float64,
    reduced_mass :: Float64
    )
    npoints = length(potential_grid)
    kinetic_operator = Array{Float64}(undef, (npoints, npoints))
    for j = 1:npoints
        for i = 1:j-1
            kinetic_operator[i, j] = (1 / (2*reduced_mass*interval^2)) * 2*( ((-1)^(i - j)) / (i - j)^2 )
        end
        kinetic_operator[j, j] = (1 / (2*reduced_mass*interval^2)) * (pi^2 / 3)
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
    Jlist :: Vector{Float64},
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
    Jlist :: Vector{Float64},
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
Counts the number of dimensions, equivalent to the sum of (2J + 1), for a specified list of
angular momentum quantum numbers. For a complete series of quantum numbers up to Jmax, this
summation can be evaluated analytically.
"""
function _countdimensions(
    Jlist :: Array{Float64, 1}
    )
    ndimensions = 0
    for J in Jlist
        ndimensions += floor(Int64, 2*J + 1)
    end
    return ndimensions
end

function _countdimensions(
    Jmax :: Float64
    )
    nJ = floor(Int64, Jmax)
    if mod(2*Jmax, 2) == 0
        return (nJ + 1)^2
    else
        return nJ*(nJ + 3) + 2
    end
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
    Jlist :: Array{Float64, 1},
    spin  :: Float64,
    reduced_mass  :: Float64,
    r     :: Float64
    )

    ndimensions = _countdimensions(Jlist)
    ndimensions_spin = floor(Int64, 2*spin + 1)
    
    identity = Matrix{Float64}(1.0I, (ndimensions_spin, ndimensions_spin))
    J2 = momentum_operator(Jlist, ndimensions)
    Jz = momentum_projection_operator(Jlist, ndimensions)
    J2SI = kron(J2, identity)
    JzSI = kron(Jz, identity)

    identity = Matrix{Float64}(I, (ndimensions, ndimensions))
    S2 = momentum_operator(spin, ndimensions_spin)
    Sz = momentum_projection_operator(spin, ndimensions_spin)
    JIS2 = kron(identity, S2)
    JISz = kron(identity, Sz)

    Jp = raising_operator(Jlist, ndimensions)
    Jm = lowering_operator(Jlist, ndimensions)
    Sp = raising_operator(spin, ndimensions_spin)
    Sm = lowering_operator(spin, ndimensions_spin)
    JpSm = kron(Jp, Sm)
    JmSp = kron(Jm, Sp)
    return ((J2SI - JzSI.^2) + (JIS2 - JISz^2) + (JpSm + JmSp))/(2*reduced_mass*r^2)
end


function integrate(bra, ket, interval)
    integrand = sum(interval * (bra + ket)/2)
    return integrand
end

end