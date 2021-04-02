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

function raising_operator(
    Jmax :: Float64
    )
    nJ = floor(Int64, Jmax)
    if mod(2*Jmax, 2) == 0
        ndimensions = (nJ + 1)^2
        Jmin = .0
    else
        ndimensions = nJ*(nJ + 3) + 2
        Jmin = .5
    end
    rows = Int64[i for i=2:ndimensions]
    cols = Int64[i for i=1:ndimensions-1]
    vals = Float64[
        (J*(J + 1) - M*(M + 1))^.5 for J=Jmin:Jmax for M=-J:J if M < Jmax
        ]
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function lowering_operator(
    Jmax :: Float64
    )
    nJ = floor(Int64, Jmax)
    if mod(2*Jmax, 2) == 0
        ndimensions = (nJ + 1)^2
        Jmin = .0
    else
        ndimensions = nJ*(nJ + 3) + 2
        Jmin = .5
    end
    rows = Int64[i for i=1:ndimensions-1]
    cols = Int64[i for i=2:ndimensions]
    vals = Float64[
        (J*(J + 1) - M*(M - 1))^.5 for J=Jmin:Jmax for M=-J:J if (J != Jmin || M !=-J)
        ]
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function momentum_operator(
    Jmax :: Float64
    )
    nJ = floor(Int64, Jmax)
    if mod(2*Jmax, 2) == 0
        ndimensions = (nJ + 1)^2
        Jmin = .0
    else
        ndimensions = nJ*(nJ + 3) + 2
        Jmin = .5
    end
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Float64[
        J*(J + 1) for J=Jmin:Jmax for M=-J:J
    ]
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function momentum_projection_operator(
    Jmax :: Float64
    )
    nJ = floor(Int64, Jmax)
    if mod(2*Jmax, 2) == 0
        ndimensions = (nJ + 1)^2
        Jmin = .0
    else
        ndimensions = nJ*(nJ + 3) + 2
        Jmin = .5
    end
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Float64[
        M for J=(Jmax - nJ):Jmax for M=-J:J
    ]
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end


function spin_raising_operator(
    spin :: Float64
    )
    ndimensions = 2*spin + 1
    rows = Int64[i for i=2:ndimensions]
    cols = Int64[i for i=1:ndimensions-1]
    vals = Float64[
        (spin*(spin + 1) - M*(M + 1))^.5 for M=-spin:spin if M < spin
        ]
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function spin_lowering_operator(
    spin :: Float64
    )
    ndimensions = 2*spin + 1
    rows = Int64[i for i=2:ndimensions]
    cols = Int64[i for i=1:ndimensions-1]
    vals = Float64[
        (spin*(spin + 1) - M*(M - 1))^.5 for M=-spin:spin if M < spin
        ]
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function spin_operator(
    spin :: Float64
    )
    ndimensions = 2*spin + 1
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Float64[spin*(spin + 1) for M=-spin:spin]
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function spin_projection_operator(
    spin :: Float64
    )
    ndimensions = 2*spin + 1
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Float64[M for M=-spin:spin]
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function hamiltonian(
    Jmax :: Float64,
    S    :: Float64,
    mass :: Float64,
    r    :: Float64
    )

    nJ = floor(Int64, Jmax)
    if mod(2*Jmax, 2) == 0
        ndimensions = (nJ + 1)^2
    else
        ndimensions = nJ*(nJ + 3) + 2
    end
    
    identity = Matrix{Float64}(1.0I, (Int64(2*S + 1), Int64(2*S + 1)))
    J2 = kron(momentum_operator(Jmax), identity)
    Jz = kron(momentum_projection_operator(Jmax), identity)

    identity = Matrix{Float64}(I, (ndimensions, ndimensions))
    S2 = kron(identity, spin_operator(S))
    Sz = kron(identity, spin_projection_operator(S))

    JpSm = kron(raising_operator(Jmax), spin_lowering_operator(S))
    JmSp = kron(lowering_operator(Jmax), spin_raising_operator(S))
    
    return ((J2 - Jz^2) + (S2 - Sz^2) + (JpSm + JmSp))/(2*mass*r^2)
end

function integrate(bra, ket, interval)
    integrand = sum(interval * (bra + ket)/2)
    return integrand
end
end