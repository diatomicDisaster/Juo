using LinearAlgebra
using SparseArrays
using Dierckx

include("operators.jl")
include("basis.jl")
include("grids.jl")

"""Hamiltonian object for basis sets and eigenproblems."""
mutable struct Hamiltonian
    grid  :: Vector{Float64}
    mu    :: Float64
    potentials :: Vector{Potential}
    couplings  :: Vector{Coupling}
    basisfrompot   :: IdDict{Potential, Vector{Vib}}
    Hamiltonian(grid, mu, args...) =
        !all(p -> p≈diff(grid[1:2])[1], diff(grid)) ? 
        error("Grid not uniformly spaced") : 
        new(grid, mu, args...)
end

Hamiltonian(grid, mu) = Hamiltonian(mu, grid, Vector{Potential}(undef, 0), Vector{Potential}(undef, 0), IdDict())

"""Add `Potential` to `Hamiltonian`, interpolates `values` onto `grid` using cubic splines."""
function (ham::Hamiltonian)(pot::Potential)
    spl = Spline1D(pot.nodes, pot.values)
    interpvals = spl(ham.grid)
    interppoten = Potential(ham.grid, interpvals, pot.lambda, pot.ess) #interpolate onto grid
    push!(ham.potentials, interppoten)
    return interppoten
end

"""Add `Coupling` to `Hamiltonian`, interpolates `values` onto `grid` using cubic splines."""
function (ham::Hamiltonian)(coupl::Coupling)
    # !!! currently a coupling is dissociated from its potentials when interpolating, might need to modify couplings and potentials in-place
    spl = Spline1D(coupl.nodes, coupl.values)
    interpvals = spl(ham.grid)
    interpcoupl = Coupling(ham.grid, interpvals, coupl.left, coupl.right)
    push!(ham.couplings, interpcoupl)
    return interpcoupl
end

Main.push!(ham::Hamiltonian, poten::Potential) = ham(poten)
Main.push!(ham::Hamiltonian, coupl::Coupling) = ham(coupl)

"""
Integrate the vibrational matrix element of a given function over the grid.
"""
function integrate(leftvib::Vib, coupl::Coupling, rightvib::Vib)
    length(leftvib.vector) == length(coupl) == length(rightvib.vector) || throw(ArgumentError("Integrands ($leftvib, $coupl, $rightvib) must be of equal length"))
    dr = diff(coupl.nodes[1:2])[1]
    integral = dr*sum(leftvib.vector.*coupl.values.*rightvib.vector)
    return integral
end

function integrate!(ham::Hamiltonian)
    rsep = diff(ham.grid[1:2])[1]
    for coupl in ham.couplings
        leftvibs = ham.basisfrompot[coupl.left]
        rightvibs = ham.basisfrompot[coupl.right]
        for i in CartesianIndices(coupl.vibop)
            leftvib = leftvibs[i[1]]
            rightvib = rightvibs[i[2]]
            coupl.vibop[i] = integrate(leftvib, coupl, rightvib)
        end
    end
end


"""Solve the vibrational problem for a given `Hamiltonian` and a list of vibrational
quantum number thresholds to discard eigenvectors above.
"""
function vibsolve!(ham::Hamiltonian; vmax_list::Vector{T}=fill(Inf, length(ham.potentials))) where {T<:Number}
    sep = diff(ham.grid[1:2])[1]
    for p=1:length(ham.potentials)
        pot = ham.potentials[p]
        vmax = vmax_list[p]
        eig = eigen(sincdvr(pot.values, diff(pot.nodes[1:2])[1], ham.mu)) #solve vibrational problem
        vmax == Inf && (vmax = length(eig.values)-1) #use all vibrational eigenstates
        vibbasis = Vector{Vib}(undef, Int(vmax)+1)
        for v=1:vmax+1
            vibbasis[v] = Vib(v-1, eig.values[v], eig.vectors[:,v])
        end
        ham.basisfrompot[pot] = vibbasis
    end
    return ham.basisfrompot
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
    @. vibronicdimen_list = (veemax_list + 1) * Int64(2*ess_list + 1)
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