using LinearAlgebra
using SparseArrays
using Dierckx


"""
    Hamiltonian(mu, grid)

Hamiltonian object for a molecule with reduced mass `mu` on the uniformly spaced
vibrational `grid` of internuclear distances.
"""
mutable struct Hamiltonian
    mu    :: Float64
    grid  :: Vector{Float64}
    potentials :: Vector{Potential}
    couplings  :: Vector{Coupling}
    vibfrompot   :: IdDict{Potential, Vector{Vib}}
    Hamiltonian(mu, grid, args...) =
        !all(p -> p≈diff(grid[1:2])[1], diff(grid)) ? 
        error("Grid not uniformly spaced") : 
        new(mu, grid, args...)
end

Hamiltonian(m::Float64, g::Vector{Float64}) = Hamiltonian(m, g, Vector{Potential}(), Vector{Potential}(), IdDict())

Main.push!(ham::Hamiltonian, poten::Potential) = ham(poten)
Main.push!(ham::Hamiltonian, coupl::Coupling) = ham(coupl)

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

function integrate!(ham::Hamiltonian)
    rsep = diff(ham.grid[1:2])[1]
    for coupl in ham.couplings
        leftvibs = ham.vibfrompot[coupl.left]
        rightvibs = ham.vibfrompot[coupl.right]
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
        pot, vmax = ham.potentials[p], vmax_list[p]
        eig = eigen(sincdvr(pot.values, diff(pot.nodes[1:2])[1], ham.mu)) #solve vibrational problem
        vmax == Inf && (vmax = length(eig.values)-1) #use all vibrational eigenstates
        vibbasis = Vector{Vib}(undef, Int(vmax)+1)
        for v=1:vmax+1
            vibbasis[v] = Vib(v-1, eig.values[v], eig.vectors[:,v])
        end
        ham.vibfrompot[pot] = vibbasis
    end
    return ham.vibfrompot
end

function spinsolve(basis::Vector{Ele})
    ndimen = length(basis)
    hamop  = Matrix{Float64}(undef, (ndimen, ndimen))
    for i in CartesianIndices(hamop)
        left, right = i[1], i[2]
        for 
end

function solve!(ham::Hamiltonian, jay::Float64)

    for poten in ham.potentials
        lambda, ess = poten.lambda, poten.ess
        for vib in ham.vibfrompot[poten]
            vee = vib.vee
            for sigma=-ess:ess
                omega = lambda + sigma
            end
        end
    end
end
