module Juo
include("functions.jl")

using LinearAlgebra
using Dierckx
include("data.jl")

abstract type AbstractGrid end
Base.length(G::AbstractGrid) = length(G.nodes)

struct VibGrid <: AbstractGrid
    dr::Float64
    nodes::Vector{Float64}
end

struct Coupling <: AbstractGrid
    dr::Float64
    nodes::Vector{Float64}
    values::Vector{Float64}
    Coupling(dr, nodes, values) = length(nodes) == length(values) ? error("Must have same number of nodes and values") : new(dr, nodes, values)
end

struct Potential <: AbstractGrid
    dr::Float64
    nodes::Vector{Float64}
    values::Vector{Float64}
    Potential(dr, nodes, values) = length(nodes) == length(values) ? error("Must have same number of nodes and values") : new(dr, nodes, values)
end

struct Diatom
    masses::Tuple{Float64, Float64}
    mu::Float64
    grid::VibGrid
    potential::Potential
    couplings::Vector{Coupling}
end
Diatom(masses::Tuple{Real, Real}) = Diatom(masses, prod(masses)/sum(masses))
Diatom(massA::Real, massB::Real) = Diatom((massA, massB))
Diatom(atoms::Tuple{String, String}) = Diatom(map(atom -> Data.atommass[atom], atoms))
Diatom(atomA::String, atomB::String) = Diatom((atomA, atomB))

function (diatom::Diatom)(poten::Potential)
    spline = spline1D(poten.nodes, poten.values)
    interpvals = spline(diatom.grid.nodes)
    interppoten = Potential(diatom.grid.dr, diatom.grid.nodes, interpvals)
    diatom.potential = interppoten
    return interppoten
end

function (diatom::Diatom)(coupl::Coupling)
    spline = spline1D(coupl.nodes, coupl.values)
    interpvals = spline(diatom.grid.nodes)
    interpcoupl = Potential(diatom.grid.dr, diatom.grid.nodes, interpvals)
    push!(diatom.couplings, interpcoupl)
    return interpcoupl
end

struct VibState
    vee::Int64
    energy::Float64
    vector::Vector{Float64}
end

function sincdvr(potential::Vector{Real}, dr, mu)
    npoints = length(potential_grid)
    kinetic_operator = Array{Real}(undef, (npoints, npoints))
    for j = 1:npoints
        for i = 1:j-1
            kinetic_operator[i, j] = (1/(mu*dr^2)) * (-1)^(i - j)/(i - j)^2
        end
        kinetic_operator[j, j] = (1/(2*mu*dr^2)) * (pi^2 / 3)
    end
    eig = eigen(Symmetric(kinetic_operator) + Diagonal(potential_grid))
    nvib = length(eig.values)
    vibbasis = Vector{VibState}(undef, nvib)
    for i = 0:nvib-1
        vibbasis[i] = VibState(i-1, val, vec)
    end
    return vibbasis
end

end



