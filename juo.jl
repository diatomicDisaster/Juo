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
    VibGrid(dr, nodes) = 
        !all(p -> p≈nodes[2]-nodes[1], diff(nodes)) ? error("Grid not uniformly spaced") :
        new(dr, nodes)
end
VibGrid() = VibGrid(0.05, collect(1.:.05:10.))
VibGrid(nodes::Vector{Float64}) = VibGrid(nodes[2]-nodes[1], nodes)

struct Coupling{N<:Number} <: AbstractGrid
    nodes::Vector{Float64}
    values::Vector{N}
    Coupling(dr, nodes, values) = 
        length(nodes) == length(values) ? error("Must have same number of nodes and values") : 
        new(dr, nodes, values)
end

struct Potential <: AbstractGrid
    nodes::Vector{Float64}
    values::Vector{Float64}
    uniform::Bool
    Potential(dr, nodes, values) = 
        length(nodes) == length(values) ? error("Must have same number of nodes and values") : 
        new(dr, nodes, values)
end

struct Diatom
    masses::Tuple{Float64, Float64}
    mu::Float64
    grid::VibGrid
    potential::Potential
    couplings::Vector{Coupling}
    vibbasis::Vector{VibState}
    Diatom(masses, mu) = new(masses, mu)
end
Diatom(masses::Tuple{Float64, Float64}) = Diatom(masses, prod(masses)/sum(masses))
Diatom(massA::Float64, massB::Float64) = Diatom((massA, massB))
Diatom(atomA::String, atomB::String) = Diatom(map(atom -> Data.atommass[atom], (atomA, atomB)))

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

function sincdvr(poten, dr, mu)
    npoints = length(poten)
    for j = 1:npoints
        for i = 1:j-1
            kinetic_operator[i, j] = (1/(mu*dr^2)) * (-1)^(i - j)/(i - j)^2
        end
        kinetic_operator[j, j] = (1/(2*mu*dr^2)) * (pi^2 / 3)
    end
    return eigen(Symmetric(kinetic_operator) + Diagonal(poten))
end

function integrate(vibbasis::Vector{VibState}, coupling::Coupling)
    for j = 
    for i, vibstate in enumerate(vibbasis)
        pre_vec = coupling.values .* vibstate.vector
        for f, vibstate in enumerate(vibbasis)
            mat[f, i] = pre_vec .* vibstate.vector
        end
    end
    return mat
end
        
function solve(diatom::Diatom, jaylist; nvee=Inf)
    vibeig = sincdvr(diatom.potential.values, diatom.vibgrid.dr, diatom.mu)
    nvee = nvee==Inf || nvee>length(vibeig.values) ? nvee : length(vibeig.values)
    for i = 1:nvee
        vee, val, vec, vee = i-1, vibeig.values[i], vibeig.vectors[i]
        diatom.vibbasis[i] = VibState(vee, val, vec)
    end
    
end

end



