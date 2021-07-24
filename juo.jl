module Juo
include("functions.jl")

using LinearAlgebra
using Dierckx
include("data.jl")

abstract type AbstractGrid end
Base.length(G::AbstractGrid) = length(G.nodes)
abstract type AbstractVibGrid <: AbstractGrid end
abstract type AbstractValGrid <: AbstractGrid end
abstract type AbstractCoupling{T<:Number} <: AbstractValGrid end
abstract type AbstractPotential{T<:Number} <: AbstractValGrid end

"""
    UniformVibGrid <: AbstractVibGrid

Uniform vibrational grid of internuclear distances `nodes` with separation `dr`.

# Examples
```jldoctest
julia> Juo.UniformVibGrid(1.:.5:5.)
Main.Juo.UniformVibGrid([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0], 0.5)
```
"""
struct UniformVibGrid <: AbstractVibGrid
    nodes::Vector{Float64}
    dr::Float64
    UniformVibGrid(nodes) = 
        !all(p -> p≈nodes[2]-nodes[1], diff(nodes)) ? throw(ArgumentError("Nodes not uniformly spaced")) :
        new(nodes, nodes[2]-nodes[1])
end
UniformVibGrid() = UniformVibGrid(1.:.05:10.)

"""
    UniformCoupling{T<:Number} <: AbstractCoupling

Coupling curve of `values` on a grid of uniformly spaced `nodes`, of type `T`. 
The curve couples electronic states with final and initial quantum numbers 
``Λ_f =`` `lambdaf`, ``S_f =`` `essf`, ``Σ_f =`` `sigmaf` and 
``Λ_i =`` `lambdai`, ``S_i =`` `essi`, ``Σ_i =`` `sigmai`, respectively.  

# Examples
```jldoctest
julia> Juo.UniformCoupling(collect(1:.5:5), collect(1:.5:5).^2)
Main.Juo.UniformCoupling{Float64}([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0], [1.0, 2.25, 4.0, 6.25, 9.0, 12.25, 16.0, 20.25, 25.0], (lambdaf = NaN, essf = NaN, sigmaf = NaN, lambdai = NaN, essi = NaN, sigmai = NaN))
```
"""
struct UniformCoupling{T} <: AbstractCoupling{T}
    nodes::Vector{T}
    values::Vector{T}
    coupling_quanta::NamedTuple{(:lambdaf, :essf, :sigmaf, :lambdai, :essi, :sigmai), Tuple{Float64, Float64, Float64, Float64, Float64, Float64}}
    function UniformCoupling(nodes, values, coupling_quanta)
        length(nodes) != length(values) && throw(ArgumentError("Must have same number of nodes and values"))
        new{T}(nodes, values)
    end
end
UniformCoupling(nodes, values) = UniformCoupling(nodes, values,
    (lambdaf=NaN, essf=NaN, sigmaf=NaN, lambdai=NaN, essi=NaN, sigmai=NaN)
)


"""
    UniformPotential{T<:Number} <: AbstractPotential

Potential energy curve of `values` on a grid of uniformly spaced `nodes`, of type `T`.
Represents an electronic state with quantum numbers ``Λ =`` `lambda`, ``S =`` `ess`, ``Σ =`` `sigma`.

# Examples
```jldoctest
julia> Juo.UniformPotential(collect(1:.5:5), collect(1:.5:5).^2)
Main.Juo.UniformPotential{Float64}([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0], [1.0, 2.25, 4.0, 6.25, 9.0, 12.25, 16.0, 20.25, 25.0])
```
"""
struct UniformPotential{T} <: AbstractPotential{T}
    nodes::Vector{T}
    values::Vector{T}
    function UniformPotential(nodes::Vector{T}, values::Vector{T}) where {T<:Number}
        length(nodes) != length(values) && throw(ArgumentError("Must have same number of nodes and values"))
        new{T}(nodes, values)
    end
end

"""
    on_grid(valgrid, vibgrid)

Check if the potential energy or coupling values `valgrid` are on the points `vibgrid`.
""" 
@inline on_grid(valgrid::AbstractValGrid, vibgrid::AbstractVibGrid) = valgrid.nodes ≈ vibgrid.nodes

"""
    interpolate!(valgrid, vibgrid)

Interpolate the potential energy or coupling values `valgrid` to a new vibrational grid `vibgrid`.
"""
function interpolate!(valgrid::AbstractValGrid, vibgrid::AbstractVibGrid)
    on_grid(valgrid, vibgrid) && return
    spline = Spline1D(valgrid.nodes, curve.values)
    interpvals = spline(vibgrid.nodes)
    valgrid.nodes = vibgrid.nodes
    valgrid.values = interpvals
end

mutable struct Diatom
    masses::Tuple{Float64, Float64}
    mu::Float64
    potential::AbstractPotential
    couplings::Vector{AbstractCoupling}
    Diatom(masses::Tuple{Float64, Float64}) = new(masses, prod(masses)/sum(masses))
end
Diatom(massA::Float64, massB::Float64) = Diatom((massA, massB))
Diatom(atomA::String, atomB::String) = Diatom(map(atom -> Data.atommass[atom], (atomA, atomB)))

function (diatom::Diatom)(poten::T) where {T<:AbstractPotential}
    if isdefined(diatom, :vibgrid)
        interpolate!(poten, diatom.vibgrid)
        diatom.potential = poten
    else
        throw(ErrorException("Diatom missing vibrational grid"))
    end
    
end

function (diatom::Diatom)(coupl::AbstractCoupling)
    if isdefined(diatom, :vibgrid)
        interpolate!(coupl, diatom.vibgrid)
        push!(diatom.couplings, coupl)
    else
        throw(ErrorException("Diatom missing vibrational grid"))
    end
end

function (diatom::Diatom)(vibgrid::AbstractVibGrid)
    if !isdefined(diatom, :vibgrid)
        diatom.vibgrid = vibgrid
    else
        throw(ErrorException("Cannot add multiple vibrational grids"))
    end
end

function sincdvr(mu, dr, poten)
    npoints = length(poten)
    for j = 1:npoints
        for i = 1:j-1
            kinetic_operator[i, j] = (1/(mu*dr^2)) * (-1)^(i - j)/(i - j)^2
        end
        kinetic_operator[j, j] = (1/(2*mu*dr^2)) * (pi^2 / 3)
    end
    return Symmetric(kinetic_operator) + Diagonal(poten)
end

struct VibBasis
    eigen::Eigen{Float64, Float64, Matrix{Float64}, Vector{Float64}}
    quanta::Vector{NamedTuple{(:vee,), Tuple{Int64}}}
end

function VibBasis(mu, dr, poten; nvee=nothing)
    nvee = isnothing(nvee) || nvee>length(poten) ? nvee : length(poten)
    vibham = sincdvr(poten, dr, mu)
    vibeig = eigen(vibham, 1:nvee)
    for i = 1:length(vibeig.values)
        vee = i-1
        quanta[i] = (vee,)
    end
    vibbasis = VibBasis(vibeig, quanta)
    return vibbasis
end

VibBasis(mu, dr, poten::AbstractPotential; kwargs...) = VibBasis(mu, dr, poten.values; kwargs...)

function coupling_matrix(vibvecs::Matrix, coupling::Vector, dr::Float64)
    for (i, veci) in enumerate(eachcol(vibvecs))
        pre_vec = coupling .* veci
        for (f, vecf) in in enumerate(eachcol(vibvecs[:, 1:veci]))
            mat[f, i] = dr * vibbasis.grid.dr * vecf⋅pre_vec
        end
    end
    return Hermitian(mat)
end

coupling_matrix(vibbasis::VibBasis, coupling::AbstractCoupling, dr::Float64) = 
    coupling_matrix(vibbasis.eigen.vectors, coupling.values, dr)

struct EleBasis
    quanta::Vector{NamedTuple{(:lambda, :ess, :sigma), Tuple{Float64, Float64, Float64}}}
    multiplicity::Float64
    term::Char
end

@inline function lambda_to_char(lambda::Real)
    char = "ΣΠΔΦ"[abs(round(Int, lambda))-1]
end

function EleBasis(lambda, ess, gerade; sym=missing)
    multi = round(Int, 2*ess + 1)
    nlam = round(Int, 2*lambda)
    quanta = Vector{NamedTuple}(undef, multi*nlam)
    i = 1
    for lam in (lambda==0. ? [lambda] : [-abs(lambda), abs(lambda)])
        for sigma = -ess:ess
            quanta[i] = (lambda=lam, ess=ess, sigma=sigma)
            i += 1
        end
    end
    #construct term symbol
    gu = gerade ? "g" : "u"
    if ismissing(sym)
        term = string(multi, lambda_to_char(lambda), gu)
    else
        pm = sym ? "+" : "-"
        term = string(multi, lambda_to_char(lambda), gu, pm)
    end
    return EleBasis(quanta, multi, term)
end

struct RotBasis
    quanta::Vector{NamedTuple{(:jay, :omega), Tuple{Float64, Float64}}}
end

function RotBasis(jay, elebasis)
    quanta = []
    for qnums in elebasis.quanta
        (lambda, ess, sigma) = qnums
        omega = lambda + sigma
        if abs(omega) > jay
            continue
        else
            push!(quanta, (jay=jay, omega=omega))
        end
    end
    return RotBasis(quanta)
end 

end

