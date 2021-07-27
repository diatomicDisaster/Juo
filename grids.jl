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
    quanta::CouplQuanta
    function UniformCoupling(nodes::Vector{T}, values::Vector{T}, quanta::CouplQuanta) where {T<:Number}
        length(nodes) != length(values) && throw(ArgumentError("Must have same number of nodes and values"))
        new{T}(nodes, values, quanta)
    end
end

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
    quanta::NamedTuple{(:lambda, :ess), Tuple{Float64, Float64}}
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
    interpolate(valgrid, vibgrid)

Interpolate the potential energy or coupling values `valgrid` to a new vibrational grid `vibgrid`.
"""
function interpolate(valgrid::T, vibgrid::AbstractVibGrid) where {T <: AbstractValGrid}
    on_grid(valgrid, vibgrid) && return T(valgrid.nodes, valgrid.values)
    spline = Spline1D(valgrid.nodes, curve.values)
    interpvals = spline(vibgrid.nodes)
    return T(vibgrid.nodes, interpvals)
end

function interpolate!(valgrid::AbstractValGrid, vibgrid::AbstractVibGrid)
    on_grid(valgrid, vibgrid) && return
    spline = Spline1D(valgrid.nodes, curve.values)
    interpvals = spline(vibgrid.nodes)
    valgrid.nodes = vibgrid.nodes
    valgrid.values = interpvals
end