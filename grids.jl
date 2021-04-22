abstract type AbstractGrid end
Base.length(A::AbstractGrid) = length(A.nodes)

struct Grid <: AbstractGrid
    nodes   :: Vector{Float64}
    values  :: Vector{Float64}
end

struct Potential <: AbstractGrid
    nodes   :: Vector{Float64}
    values  :: Vector{Float64}
    lambda  :: Float64
    ess     :: Float64
end

struct Coupling <: AbstractGrid
    nodes   :: Vector{Float64}
    values  :: Vector{Float64}
    left    :: Potential
    right   :: Potential
    vibop   :: Matrix{Float64}
end

Coupling(nodes, values, left, right) = 
    Coupling(nodes, values, left, right, Matrix{Float64}(undef, length(nodes), length(nodes)))
#Base.iterate(coupl::Coupling, state=1) = state > length(coupl.nodes) ? nothing : (coupl.nodes[state], state+1)
