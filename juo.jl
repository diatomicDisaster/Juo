module Juo
include("functions.jl")
include("data.jl")
include("grids.jl")
include("basis.jl")
using LinearAlgebra
using Dierckx
using .Grids, .Basis

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
    diatom.potential = poten
end

function (diatom::Diatom)(coupl::AbstractCoupling)
    push!(diatom.couplings, coupl)
end

# function (diatom::Diatom)(poten::T) where {T<:AbstractPotential}
#     if isdefined(diatom, :vibgrid)
#         interpolate!(poten, diatom.vibgrid)
#         diatom.potential = poten
#     else
#         throw(ErrorException("Diatom missing vibrational grid"))
#     end
    
# end

# function (diatom::Diatom)(coupl::AbstractCoupling)
#     if isdefined(diatom, :vibgrid)
#         interpolate!(coupl, diatom.vibgrid)
#         push!(diatom.couplings, coupl)
#     else
#         throw(ErrorException("Diatom missing vibrational grid"))
#     end
# end

# function (diatom::Diatom)(vibgrid::AbstractVibGrid)
#     if !isdefined(diatom, :vibgrid)
#         diatom.vibgrid = vibgrid
#     else
#         throw(ErrorException("Cannot add multiple vibrational grids"))
#     end
# end


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

end