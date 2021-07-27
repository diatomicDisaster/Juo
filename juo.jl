module Juo
using LinearAlgebra
using SparseArrays
using Dierckx
include("data.jl")
include("functions.jl")
include("operators.jl")
include("grids.jl")
include("basis.jl")

mutable struct Diatom
    masses::Tuple{Float64, Float64}
    mu::Float64
    potential::AbstractPotential
    couplings::Vector{AbstractCoupling}
    Diatom(masses::Tuple{Float64, Float64}) = new(masses, prod(masses)/sum(masses))
end
Diatom(massA::Float64, massB::Float64) = Diatom((massA, massB))
Diatom(atomA::String, atomB::String) = Diatom(map(atom -> atommass[atom], (atomA, atomB)))

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


function vibrational_matrix(vibvecs::Matrix, coupling::Vector, dr::Float64)
    for (i, veci) in enumerate(eachcol(vibvecs))
        pre_vec = coupling .* veci
        for (f, vecf) in in enumerate(eachcol(vibvecs[:, 1:veci]))
            mat[f, i] = dr * vibbasis.grid.dr * vecf⋅pre_vec
        end
    end
    return Hermitian(mat)
end

vibrational_matrix(vibbasis::VibBasis, coupling::AbstractCoupling, dr::Float64) = 
    vibrational_matrix(vibbasis.eigen.vectors, coupling.values, dr)

function vibronic_matrix(elebasis, vibbasis, coupls)
    for coupl in coupls
        elemat = zeros(length(elebasis), length(elebasis))
        f = findall(e->e=(coupl.quanta[1:3]), elebasis)
        i = findall(e->e=(coupl.quanta[4:6]), elebasis)
        elemat[f, i] = 1.0
        vibmat = vibrational_matrix(vibbasis, coupl, dr)
        vibronic_matrix += kron(elemat, vibmat)
    end
    return vibronic_matrix
end

function solve(jay, mu, dr, vibgrid, poten, poten_quanta, coupls, coupl_quantas)
    elebasis = Juo.EleBasis(poten_quanta...)
    vibbasis = Juo.VibBasis(mu, dr, poten; nvee=5)
    rotbasis = Juo.RotBasis(jay, elebasis)
    rovibronicbasis = Juo.RovibronicBasis(elebasis, vibbasis, rotbasis)
    coupl_mats = vibrational_matrix.(vibbasis, coupls, dr)
    for statei in rovibronicbasis
        for statef in rovibronicbasis
            coupl_mels = map(c -> c[statef.vee, statei.vee], coupl_mats)

        end
    end
end

end