include("Juo.jl")
include("functions.jl")
using LinearAlgebra
using Plots
ENV["GKSwstype"]="nul"

geoms = collect(1.:0.02:10.0)

poten1 = Potentials.morse.(geoms; D_e=0.1916182, T_e=0.0, r_e=2.281844296, a=1.471112646)
poten2 = Potentials.morse.(geoms; D_e=0.1555084, T_e=0.036077564, r_e=2.29734005, a=1.497571507)

potentials = [
    Potential("1", poten1, geoms, 0, 1, false, true),
    Potential("2", poten2, geoms, 0, 0, true, true)
]

vibbasis = sincdvr_basis(29166.2, potentials)

println(length(geoms))
for potential in potentials
    println("New poten")
    basis = vibbasis[map(v -> v.potential, vibbasis).==potential]
    basis = basis[1:100]
    geometries = basis[1].potential.geometries
    levels = map(v -> v.energy, basis)
    matrix = Matrix(undef, length(levels), length(geometries))
    for i = 1:size(levels, 1)
        matrix[i,:] = basis[i].vec
    end
    savefig(heatmap(geometries, levels, matrix.^2, reuse=false, ylims=[0, last(basis).energy]), "Potential $(potential.label)")
end
