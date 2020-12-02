module Juo
include("input.jl")
include("functions.jl")
include("schrödinger.jl")
using Plots, LinearAlgebra
using BenchmarkTools
#Read input file
#(fIn, fOut) = ReadInput.read_input()

#Define grid
rmin = 0.1
rmax = 10.1
dr = 0.01

_vibrational_grid = collect(rmin:dr:rmax)
npoints = length(_vibrational_grid)
println("No. initial grid points: $npoints")
_potential_grid = Potentials.morse.(_vibrational_grid; 
    D_e=0.1743608,
    r_e=1.401420894,
    mu= 918.5717371,
    a=1.440558
)
areless = isless.(_potential_grid, 0.173)
vibrational_grid = _vibrational_grid[areless]
potential_grid = _potential_grid[areless]
ngoodpoints = length(vibrational_grid)
println("No. retained grid points: $ngoodpoints")

potential_matrix = Diagonal(potential_grid)
kinetic_matrix = Vibrational.sinc_dvr_kinetic(vibrational_grid, dr, 918.5717371)
hamiltonian_matrix = potential_matrix + kinetic_matrix
vibrational_eigens = eigen(hamiltonian_matrix)

# println("Eigenvalues and Potentials:")
# for i =1:length(vibrational_grid)
#     println("E$i = ", vibrational_eigens.values[i])
#     println("v = ", findmax(vibrational_eigens.vectors[:, i]))
# end

plot(vibrational_grid, potential_grid)
savefig("potential.png")

heatmap(1:ngoodpoints, 1:ngoodpoints, hamiltonian_matrix)
savefig("hamiltonian.png")

heatmap(1:ngoodpoints, 1:ngoodpoints, vibrational_eigens.vectors)
savefig("eigenvectors.png")
end