module Juo

include("input.jl")
include("functions.jl")
include("schrödinger.jl")
include("test.jl")

using Plots, LinearAlgebra
using BenchmarkTools

intervals = [1.0, 0.75, 0.5, 0.25, 0.1, 0.075, 0.05, 0.025, 0.01]
diffs_list = SchrodingerTest.test_vibrational.(intervals)
for v in [1, 6, 11]
    plot!(intervals, abs.([diffs_list[i][v] for i=1:length(intervals)]), yaxis=:log, xaxis=:log)
end
savefig("vibtest.png")

#Define grid
rmin = -10
rmax = 10
dr = 0.1
thresh = 50.
show_levels = 30

vibrational_grid = collect(rmin:dr:rmax)
npoints = length(vibrational_grid)
println("No. initial grid points: $npoints")
potential_grid = Potentials.qho.(vibrational_grid)

(vibrational_eigen, vibrational_hamiltonian) = Vibrational.sincdvr(potential_grid=potential_grid, interval=dr, threshold=thresh, mass=1.0)
ngoodpoints = length(vibrational_eigen.values)
println("No. retained grid points: $ngoodpoints")

plot(vibrational_grid, potential_grid)
savefig("potential.png")

heatmap(1:ngoodpoints, 1:ngoodpoints, vibrational_hamiltonian)
savefig("hamiltonian.png")

heatmap(1:ngoodpoints, 1:show_levels, transpose(vibrational_eigen.vectors[:,1:show_levels].^2))
savefig("eigenvectors.png")
end