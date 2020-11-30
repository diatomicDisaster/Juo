module Juo
include("Input.jl")
include("Functions.jl")
include("Schrödinger.jl")
using Plots
#Read input file
(fIn, fOut) = ReadInput.read_input()

#Define grid
const N = 501
vib_grid = collect(LinRange(0.1, 5.0, N))
@time pot_grid = PotentialCurves.mor_pot(vib_grid; a=1.0, r_e=1.0, D_e=1.0)
@time kinetic_energy_mat = Vibrational.sinc_dvr(N)
heatmap(1:N, 1:N, kinetic_energy_mat)
plot(vib_grid, pot_grid)
savefig("fig.png")
end