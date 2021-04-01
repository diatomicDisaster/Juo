module SchrodingerTest

include("functions.jl")
include("schrödinger.jl")

using Plots, LinearAlgebra
using BenchmarkTools

function test_vibrational()
    #Morse potential parameters for CO
    # A. Roy (2013) - https://doi.org/10.1016/j.rinp.2013.06.001
    D_e = 0.412533191
    mu  = 12506.23981
    r_e = 2.132177990 
    a   = 2.59441

    #Define grid (same as original Duo paper)
    # S. Yurchenko et al. (2016) - https://doi.org/10.1016/j.cpc.2015.12.021
    rmin = 1.322808287
    rmax = 3.77845225
    dr = 0.01
    thresh = 50.
    show_levels = 30

    vibrational_grid = collect(rmin:dr:rmax)
    npoints = length(vibrational_grid)
    println("No. initial grid points: $npoints")
    potential_grid = Potentials.morse.(vibrational_grid; r_e=r_e, a=a, D_e=D_e)

    (vibrational_eigen, vibrational_hamiltonian) = Vibrational.sincdvr(potential_grid=potential_grid, interval=dr, threshold=thresh, mass=mu)
    ngoodpoints = length(vibrational_eigen.values)
    println("No. retained grid points: $ngoodpoints")

    plot(vibrational_grid, potential_grid)
    savefig("potential.png")

    heatmap(1:ngoodpoints, 1:ngoodpoints, vibrational_hamiltonian)
    savefig("hamiltonian.png")

    heatmap(1:ngoodpoints, 1:show_levels, transpose(vibrational_eigen.vectors[:,1:show_levels].^2))
    savefig("eigenvectors.png")
end
end