module SchrodingerTest

include("Functions.jl")
include("Schrödinger.jl")

using Plots, LinearAlgebra
using BenchmarkTools
using Arpack
using ColorSchemes

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

    vib_grid = collect(rmin:dr:rmax)
    npoints = length(vib_grid)
    println("No. initial grid points: $npoints")
    pot_grid = Potentials.morse.(vib_grid; r_e=r_e, a=a, D_e=D_e)

    vib_hamiltonian = Vibrational.hamiltonian(pot_grid, dr, mu)
    vibrational_eigen = eigen(vib_hamiltonian)
    ngoodpoints = length(vibrational_eigen.values[1:show_levels])
    println("No. retained grid points: $ngoodpoints")

    plot(vibrational_grid, potential_grid)
    savefig(joinpath("images", "potential.png"))

    heatmap(1:ngoodpoints, 1:ngoodpoints, vibrational_hamiltonian)
    savefig(joinpath("images", "hamiltonian.png"))

    heatmap(1:ngoodpoints, 1:show_levels, transpose(vibrational_eigen.vectors[:,1:show_levels].^2))
    savefig(joinpath("images", "eigenvectors.png"))
end

function test_rotational()
    Js = [i for i=0.:20.]
    rot_hamiltonian = Rotational.hamiltonian(Js, 0.5, 2.0, 1.0)
    rot_eigen = Rotational.solve(rot_hamiltonian)
    xs = [0, 1]
    p = plot(xs, [0., 0.])
    for val in real(rot_eigen.values)
        plot!(p, xs, [val, val])
    end
    savefig(joinpath("images", "rotational.png"))
end

function test_rovibrational()
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
    vmax = 2

    ndimen_vib = vmax + 1
    vib_grid = collect(rmin:dr:rmax)
    npoints = length(vib_grid)
    println("No. vibrational points: $npoints")
    pot_grid = Potentials.morse.(vib_grid; r_e=r_e, a=a, D_e=D_e)

    vib_ham = Vibrational.sincdvr(pot_grid, dr, mu)
    vib_eigen = eigen(vib_ham)
    vib_hamiltonian = Diagonal(vib_eigen.values[1:ndimen_vib])

    Js = [i for i=0.:5.]
    rot_ham = Rotational.hamiltonian(Js, 0., mu, r_e)
    rot_hamiltonian = Matrix(rot_ham)
    ndimen_rot = size(rot_hamiltonian)[1]

    identity = Matrix{Float64}(I, (ndimen_rot, ndimen_rot))
    vibrational_hamiltonian = kron(vib_hamiltonian, identity)
    identity = Matrix{Float64}(I, (ndimen_vib, ndimen_vib))
    rotational_hamiltonian  = kron(identity, rot_hamiltonian)

    rovib_hamiltonian = vibrational_hamiltonian + rotational_hamiltonian
    rovib_eigen = eigen(rovib_hamiltonian)
    
    reals = real(rovib_eigen.values)
    xs = [0, 1]
    p = plot(xs, [reals[1], reals[1]])
    for (ind, val) in enumerate(reals)
        if ind != 1
            plot!(p, xs, [val, val])
        end
    end
    savefig(joinpath("images", "rovibrational.pdf"))
    plot!(p, ylims=[0.01, 0.011])
    savefig(joinpath("images", "rovibrational_v0.pdf"))

end

end