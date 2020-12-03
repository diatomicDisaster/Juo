module SchrodingerTest

include("functions.jl")
include("schrödinger.jl")

function test_vibrational(interval::Float64)
    (rmin, rmax) = (-10, 10)
    (r_e, w, m) = (0., 1., 1.)

    vib_grid = collect(rmin:interval:rmax)
    pot_grid = Potentials.qho.(vib_grid; r_e=r_e, w=w, m=m)
    (vib_eigen, vib_hamiltonian) = Vibrational.sincdvr(threshold=20., potential_grid=pot_grid, interval=interval, mass=m)
    
    exact_eigenvalues = w.*([i for i=0:length(vib_eigen.values)-1] .+ .5)
    diffs = vib_eigen.values - exact_eigenvalues
    return diffs
end

end