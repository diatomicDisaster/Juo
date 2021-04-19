module JuoTest

include("input.jl")
include("functions.jl")
include("schrödinger.jl")
include("grids.jl")

using LinearAlgebra
using CSV, DataFrames

function full()
    println("Testing Juo")

    #Build potential energy grids
    reducedmass = 14583.10789
    (rmin, rmax) = (1.5, 7.5)
    gridsep = 0.5
    vibgrid = [i for i=rmin:gridsep:rmax]

    oxygen_data = CSV.read("oxygen.csv", DataFrame)
    couplings = [Grids.Interpolate(Grids.Grid(oxygen_data.R, oxygen_data."Q11"), vibgrid)]

    print("Building potential grids... ")
    potens  = [#Morse potential parameters for O2
        Grid(vibgrid,
            Potentials.morse.(vibgrid.values; D_e=0.1916182, T_e=0.0, r_e=2.281844296, a=1.471112646)
        ),
        Grid(vibgrid,
            Potentials.morse.(vibgrid.values; D_e=0.1555084, T_e=0.036077564, r_e=2.29734005, a=1.497571507)
        )
    ]
    println("Done.")
    
    #Solve vibrational problem
    veemax_s = [30, 30]
    print("Building vibronic basis... ")
    vibronic_basis = Hamiltonian.vibronic_basis(potens, veemax_s, gridsep, reducedmass)
    println("Done.")
    
    lambda_s = [0., 2.]
    ess_s    = [1., 0.]
    jay_list = [i for i=0.:10.]
    for jay in jay_list
        rovibronic_basis = Hamiltonian.rovibronic_basis(jay, veemax_s, lambda_s, ess_s)
        #calculate rovibronic hamiltonian
    end
end

end