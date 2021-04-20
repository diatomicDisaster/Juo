
include("Juo.jl")

using LinearAlgebra
using CSV, DataFrames

function full()
    println("Testing Juo")

    #Build potential energy grids
    reducedmass = 14583.10789
    (rmin, rmax) = (1.5, 7.5)
    gridsep = 0.5
    vibgrid = [i for i=rmin:gridsep:rmax]

    oxygen_data = CSV.read(joinpath("files", "oxygen.csv"), DataFrame)
    couplings = [Juo.Interpolate(Juo.Grid(oxygen_data.R, oxygen_data.Q), vibgrid)]

    print("Building potential grids... ")
    potens  = [#Morse potential parameters for O2
        Juo.Grid(vibgrid,
            Juo.Potentials.morse.(vibgrid; D_e=0.1916182, T_e=0.0, r_e=2.281844296, a=1.471112646)
        ),
        Juo.Grid(vibgrid,
            Juo.Potentials.morse.(vibgrid; D_e=0.1555084, T_e=0.036077564, r_e=2.29734005, a=1.497571507)
        )
    ]
    println("Done.")
    
    #Solve vibrational problem
    veemax_s = [3, 3]
    print("Building vibronic basis... ")
    vibronic_basis = Juo.vibronic_basis.(veemax_s, [p.values for p in potens], gridsep, reducedmass)
    println("Done.")
    
    lambda_s = [0., 2.]
    ess_s    = [1., 0.]
    jay_list = [i for i=0.:10.]
    for jay in jay_list
        #rovibronic_basis = Juo.rovibronic_basis(jay, veemax_s, lambda_s, ess_s)
        #calculate rovibronic hamiltonian
    end
end

full()
