module Test

include("input.jl")
include("functions.jl")
include("schrödinger.jl")

using LinearAlgebra

function full()
    println("Testing Juo")

    #Build potential energy grids
    reducedmass = 14583.10789
    gridsep = 0.05
    Rvals   = [i for i=1.5:gridsep:7.5]
    print("Building potential grids... ")
    potens  = [#Morse potential parameters for O2
        Potentials.morse.(Rvals; D_e=0.1916182, T_e=0.0, r_e=2.281844296, a=1.471112646),
        Potentials.morse.(Rvals; D_e=0.1555084, T_e=0.036077564, r_e=2.29734005, a=1.497571507)
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