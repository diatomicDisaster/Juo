include("Juo.jl")

using CSV, DataFrames

oxygen_data = CSV.read(joinpath("files", "oxygen.csv"), DataFrame)

potgrid = Vector(.5:.1:11.)
X3Σ = Juo.Potentials.morse.(potgrid; D_e=0.1916182, T_e=0.0, r_e=2.281844296, a=1.471112646)
a1Δ = Juo.Potentials.morse.(potgrid; D_e=0.1555084, T_e=0.036077564, r_e=2.29734005, a=1.497571507)
b1Σ = Juo.Potentials.morse.(potgrid; D_e=0.108659, T_e=0.0596876277, r_e=2.3164263, a=1.5610728)
potens = [
    #Juo.Potential(potgrid, b1Σ, 0., 0.),
    #Juo.Potential(potgrid, a1Δ, 2., 0.),
    Juo.Potential(potgrid, X3Σ, 0., 1.)
]


vibgrid = Vector(.8:.1:10.)
hamil = Juo.Hamiltonian(1.0, vibgrid)
hampotens = hamil.(potens)

couples = [
    Juo.Coupling(oxygen_data.R, oxygen_data.QX3Σ, hampotens[1], hampotens[1])
]

hamil.(couples)

Juo.vibsolve!(hamil)
Juo.integrate!(hamil)
