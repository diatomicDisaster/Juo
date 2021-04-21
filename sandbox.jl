include("Juo.jl")
include("Functions.jl")

vibgrid = Vector(.8:.1:10.)
hamil = Juo.Hamiltonian(1.0, vibgrid)

potgrid = Vector(.5:.2:11.)
X3Sigma = Potentials.morse.(potgrid; D_e=0.1916182, T_e=0.0, r_e=2.281844296, a=1.471112646)
poten1 = Juo.Potential(potgrid, X3Sigma, 0., 0.)

hamil(poten1)
Juo.vibsolve!(hamil)

ham