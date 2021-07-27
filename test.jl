include("juo.jl")

diat = Juo.Diatom("O", "O")
dr = .05
rs = collect(1:dr:10)
morse_pot = Juo.Potentials.morse.(rs)
poten = Juo.UniformPotential(rs, morse_pot)
coupl_quanta = (lambdaf=0, essf=1, sigmaf=0, lambdai=0, essi=1, sigmai=1)
quad_coupl = Juo.UniformCoupling(rs, 0.1.*rs .- 5, coupl_quanta)

elebasis = Juo.EleBasis(1, 0, true, false)
vibbasis = Juo.VibBasis(diat.mu, dr, poten; nvee=5)
rotbasis = Juo.RotBasis(0.0, elebasis)
rovibronicbasis = Juo.RovibronicBasis(elebasis, vibbasis, rotbasis)

println("Done.")