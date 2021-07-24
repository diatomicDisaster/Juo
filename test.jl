include("juo.jl")

coupl = Juo.UniformCoupling([1,2,3],[1,2,3])
println(coupl.coupling_quanta)

# x = collect(1.:1.:10.)
# y = x .^2
# o2 = Juo.Diatom(16.0, 16.0)
# pot = Juo.UniformPotential(x, y)
# o2(pot)
# println(o2)