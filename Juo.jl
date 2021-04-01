module Juo
using BenchmarkTools
include("input.jl")
include("functions.jl")
include("schrödinger.jl")
include("test.jl")

#SchrodingerTest.test_vibrational()
println("Working")
@benchmark Vibrational.kinetic_operator(100, 0.1)
@benchmark Rotational.raising_operator(Jmax=50.0)

end