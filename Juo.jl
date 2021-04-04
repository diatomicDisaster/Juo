module Juo

using BenchmarkTools
using DataStructures

include("Input.jl")
include("Functions.jl")
include("Schrödinger.jl")
include("test.jl")

Jlist = [i for i in 0:10.0]


ndimensions = dimensions(Jlist)

Rotational.raising_operator(Jlist, ndimensions)

end