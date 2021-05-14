include("types.jl")

function buildbasis(jay, poten::Potential, vib::Array{Vibrational})
    for sigma=-poten.ess:poten.ess
        ele = Electronic(poten, poten.ess, sigma, poten.lambda)
        for 