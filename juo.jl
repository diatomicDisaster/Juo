module Juo

using LinearAlgebra

struct VibState
    vee :: Int64
    energy :: Float64
    vector :: Vector{Float64}
end

struct EleState
    ele :: Int64
    lambda :: Int64
    ess :: Float64
    sigma :: Float64
end

function sincdvr(potential_grid::Vector, dr, mu)
    npoints = length(potential_grid)
    kinetic_operator = Array{Float64}(undef, (npoints, npoints))
    for j = 1:npoints
        for i = 1:j-1
            kinetic_operator[i, j] = (1 / (2*mu*dr^2)) * 2*( ((-1)^(i - j)) / (i - j)^2 )
        end
        kinetic_operator[j, j] = (1 / (2*mu*dr^2)) * (pi^2 / 3)
    end
    eig = eigen(Symmetric(kinetic_operator) + Diagonal(potential_grid))
    nvib = length(eig.values)
    vibeigen = Vector{VibState}(undef, nvib)
    for i = 0:nvib-1
        vibeigen[i] = VibState(i-1, val, vec)
    end
    return vibeigen
end

end