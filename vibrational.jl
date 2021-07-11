module Vibrational

function sincdvr(potential_grid::Vector, dr, mu)
    npoints = length(potential_grid)
    kinetic_operator = Array{Float64}(undef, (npoints, npoints))
    for j = 1:npoints
        for i = 1:j-1
            kinetic_operator[i, j] = (1 / (2*mu*dr^2)) * 2*( ((-1)^(i - j)) / (i - j)^2 )
        end
        kinetic_operator[j, j] = (1 / (2*mu*dr^2)) * (pi^2 / 3)
    end
    ham = Symmetric(kinetic_operator) + Diagonal(potential_grid)
    eig = eigen(ham)
    return eig
end

end