module Vibronic

function grid(rmin::Float64, rmax::Float64, npoints::Int64)
    sep = (rmax - rmin)/npoints
    vib_grid :: Array{Float64, 1}(rmin:sep:rmax)
end

end