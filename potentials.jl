module Potentials

function morse(r::Float64; r_e::Float64=1.0, a::Float64=1.0, D_e::Float64=1.0, T_e::Float64=0.0)
    return T_e + D_e*(1 - exp(-a*(r - r_e)))^2
end

end