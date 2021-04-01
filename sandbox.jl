using BenchmarkTools


function vect(k)
    return sum(k .* 2 .^(UnitRange(1, size(k)[1]) .+ 3))
end

function loop(k)
    x = 0.
    for (i, j) in enumerate(k)
        x += k[i]*2^(j+3)
    end
    return x
end

k = [.5, .6, .7, .4, .3, .2, .6, .7, .4, .3, .2, .6, .7, .4, .3, .2, .6, .7, .4, .3, .2]

@time vect(k)
@time loop.(k)
@time sum(k .* 2 .^(UnitRange(1, size(k)[1]) .+ 3))
@btime vect($k)
@btime loop($k)
@btime sum(k .* 2 .^(UnitRange(1, size(k)[1]) .+ 3))