struct Vibronic
    state  :: Int64
    energy :: Float64
    vector :: Array{Float64}
end

struct Rovibronic
    state  :: Int64
    vee    :: Int64
    lambda :: Float64
    ess    :: Float64
    sigma  :: Float64
    jay    :: Float64
    omega  :: Float64
end