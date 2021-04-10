struct RovibronicBasis
    xi     :: Int64
    lambda :: Int64
    vee    :: Int64
    spin   :: Int64
    sigma  :: Int64
    jay    :: Float64
    em     :: Float64
    omega  :: Float64
    vibrational :: Vector{Float64}
end

struct VibrationalGrid
    rmin :: Float64
    rmax :: Float64
    sep  :: Float64
    num  :: Int64
    vals :: Vector{Float64}
end

struct BasisSet
    nstates :: Int64
    veemaxlist :: Vector{Int64}
    lambdalist :: Vector{Float64}
    esslist    :: Vector{Float64}
    jaylist    :: Vector{Float64}
    ndimen     :: Int64
    rotdimen  :: Int64
    vibdimen  :: Int64
    elecdimen :: Int64
end