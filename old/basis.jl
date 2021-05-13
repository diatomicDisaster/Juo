abstract type AbstractBasis end
Base.iterate(x::AbstractBasis) = (x, nothing)
Base.iterate(x::AbstractBasis, ::Any) = nothing
Base.isempty(x::AbstractBasis) = false
in(x::AbstractBasis, y::AbstractBasis) = x == y

"""Vibrational basis state for a given vibrational quantum number `vee` represented
by a `vector` on the `Hamiltonian` `grid`, with eigenvalue `energy`.
"""
struct Vib <: AbstractBasis
    vee    :: Int64
    energy :: Float64
    vector :: Vector{Float64}
end

"""Electronic basis state for a given `state` with electronic angular momentum 
projection `lambda`, spin `ess` and vibrational basis `vib`.
"""
struct Ele <: AbstractBasis
    state   :: Int64
    lambda  :: Float64
    ess     :: Float64
    sigma   :: Float64
end

struct Vibronic <: AbstractBasis
    ele :: Ele
    vee :: Int64
end

struct Rot <: AbstractBasis
    jay   :: Float64
    omega :: Float64
end

struct Rovibronic <: AbstractBasis
    rot :: Rot
    vibronic :: Vibronic
end

function Rovibronic(
    state::Int64, lambda::Float64, 
    ess::Float64, sigma::Float64, 
    vee::Int64, 
    jay::Float64, omega::Float64)
    
    vib = Vib(vee)
end


"""
Build the rovibronic basis for a given total angular momentum quantum number and list of electronic quantum numbers.
"""
function rovibronic_basis(
    jay         :: Vector{Float64},
    veemax_list :: Vector{Int64},
    lambda_list :: Vector{Float64},
    ess_list    :: Vector{Float64}
)
    length(veemax_list) == length(lambda_list) == length(ess_list) || throw(ArgumentError("Length of veemax_list (=$(length(veemax_list))), lambda_list (=$(length(lambda_list))), and ess_list (=$(length(ess_list))) should be equal."))
    vibronicdimen_list = Array{Int64}(undef, length(veemax_list))
    @. vibronicdimen_list = (veemax_list + 1) * Int64(2*ess_list + 1)
    basis = Vector{Rovibronic}(undef, sum(vibronicdimen_list))
    for state = 1:nstates
        lambda = lambda_list[state]
        ess    = ess_list[state]
        veemax = veemax_list[state]
        for (s, sigma) in enumerate(-ess:ess)
            omega = lambda + sigma
            for vee = 0:veemax
                i = (j - 1) * (sum(vibronicdimen_list)) + sum(vibronicdimen_list[1:(state-1)]) + (s-1)*(veemax+1) + (vee+1)
                basis[i] = Rovibronic(state, vee, lambda, ess, sigma, jay, omega)
            end
        end
    end
    return basis
end