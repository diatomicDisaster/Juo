abstract type AbstractBasis end
Base.iterate(B::AbstractBasis, state...) = iterate(B.quanta, state...)
Base.length(B::AbstractBasis) = length(B.quanta)
Base.size(B::AbstractBasis) = size(B.quanta)

"""
VibBasis

Vibrational basis set, with a vector `quanta` of vibrational quantum numbers `vee`
representing vibrational eigenstates. 
"""
struct VibBasis <: AbstractBasis
eigen::Eigen{Float64, Float64, Matrix{Float64}, Vector{Float64}}
quanta::Vector{NamedTuple{(:vee,), Tuple{Int64}}}
end

"""
VibBasis(mu, dr, poten; nvee=nothing)

Create a sinc-DVR vibrational basis set, with a vector `quanta` of vibrational quantum numbers `vee`
representing vibrational eigenstates of the potential `poten` for a molecule with reduced mass `mu`
on a vibrational grid with spacing `dr`.

# Examples
```jldoctest
julia> VibBasis(1.0, 0.5, ([1., 1.5, 2., 2.5, 3] .- 2).^2)
````
"""
function VibBasis(mu, dr, poten; nvee=nothing)
    nvee = (isnothing(nvee) || nvee>length(poten)) ? length(poten) : nvee
    vibham = sinc_dvr(mu, dr, poten)
    vibeig = eigen(vibham, 1:nvee)
    quanta = Vector{NamedTuple{(:vee,), Tuple{Int64}}}(undef, nvee)
    for i = 1:length(vibeig.values)
        vee = i-1
        quanta[i] = (vee=vee,)
    end
    return VibBasis(vibeig, quanta)
end
VibBasis(mu, dr, poten::AbstractPotential; nvee=nothing) = VibBasis(mu, dr, poten.values; nvee=nvee)

"""
    EleBasis(lambda, ess, gerade, symmetric)

Create the electronic basis for a given diatomic state with orbital angular momentum projection `lambda`,
spin angular momentum `ess`, g/u symmetry, and optionally +/- symmetry for ``Σ`` electronic states.

# Examples
```jldoctest
julia> EleBasis(1.0, 1.0, true)
julia> EleBasis(0, 0, true, false)
````
"""
struct EleBasis <: AbstractBasis
    quanta::Vector{NamedTuple{(:lambda, :ess, :sigma), Tuple{Float64, Float64, Float64}}}
    multiplicity::Int64
    term::String
end

@inline lambda_to_char(lambda::Real) = "ΣΠΔΦΓΗΙ"[2*abs(round(Int, lambda))+1]

function EleBasis(ess, lambda, gerade::Bool)
    lambda == 0 && throw(ArgumentError("± symmetry (`symmetric = true` or `false`) is required for ``Σ``` electronic states.")) 
    multi = round(Int, 2*ess + 1)
    quanta = []
    for lam in (lambda==0. ? [lambda] : [-abs(lambda), abs(lambda)])
        for sigma = -ess:ess
            push!(quanta, (lambda=lam, ess=ess, sigma=sigma))
        end
    end
    #construct term symbol
    gu = gerade ? "g" : "u"
    term = string(multi, lambda_to_char(lambda), gu)
    return EleBasis(quanta, multi, term)
end

function EleBasis(ess, lambda, gerade::Bool, symmetric::Bool)
    if lambda != 0
        #if +/- given for non-sigma state, default to sigma constructor
        @warn "± symmetry not required for ``Λ > 0`` electronic states"
        return EleBasis(ess, lambda, gerade)
    else
        multi = round(Int, 2*ess + 1)
        quanta = []
        for lam in (lambda==0. ? [lambda] : [-abs(lambda), abs(lambda)])
            for sigma = -ess:ess
                push!(quanta, (lambda=lam, ess=ess, sigma=sigma))
            end
        end
        #construct term symbol
        gu = gerade ? "g" : "u"
        pm = symmetric ? "+" : "-"
        term = string(multi, lambda_to_char(lambda), gu, pm)
        return EleBasis(quanta, multi, term)
    end
end

"""
    RotBasis(jay, elebasis)

Create a rotational basis for a given total angular momentum quantum number `jay`
according to Hund's case a) from the electronic basis `elebasis`.
"""
struct RotBasis <: AbstractBasis
    quanta::Vector{NamedTuple{(:jay, :omega), Tuple{Float64, Float64}}}
end

function RotBasis(jay, elebasis)
    quanta = []
    for qnums in elebasis.quanta
        (lambda, ess, sigma) = qnums
        omega = lambda + sigma
        if abs(omega) > jay
            continue
        else
            push!(quanta, (jay=jay, omega=omega))
        end
    end
    return RotBasis(quanta)
end 

"""
    RovibronicBasis(elebasis, vibbasis, rotbasis)

Create a rovibronic basis from a given set of electronic, vibrational and rotational bases.
"""
struct RovibronicBasis <: AbstractBasis
    quanta::RovibronicQuanta
end

function RovibronicBasis(elebasis, vibbasis, rotbasis)
    quanta = []
    for ele in elebasis.quanta
        for vib in vibbasis.quanta
            for rot in rotbasis.quanta
                ele.lambda + ele.sigma != rot.omega && continue
                push!(quanta, (
                    jay=rot.jay, omega=rot.omega, 
                    lambda=ele.lambda, ess=ele.ess, sigma=ele.sigma, 
                    vee=vib.vee
                    ))
            end
        end
    end
    return RovibronicBasis(quanta)
end