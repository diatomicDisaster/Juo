module Basis
include("grids.jl")
using LinearAlgebra
using .Grids

export VibBasis, EleBasis, RotBasis, RovibronicBasis

"""
    sinc_dvr(mu, dr, poten)

Build the vibrational Hamiltonian using the sinc-DVR method for a molecule with reduced mass `mu`, 
and a vibrational grid spaceing `dr` with a series of potential points `poten`.

# Examples
```jldoctest
julia> sinc_dvr(1.0, 0.5, ([1., 1.5, 2., 2.5, 3.].- 2.).^2)
5×5 Symmetric{Float64, Matrix{Float64}}:
  7.57974   -4.0        1.0      -0.444444   0.25
 -4.0        6.82974   -4.0       1.0       -0.444444
  1.0       -4.0        6.57974  -4.0        1.0
 -0.444444   1.0       -4.0       6.82974   -4.0
  0.25      -0.444444   1.0      -4.0        7.57974
```
"""
function sinc_dvr(mu, dr, poten)
    npoints = length(poten)
    kinetic_operator = Matrix{Float64}(undef, npoints, npoints)
    for j = 1:npoints
        for i = 1:j-1
            kinetic_operator[i, j] = (1/(mu*dr^2)) * (-1)^(i - j)/(i - j)^2
        end
        kinetic_operator[j, j] = (1/(2*mu*dr^2)) * (pi^2 / 3)
    end
    return Symmetric(kinetic_operator) + Diagonal(poten)
end

"""
VibBasis

Vibrational basis set, with a vector `quanta` of vibrational quantum numbers `vee`
representing vibrational eigenstates. 
"""
struct VibBasis
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
VibBasis(mu, dr, poten::AbstractPotential; kwargs...) = VibBasis(mu, dr, poten.values; kwargs...)

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
struct EleBasis
    quanta::Vector{NamedTuple{(:lambda, :ess, :sigma), Tuple{Float64, Float64, Float64}}}
    multiplicity::Int64
    term::String
end

@inline lambda_to_char(lambda::Real) = "ΣΠΔΦΓΗΙ"[2*abs(round(Int, lambda))+1]

function EleBasis(ess, lambda, gerade::Bool)
    lambda == 0 && throw(ArgumentError("± symmetry (`symmetric = true` or `false`) is required for ``Σ``` electronic states.")) 
    multi = round(Int, 2*ess + 1)
    quanta = Vector{NamedTuple}(undef, multi*round(Int, 2*lambda))
    i = 1
    for lam in (lambda==0. ? [lambda] : [-abs(lambda), abs(lambda)])
        for sigma = -ess:ess
            quanta[i] = (lambda=lam, ess=ess, sigma=sigma)
            i += 1
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
        quanta = Vector{NamedTuple}(undef, multi*round(Int, 2*lambda))
        i = 1
        for lam in (lambda==0. ? [lambda] : [-abs(lambda), abs(lambda)])
            for sigma = -ess:ess
                quanta[i] = (lambda=lam, ess=ess, sigma=sigma)
                i += 1
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
struct RotBasis
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

struct RovibronicBasis
    quanta::Vector{NamedTuple{(:jay, :omega, :lambda, :ess, :sigma, :vee), Tuple{Float64, Float64, Float64, Float64, Float64, Int64}}}
end

function RovibronicBasis(elebasis, vibbasis, rotbasis)
    quanta = []
    for ele in elebasis.quanta
        for vib in vibbasis.quanta
            for rot in rotbasis.quanta
                ele.omega != rot.omega && continue
                push!(quanta, (jay=rot.jay, omega=rot.omega, ele.lambda, ele.ess, ele.sigma, vib.vee))
            end
        end
    end
    return RovibronicBasis(quanta)
end

end