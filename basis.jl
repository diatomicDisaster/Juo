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
    quanta::Vector{EleQuanta}
    ess::Int64
    abslambda::Float64
    term::String
end

@inline lambda_to_char(lambda::Real) = "ΣΠΔΦΓΗΙ"[2*abs(round(Int, lambda))+1]

function EleBasis(ess, lambda, gerade::Bool)
    lambda == 0 && throw(ArgumentError("± symmetry (`symmetric = true` or `false`) is required for ``Σ``` electronic states.")) 
    #construct term symbol
    gu = gerade ? "g" : "u"
    multi = round(Int, 2*ess + 1)
    term = string(multi, lambda_to_char(lambda), gu)quanta = []
    for lam in (lambda==0. ? [lambda] : [-abs(lambda), abs(lambda)])
        for sigma = -ess:ess
            push!(quanta, (lambda=lam, ess=ess, sigma=sigma))
        end
    end
    return EleBasis(quanta, ess, abs(lambda), term)
end

function EleBasis(ess, lambda, gerade::Bool, symmetric::Bool)
    if lambda != 0
        #if +/- given for non-sigma state, default to sigma constructor
        @warn "± symmetry not required for ``Λ > 0`` electronic states"
        return EleBasis(ess, lambda, gerade)
    else
        #construct term symbol
        gu = gerade ? "g" : "u"
        pm = symmetric ? "+" : "-"
        multi = round(Int, 2*ess + 1)
        term = string(multi, lambda_to_char(lambda), gu, pm)
        quanta = []
        for lam in (lambda==0. ? [lambda] : [-abs(lambda), abs(lambda)])
            for sigma = -ess:ess
                push!(quanta, (lambda=lam, ess=ess, sigma=sigma))
            end
        end
        return EleBasis(quanta, ess, abs(lambda), term)
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
    range = elebasis.abslambda + (elebasis.multiplicity - 1)/2
    for omega = -range:range
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
    elebasis::EleBasis
    vibbasis::VibBasis
    rotbasis::RotBasis
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

function vibmatrix(vibvecs::Matrix, coupling::Vector, dr::Float64)
    mat = Matrix{Float64}(undef, size(vibvecs, 2), size(vibvecs, 2))
    for (i, veci) in enumerate(eachcol(vibvecs))
        pre_vec = coupling .* veci
        for (f, vecf) in enumerate(eachcol(vibvecs))
            mat[f, i] = dr * vecf⋅pre_vec
        end
    end
    return Hermitian(mat)
end

vib_matrix(vibbasis::VibBasis, coupling::AbstractCoupling, dr::Float64) = 
    vib_matrix(vibbasis.eigen.vectors, coupling.values, dr)

@inline couplhasstate(state::EleQuanta, coupl::CouplQuanta, fi::Char) =
    all(map(k->getindex(state, k)==getindex(coupl, Symbol(string(k, fi))), keys(state)))

function hamiltonian(elebasis::EleBasis, coupling::AbstractCoupling)
    ele_hamiltonian = zeros(length(elebasis), length(elebasis))
    f = findall(e->couplhasstate(e, coupling.quanta, 'f'), elebasis.quanta)
    i = findall(e->couplhasstate(e, coupling.quanta, 'i'), elebasis.quanta)
    if !(length(f) == 1 && length(i) == 1)
        throw(DomainError("wrong number of (f, i) states matching coupling quanta, (expected (1, 1), got ($(length(f)), $(length(i))))"))
    end
    electrontic_matrix[f[1], i[1]] = 1.0
    return ele_hamiltonian
end

function hamiltonian(elebasis::EleBasis, vibbasis::VibBasis, coupls::AbstractCoupling, dr)
    ndimen = size(elebasis) * size(vibbasis)
    vibronic_hamiltonian = zeros(ndimen, ndimen)
    for coupl in coupls
        vibronic_hamiltonian += 
            kron(ele_hamiltonian(elebasis, coupl), vib_hamiltonian(vibbasis, coupl, dr))
    end
    vibronic_hamiltonian += kron(Diagonal(ones(size(elebasis))), Diagonal(vibbasis.eigen.values))
    return vibronic_hamiltonian
end

function hamiltonian(rovibronicbasis::RovibronicBasis)
    hamiltonian = zeros(length(rovibronicbasis), length(rovibronicbasis))
    for (i, qi) in enumerate(rotbasis.quanta)
        for (f, qf) in enumerate(rotbasis.quanta)
            if qf.omega - qi.omega == 0
                rot_hamiltonian[i, f] = qi.omega + qi.jay*(qi.jay + 1)
            elseif abs(qf.omega - qi.omega) == 1
                rot_hamiltonian[i, f] = (qi.jay*(qi.jay + 1) - qi.omega*qf.omega)
            end
        end
    end
end
