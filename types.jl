struct Potential
    values   ::Vector{Float64}
    lambda   ::Float64
    ess      ::Float64
    symmetry ::Bool
    parity   ::Bool
end

Base.broadcastable(f::Potential) = Ref(f)

struct Coupling
    states::Tuple{Potential, Potential}
    values::Vector{Float64}
end

abstract type AbstractBasis end

struct Electronic <: AbstractBasis
    poten::Potential
    ess::Float
    sigma::Float64
    lambda::Float64
end

struct Vibrational <: AbstractBasis
    vee::Int
    vec::Vector{Float64}
end

struct Rotational <: AbstractBasis
    jay::Float64
    omega::Float64
    em::Float64
end

struct Rovibronic <: AbstractBasis
    e :: Electronic
    v :: Vibrational
    r :: Rotational
end

struct Hamiltonian
    vibsep::Float64
    vibgrid::Vector{Float64}
    reducedmass::Float64
    potentials::Tuple{Potential}
    couplings::Tuple{Coupling}
    matrix::Matrix{Float64}
end