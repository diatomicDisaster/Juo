abstract type AbstractGrid end

struct Potential <: AbstractGrid
    label      ::String
    values     ::Vector{Float64}
    geometries ::Vector{Float64}
    lambda     ::Float64
    ess        ::Float64
    symmetry   ::Bool
    parity     ::Bool
end

Base.broadcastable(f::Potential) = Ref(f)

struct Coupling <: AbstractGrid
    potentials ::Tuple{Potential, Potential}
    values     ::Vector{Float64}
    geometries ::Vector{Float64}
    lambdas    ::Tuple{Float64, Float64}
    ess        ::Tuple{Float64, Float64}
    symmetry   ::Tuple{Bool, Bool}
    parity     ::Tuple{Bool, Bool}
end

Base.size(G::AbstractGrid) = size(G.values)

abstract type AbstractState end

struct Electronic <: AbstractState
    potential::Potential
    lambda::Float64
    ess::Float64
    sigma::Float64
end

Electronic(potential::Potential) = Electronic.(potential, potential.lambda, potential.ess, collect(-potential.ess:potential.ess))  

struct Vibrational <: AbstractState
    potential::Potential
    vee::Int64
    energy::Float64
    vec::Vector{Float64}
end

struct Rotational <: AbstractState
    jay::Float64
    omega::Float64
    em::Float64
end

struct Rovibronic <: AbstractState
    e :: Electronic
    v :: Vibrational
    r :: Rotational
end

struct Molecule
    reducedmass::Float64
    potentials::Vector{Potential}
end

function sincdvr_hamiltonian(reducedmass::Real, gridgeometries::Vector, potentialvalues::Vector)
    dr = gridgeometries[2]-gridgeometries[1]
    all(v -> v≈gridgeometries[2]-gridgeometries[1], diff(gridgeometries)) || error("Grid not uniformly spaced")
    npoints = size(potentialvalues, 1)
    kin_oper = Array{Float64}(undef, (npoints, npoints))
    for j = 1:npoints
        for i = 1:j-1
            kin_oper[i, j] = (1 / (2*reducedmass*dr^2)) * 2*( ((-1)^(i - j)) / (i - j)^2 )
        end
        kin_oper[j, j] = (1 / (2*reducedmass*dr^2)) * (pi^2 / 3)
    end
    return Symmetric(kin_oper) + Diagonal(potentialvalues)
end

function sincdvr_basis(reducedmass::Real, poten::Potential; vmax=missing)
    vibhamiltonian = sincdvr_hamiltonian(reducedmass, poten.geometries, poten.values)
    eig = eigen(vibhamiltonian)
    vibbasis = Vector{Vibrational}(undef, length(eig.values))
    for v=1:length(eig.values)
        vee = v - 1
        vibbasis[v] = Vibrational(poten, vee, eig.values[v], eig.vectors[:, v])
    end
    return vibbasis
end

sincdvr_basis(reducedmass::Float64, potentials::Vector{Potential}) = sort(reduce(vcat, sincdvr_basis.(reducedmass, potentials)), by = vibstate -> vibstate.energy)

function electronic_basis(potential::Potential)
    ndimen = Int(2*(2*potential.ess + 1))
    elebasis = Vector{Electronic}(undef, ndimen)
    elebasis[1:Int(ndimen/2)] = Electronic.(potential, -potential.lambda, potential.ess, -potential.ess:potential.ess)
    elebasis[Int(ndimen/2)+1:ndimen] = Electronic.(potential, potential.lambda, potential.ess, -potential.ess:potential.ess)
    return elebasis
end

electronic_basis(potentials::Vector{Potential}) = reduce(vcat, electronic_basis.(potentials))

function vibrational_matrixelements(coupling::Coupling, vibbasis::Vector{Vibrational})
    matrixelements = zeros(Float64, size(vibbasis, 1), size(vibbasis, 1))
    gridspacing = diff(coupling.geometries)
    for i in CartesianIndices(matrixelements)
        finalvib = vibbasis[i[1]]
        initialvib = vibbasis[i[2]]
        if (finalvib.potential == coupling.potentials[1]) && (initialvib.potential == coupling.potentials[2]) 
            couplingintegral = finalvib.vector .* coupling.values .* initialvib.vector
            matrixelements[i] = sum(gridspacing .* (couplingintegral[1:size(coupling, 1)-1] .+ couplingintegral[2:size(coupling, 1)]) / 2)
        else
            matrixelements[i] = 0
        end
    end
end