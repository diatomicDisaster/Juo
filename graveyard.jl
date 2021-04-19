
function vibronicbasis(
    veemax_s :: Vector{Int64},
    lambda_s :: Vector{Float64},
    ess_s    :: Vector{Float64},
    jay_list :: Vector{Float64}
    )
    length(veemax_s) == length(lambda_s) == length(ess_s) || throw(ArgumentError("Length of veemax_s (=$(length(veemax_s))), lambda_s (=$(length(lambda_s))), and ess_s (=$(length(ess_s))) should be equal."))
    njays    = length(jay_list)
    nstates  = length(veemax_s)
    essdimen_s = Array{Int64}(undef, nstates)
    vibronicdimen_s = Array{Int64}(undef, nstates)
    @. essdimen_s = floor(Int64, 2*ess_s + 1)
    @. vibronicdimen_s = (veemax_s + 1) * essdimen_s
    vecs = Vector{BasisVector}(undef, sum(vibronicdimen_s)*njays)
    for (j, jay) in enumerate(jay_list)
        for state = 1:nstates
            lambda = lambda_s[state]
            ess    = ess_s[state]
            veemax = veemax_s[state]
            for (s, sigma) in enumerate(-ess:ess)
                omega = lambda + sigma
                for vee = 0:veemax
                    i = (j - 1) * (sum(vibronicdimen_s)) + sum(vibronicdimen_s[1:(state-1)]) + (s-1)*(veemax+1) + (vee+1)
                    vecs[i] = BasisVector(state, vee, lambda, ess, sigma, jay, omega)
                end
            end
        end
    end
    return vecs
end


"""
Builds the raising operator matrix for the basis of angular momentum quantum numbers provided.
    Arguments
        Jlist :: Array{Float64, 1}
            List of angular moment quantum numbers to build the matrix for, can be non-sequential.
        ndimensions :: Int64
            The dimension of the rotational basis. For a complete set of rotational quantum 
            numbers, this is equal to (J_max + 1)^2 if J is integer, or Jmax(Jmax + 3) + 2 
            if J is half-integer.
    Returns
        SparseArrays.AbstractSparseMatrixCSC{Float64, 2}
            CSC sparse matrix representing the raising operator for the rotational basis.
"""
function raising_operator(
    Jlist :: Array{Float64, 1},
    ndimensions :: Int64
    )
    rows = Int64[i for i=2:ndimensions]
    cols = Int64[i for i=1:ndimensions-1]
    vals = Array{Float64, 1}(undef, ndimensions-1)
    i = 1
    for J in Jlist
        for M = -J:J-1
            vals[i] = (J*(J + 1) - M*(M + 1))^.5
            i += 1
        end
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function raising_operator(
    J :: Float64,
    ndimensions :: Int64
    )
    rows = Int64[i for i=2:ndimensions]
    cols = Int64[i for i=1:ndimensions-1]
    vals = Array{Float64, 1}(undef, ndimensions-1)
    i = 1
    for M = -J:J-1
        vals[i] = (J*(J + 1) - M*(M + 1))^.5
        i += 1
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end


"""
Builds the lowering operator matrix for the basis of angular momentum quantum numbers provided.
    Arguments
        Jlist :: Array{Float64, 1}
            List of angular moment quantum numbers to build the matrix for, can be non-sequential.
        ndimensions :: Int64
            The dimension of the rotational basis. For a complete set of rotational quantum 
            numbers, this is equal to (J_max + 1)^2 if J is integer, or Jmax(Jmax + 3) + 2 
            if J is half-integer.
    Returns
        SparseArrays.AbstractSparseMatrixCSC{Float64, 2}
            CSC sparse matrix representing the lowering operator for the rotational basis.
"""
function lowering_operator(
    Jlist :: Array{Float64, 1},
    ndimensions :: Int64
    )
    rows = Int64[i for i=1:ndimensions-1]
    cols = Int64[i for i=2:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions-1)
    i = 0
    for J in Jlist
        for M = -J:J
            if i > 0
                vals[i] = (J*(J + 1) - M*(M - 1))^.5
                i += 1
            end
        end
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function lowering_operator(
    J :: Float64,
    ndimensions :: Int64
    )
    rows = Int64[i for i=2:ndimensions]
    cols = Int64[i for i=1:ndimensions-1]
    vals = Array{Float64, 1}(undef, ndimensions-1)
    i = 0
    for M = -J:J
        if i > 0
            vals[i] = (J*(J + 1) - M*(M - 1))^.5
            i += 1
        end
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end


"""
Builds the squared momentum operator matrix for the basis of angular momentum quantum numbers provided.
    Arguments
        Jlist :: Array{Float64, 1}
            List of angular moment quantum numbers to build the matrix for, can be non-sequential.
        ndimensions :: Int64
            The dimension of the rotational basis. For a complete set of rotational quantum 
            numbers, this is equal to (J_max + 1)^2 if J is integer, or Jmax(Jmax + 3) + 2 
            if J is half-integer.
    Returns
        SparseArrays.AbstractSparseMatrixCSC{Float64, 2}
            CSC sparse matrix representing the squared momentum operator for the rotational basis.
"""
function momentum_operator(
    Jlist :: Vector{Float64},
    ndimensions :: Int64
    )
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions)
    i = 1
    for J in Jlist
        for M = -J:J
            vals[i] = J*(J + 1)
            i += 1
        end
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function momentum_operator(
    J :: Float64,
    ndimensions :: Int64
    )
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions)
    i = 1
    for M = -J:J
        vals[i] = J*(J + 1)
        i += 1
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end


"""
Builds the momentum projection operator matrix for the basis of angular momentum quantum numbers provided.
    Arguments
        Jlist :: Array{Float64, 1}
            List of angular moment quantum numbers to build the matrix for, can be non-sequential.
        ndimensions :: Int64
            The dimension of the rotational basis. For a complete set of rotational quantum 
            numbers, this is equal to (J_max + 1)^2 if J is integer, or Jmax(Jmax + 3) + 2 
            if J is half-integer.
    Returns
        AbstractSparseMatrixCSC{Float64, 2}
            CSC sparse matrix representing the momentum projection operator for the rotational basis.
"""
function momentum_projection_operator(
    Jlist :: Array{Float64, 1},
    ndimensions :: Int64
    )
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions)
    i = 1
    for J in Jlist
        for M = -J:J
            vals[i] = M
            i += 1
        end
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end

function momentum_projection_operator(
    J :: Float64,
    ndimensions :: Int64
    )
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions)
    i = 1
    for M = -J:J
        vals[i] = M
        i += 1
    end
    return sparse(rows, cols, vals, ndimensions, ndimensions)
end


"""
Builds the rigid-rotor Hamiltonian for the specified list of angular momentum quantum numbers,
and given value of total spin angular momentum.
    Arguments
        Jlist :: Array{Float64, 1}
            List of angular moment quantum numbers to build the matrix for, can be non-sequential.
        spin :: Float64
            Total spin angular momentum for the system.
        mass :: Float64
            Reduced mass of the molecule.
        r :: Float64
            Rigid rotor length (internuclear distance).
    Returns
        SparseArrays.AbstractSparseMatrixCSC{Float64, 2}
            CSC sparse matrix representing the momentum projection operator for the rotational basis.
"""
function hamiltonian(
    jaylist :: Array{Float64, 1},
    spin  :: Float64,
    reduced_mass  :: Float64,
    r     :: Float64
    )

    ndimensions = _countdimensions(jaylist)
    ndimensions_spin = floor(Int64, 2*spin + 1)
    
    #Build the J, M operators in the tensor product basis
    identity = Matrix{Float64}(1.0I, (ndimensions_spin, ndimensions_spin))
    jaysq = momentum_operator(jaylist, ndimensions)
    em = momentum_projection_operator(jaylist, ndimensions)
    jaysq_spinid = kron(jaysq, identity)
    em_spinid = kron(em, identity)

    #Build the S, Σ operators in the tensor product basis
    identity = Matrix{Float64}(I, (ndimensions, ndimensions))
    spinsq = momentum_operator(spin, ndimensions_spin)
    sigma = momentum_projection_operator(spin, ndimensions_spin)
    jayid_spinsq = kron(identity, spinsq)
    jayid_sigma = kron(identity, sigma)

    #Build S-uncoupling operators in the tensor product basis
    jaypl = raising_operator(jaylist, ndimensions)
    jaymi = lowering_operator(jaylist, ndimensions)
    spinpl = raising_operator(spin, ndimensions_spin)
    spinmi = lowering_operator(spin, ndimensions_spin)
    jaypl_spinmi = kron(jaypl, spinmi)
    jaymi_spinpl = kron(jaypl, spinmi)

    #Return sum of Hamiltonian terms
    return ((jaysq_spinid - em_spinid.^2) + (jayid_spinsq - jayid_sigma.^2) + (jaypl_spinmi + jaymi_spinpl))/(2*reduced_mass*r^2)
end


function test_vibrational()
    #Morse potential parameters for CO
    # A. Roy (2013) - https://doi.org/10.1016/j.rinp.2013.06.001
    D_e = 0.412533191
    mu  = 12506.23981
    r_e = 2.132177990 
    a   = 2.59441

    #Define grid (same as original Duo paper)
    # S. Yurchenko et al. (2016) - https://doi.org/10.1016/j.cpc.2015.12.021
    rmin = 1.322808287
    rmax = 3.77845225
    dr = 0.01
    thresh = 50.
    show_levels = 30

    vib_grid = collect(rmin:dr:rmax)
    npoints = length(vib_grid)
    println("No. initial grid points: $npoints")
    pot_grid = Potentials.morse.(vib_grid; r_e=r_e, a=a, D_e=D_e)

    vib_hamiltonian = Vibrational.hamiltonian(pot_grid, dr, mu)
    vibrational_eigen = eigen(vib_hamiltonian)
    ngoodpoints = length(vibrational_eigen.values[1:show_levels])
    println("No. retained grid points: $ngoodpoints")

    plot(vibrational_grid, potential_grid)
    savefig(joinpath("images", "potential.png"))

    heatmap(1:ngoodpoints, 1:ngoodpoints, vibrational_hamiltonian)
    savefig(joinpath("images", "hamiltonian.png"))

    heatmap(1:ngoodpoints, 1:show_levels, transpose(vibrational_eigen.vectors[:,1:show_levels].^2))
    savefig(joinpath("images", "eigenvectors.png"))
end

function test_rotational()
    Js = [i for i=0.:20.]
    rot_hamiltonian = Rotational.hamiltonian(Js, 0.5, 2.0, 1.0)
    rot_eigen = Rotational.solve(rot_hamiltonian)
    xs = [0, 1]
    p = plot(xs, [0., 0.])
    for val in real(rot_eigen.values)
        plot!(p, xs, [val, val])
    end
    savefig(joinpath("images", "rotational.png"))
end

function test_rovibrational()
    #Morse potential parameters for CO
    # A. Roy (2013) - https://doi.org/10.1016/j.rinp.2013.06.001
    D_e = 0.412533191
    mu  = 12506.23981
    r_e = 2.132177990 
    a   = 2.59441

    #Define grid (same as original Duo paper)
    # S. Yurchenko et al. (2016) - https://doi.org/10.1016/j.cpc.2015.12.021
    rmin = 1.322808287
    rmax = 3.77845225
    dr = 0.01
    thresh = 50.
    vmax = 2

    ndimen_vib = vmax + 1
    vib_grid = collect(rmin:dr:rmax)
    npoints = length(vib_grid)
    println("No. vibrational points: $npoints")
    pot_grid = Potentials.morse.(vib_grid; r_e=r_e, a=a, D_e=D_e)

    vib_ham = Vibrational.sincdvr(pot_grid, dr, mu)
    vib_eigen = eigen(vib_ham)
    vib_hamiltonian = Diagonal(vib_eigen.values[1:ndimen_vib])

    Js = [i for i=0.:5.]
    rot_ham = Rotational.hamiltonian(Js, 0., mu, r_e)
    rot_hamiltonian = Matrix(rot_ham)
    ndimen_rot = size(rot_hamiltonian)[1]

    identity = Matrix{Float64}(I, (ndimen_rot, ndimen_rot))
    vibrational_hamiltonian = kron(vib_hamiltonian, identity)
    identity = Matrix{Float64}(I, (ndimen_vib, ndimen_vib))
    rotational_hamiltonian  = kron(identity, rot_hamiltonian)

    rovib_hamiltonian = vibrational_hamiltonian + rotational_hamiltonian
    rovib_eigen = eigen(rovib_hamiltonian)
    
    reals = real(rovib_eigen.values)
    xs = [0, 1]
    p = plot(xs, [reals[1], reals[1]])
    for (ind, val) in enumerate(reals)
        if ind != 1
            plot!(p, xs, [val, val])
        end
    end
    savefig(joinpath("images", "rovibrational.pdf"))
    plot!(p, ylims=[0.01, 0.011])
    savefig(joinpath("images", "rovibrational_v0.pdf"))

end
