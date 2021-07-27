"""
    sinc_dvr(mu, dr, poten)

Build the vibrational Hamiltonian using the sinc-DVR method for a molecule with reduced mass `mu`, 
and a vibrational grid spacing `dr` with a series of potential points `poten`.

For details of the sinc-DVR method, see: 
D. Colbert and W. Miller, J. Chem. Phys., 96 (1992), doi:10.1063/1.462100

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
Builds the raising operator matrix for a given angular momentum quantum number.
"""
function raising(J)
    ndimensions = Int(2*J + 1)
    rows = Int[i for i=2:ndimensions]
    cols = Int[i for i=1:ndimensions-1]
    vals = Array{Float64, 1}(undef, ndimensions-1)
    for (i, M) in enumerate(-J:J-1)
        vals[i] = (J*(J + 1) - M*(M + 1))^.5
    end
    sparse_raising = sparse(rows, cols, vals, ndimensions, ndimensions)
    return sparse_raising
end

"""
Builds the lowering operator matrix for a given angular momentum quantum number.
"""
function lowering(J)
    ndimensions = Int(2*J + 1)
    rows = Int[i for i=2:ndimensions]
    cols = Int[i for i=1:ndimensions-1]
    vals = Array{Float64, 1}(undef, ndimensions-1)
    for (i, M) in enumerate(-J+1:J)
        vals[i] = (J*(J + 1) - M*(M - 1))^.5
    end
    sparse_lowering = sparse(rows, cols, vals, ndimensions, ndimensions)
    return sparse_lowering
end

"""
Builds the total momentum operator matrix for a given angular momentum quantum number.
"""
function momentum_squared(J)
    ndimensions = Int(2*J + 1)
    rows = Int[i for i=1:ndimensions]
    cols = Int[i for i=1:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions)
    for (i, M) in enumerate(-J:J)
        vals[i] = J*(J + 1)
    end
    sparse_momentum_squared = sparse(rows, cols, vals, ndimensions, ndimensions)
    return sparse_momentum_squared
end

"""
Builds the momentum projection operator matrix for a given angular momentum quantum number.
"""
function momentum_projection(J)
    ndimensions = Int(2*J + 1)
    rows = Int[i for i=1:ndimensions]
    cols = Int[i for i=1:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions)
    for (i, M) in enumerate(-J:J)
        vals[i] = M
    end
    sparse_momentum_projection = sparse(rows, cols, vals, ndimensions, ndimensions)
    return sparse_momentum_projection
end