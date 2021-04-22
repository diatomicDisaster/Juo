"""Build and diagonalise the vibrational Hamiltonian for a uniformly spaced 
potential energy grid using the sinc-DVR kinetic energy operator.

For details of the sinc-DVR method, see: 
D. Colbert and W. Miller, J. Chem. Phys., 96 (1992), doi:10.1063/1.462100
"""
function sincdvr(potential_grid::Array{Float64, 1}, dr::Float64, mu::Float64)
    npoints = length(potential_grid)
    kinetic_operator = Array{Float64}(undef, (npoints, npoints))
    for j = 1:npoints
        for i = 1:j-1
            kinetic_operator[i, j] = (1 / (2*mu*dr^2)) * 2*( ((-1)^(i - j)) / (i - j)^2 )
        end
        kinetic_operator[j, j] = (1 / (2*mu*dr^2)) * (pi^2 / 3)
    end
    return Symmetric(kinetic_operator) + Diagonal(potential_grid)
end

"""
Builds the raising operator matrix for a given angular momentum quantum number.
"""
function raising(J :: Float64)
    ndimensions = Int64(2*J + 1)
    rows = Int64[i for i=2:ndimensions]
    cols = Int64[i for i=1:ndimensions-1]
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
function lowering(J :: Float64, ndimensions :: Int64)
    rows = Int64[i for i=2:ndimensions]
    cols = Int64[i for i=1:ndimensions-1]
    vals = Array{Float64, 1}(undef, ndimensions-1)
    for (i, M) in enumerate(-J+1:J)
        vals[i] = (J*(J + 1) - M*(M - 1))^.5
    end
    sparse_lowering = sparse(rows, cols, vals, ndimensions, ndimensions)
    return sparse_lowering
end

lowering(J::Float64) = lowering(J, Int64(2*J + 1))
lowering(J::Int64) = lowering(float(J), 2*J + 1)

"""
Builds the total momentum operator matrix for a given angular momentum quantum number.
"""
function momentum_squared(J :: Float64)
    ndimensions = Int64(2*J + 1)
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
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
function momentum_projection(J :: Float64)
    ndimensions = Int64(2*J + 1)
    rows = Int64[i for i=1:ndimensions]
    cols = Int64[i for i=1:ndimensions]
    vals = Array{Float64, 1}(undef, ndimensions)
    for (i, M) in enumerate(-J:J)
        vals[i] = M
    end
    sparse_momentum_projection = sparse(rows, cols, vals, ndimensions, ndimensions)
    return sparse_momentum_projection
end
