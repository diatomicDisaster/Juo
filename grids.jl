module Grids

include("data.jl")

using Interpolations

abstract type AbstractGrid end

struct Grid <: AbstractGrid
    rmin    :: Float64
    rmax    :: Float64
    rsep    :: Array{Float64}
    npoints :: Int64
    points  :: Array{Float64}
    values  :: Array{Float64}
end

struct UniformGrid <: AbstractGrid
    rmin    :: Float64
    rmax    :: Float64
    rsep    :: Float64
    npoints :: Int64
    points  :: Array{Float64}
    values  :: Array{Float64}
end

function Grid(points::Array{Float64}, values::Array{Float64})
    npoints = length(points)
    sortinds = sortperm(points)
    points_ = points[sortinds]
    if all(p -> p==points_[1], diff(points_))
        return UniformGrid(points_[1], points_[npoints], diff(points_[1:2])[1], npoints, points_, values[sortinds])
    else
        return Grid(points_[1], points_[npoints], diff(points_), npoints, points_, values[sortinds])
    end
end

function Interpolate(grid::UniformGrid, to_points::Array{Float64})
    itp = interpolate(grid.values, BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp, grid.points)
    return Grid(topoints, sitp.(to_points))
end

# function Interpolate(grid::Grid, to_points::Array{Float64})
#     itp = interpolate((grid.points,), grid.values, Gridded(Cubic(Line())))
#     return Grid(topoints, itp.(to_points))
# end

end