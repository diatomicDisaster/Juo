module Grids

include("data.jl")

using Dierckx

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
    diffs = diff(points_)
    if all(p -> p≈first(diffs), diffs)
        return UniformGrid(points_[1], points_[npoints], diff(points_[1:2])[1], npoints, points_, values[sortinds])
    else
        return Grid(points_[1], points_[npoints], diff(points_), npoints, points_, values[sortinds])
    end
end


function Interpolate(grid::AbstractGrid, to_points::Array{Float64})
    spline = Spline1D(grid.points, grid.values)
    return Grid(to_points, spline(to_points))
end

function Grid(points::AbstractArray)
    sort!(points)
    
end

end