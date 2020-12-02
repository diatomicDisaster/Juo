module Juo
include("Input.jl")
include("Functions.jl")

#Read input file
(fIn, fOut) = ReadInput.read_input()

#Define grid
const N = 501
vib_grid = LinRange(0.1, 5.0, N)
vib_grid_a = collect(LinRange(0.1, 5.0, N))

println("Internal vectorization as array")
@time pot_grid = PotentialCurves.mor_pot(vib_grid_a; a=1.0, r_e=1.0, D_e=1.0)
@time pot_grid = PotentialCurves.mor_pot(vib_grid_a; a=1.0, r_e=1.0, D_e=1.0)
println("Internal vectorization as LinRange")
@time pot_grid = PotentialCurves.mor_pot(vib_grid; a=1.0, r_e=1.0, D_e=1.0)
@time pot_grid = PotentialCurves.mor_pot(vib_grid; a=1.0, r_e=1.0, D_e=1.0)
println("Internal vectorization as loop")
@time pot_grid = PotentialCurves.mor_pot_loop(vib_grid_a; a=1.0, r_e=1.0, D_e=1.0)
@time pot_grid = PotentialCurves.mor_pot_loop(vib_grid_a; a=1.0, r_e=1.0, D_e=1.0)
println("External vectorization")
@time pot_grid = @. PotentialCurves.mor_pot(vib_grid_a; a=1.0, r_e=1.0, D_e=1.0)
@time pot_grid = @. PotentialCurves.mor_pot(vib_grid_a; a=1.0, r_e=1.0, D_e=1.0)

# global bounds = [0, 0]
# global counter = 1
# global mode = "under"
# for (i, V) in enumerate(potential_grid)
#     if V < potential_grid[length(potential_grid)]
#         if mode == "under"
#             global bounds[1] = i
#             global mode = "over"
#         else
#             continue
#         end
#     else
#         if mode == "over"
#             global bounds[2] = i
#             break
#         else
#             continue
#         end
#     end
# end

end