module Juo
include("Input.jl")

input_f  = stdin
output_f = stdout

if size(ARGS, 1) !== 0
  (input_f, output_f) = ReadInput.parse_args()
end

for (l, line) in enumerate(eachline(input_f))
  line = uppercase(line)
  println("Found $line")
end

end