module ReadInput

export read_input

struct InputError <: Exception
  var :: String
end

Base.showerror(io::IO, e::InputError) = print(io, "Unrecognised input ", e.var)

"""Parse the command line arguments for the input and output files, if there
are any, and set global I/O variables accordingly.
"""
function parse_args()
  
  fin=stdin
  fout=stdout
  next = iterate(ARGS)
  while next !== nothing
    (arg, state) = next
    if arg in ("--input", "-i")
      (arg, state) = iterate(ARGS, state)
      fin = open(arg, "r")
      next = iterate(ARGS, state)
    elseif arg in ("--output", "-o")
      (arg, state) = iterate(ARGS, state)
      fout = open(arg, "w")
      next = iterate(ARGS, state)
    else
      throw(InputError(arg))
      next = iterate(ARGS, state)
    end
  end
  return (fin, fout)
end

function read_input()
  (fin, fout) = parse_args()
  for (l, line) in enumerate(eachline(fin))
    line = uppercase(line)
    write(fout, line)
  end
  return (fin, fout)
end

end