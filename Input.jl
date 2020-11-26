module ReadInput

export parse_args

struct InputError <: Exception
  var :: String
end

Base.showerror(io::IO, e::InputError) = print(io, "Unrecognised input ", e.var)

function parse_args()
  """Parse the command line arguments for the input and output files, if there
  are any, and set global I/O variables accordingly.
  """
  next = iterate(ARGS)
  while next !== nothing
    (arg, state) = next
    if arg in ("--input", "-i")
      (arg, state) = iterate(ARGS, state)
      input_f = open(arg, "r")
      next = iterate(ARGS, state)
    elseif arg in ("--output", "-o")
      (arg, state) = iterate(ARGS, state)
      output_f = open(arg, "w")
      next = iterate(ARGS, state)
    else
      throw(InputError(arg))
      next = iterate(ARGS, state)
    end
  end
  return (input_f, output_f)
end

end