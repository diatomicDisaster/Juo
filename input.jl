module ReadInput

using ArgParse

function parse_commandline()
  cla = ArgParseSettings()
  @add_arg_table! cla begin
    "input"
      help = "Input file, mandatory positional argument."
      required = true
    "--output", "-o"
      help = "Output file"
  end
  return parse_args(cla)
end

args = parse_commandline()
for (arg, val) in args
  println(arg, val)
end

end