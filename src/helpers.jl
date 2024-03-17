using JSON

function output_to_string(
  meas::Dict{String,Tuple{Any,Any}},
  obs::Vector{Tuple{String,String}},
)
  """
  Convert the raw measurments (dictionary) into strings
  """
  data = String[]
  for pair in obs
    push!(data, parseObs(meas, pair...))
  end
  return join(data, ",")
end

function parseObs(meas::Dict{String,Tuple{Any,Any}}, name::String, type="mean")
  value = meas[name]
  value = type != "mean" ? getfield(Main, Symbol(type)).(value) : value
  return join(json.(json.(value)), ",")
end

function parseObs(meas::Dict{String,Any}, name::String, type="mean")
  value = meas[name]
  value = type != "mean" ? getfield(Main, Symbol(type)).(value) : value
  if eltype(meas[name]) <: Number
    return string(value)
  else
    return join(string.(value), " ")
  end
end

function file_to_dict(file)
  """
  Parse the option file into a dictionary
  """
  #float parameters
  floats = [
    "lambda1",
    "lambda2",
    "lambda3",
    "eta",
    "delta_theta",
    "chem_potential",
    "chem_disorder",
    "bond_disorder_villain",
    "bond_disorder_loop",
    "prob",
    "Lt_ratio",
    "Lt_power",
  ]
  #integer parameters
  ints = [
    "seed",
    "bond_type",
    "L",
    "num_of_thermal",
    "num_of_measure",
    "num_of_sweeps",
    "num_of_relative_sweeps",
    "dual_lambda",
    "nbin",
    "dual_L",
  ]

  params = Dict{String,Real}()

  open(file) do f
    lines = readlines(f)
    for line in lines
      k, v = split(line, ",")
      try
        if k in floats
          params[k] = parse(Float64, v)
        else
          params[k] = parse(Int64, v)
        end
      catch e
        println("Error parsing $k:$v")
        rethrow(e)
      end
    end
  end

  if params["dual_lambda"] == 1
    params["lambda2"] = params["lambda1"]
  end


  if params["dual_L"] == 1
    params["Lt"] = params["L"]
  else
    params["Lt"] = Int64(params["Lt_ratio"] * (params["L"]^params["Lt_power"]))
  end

  to_delete = ["Lt_ratio", "Lt_power", "dual_lambda", "dual_L"]

  for key in to_delete
    pop!(params, key)
  end

  return params
end


function print_header()
  """
  Print out the headers for the raw data in the CSV format
  """
  #All Observables with the properties to extract
  obs = [
    ("Energy", "mean"),
    ("Complex Hall Conductivity", "real"),
    ("Complex Hall Conductivity", "imag"),
    ("Stiffness Loop All", "mean"),
    ("Current Time Loop", "mean"),
    ("Current Space Loop", "mean"),
    ("Current Loop FT", "real"),
    ("Current Loop FT", "imag"),
    ("Compressibility Loop", "mean"),
    ("Green Loop Time", "mean"),
    ("Green Loop Space", "mean"),
    ("ZRatio", "mean"),
    ("Green Villain Time", "real"),
    ("Green Villain Time", "imag"),
    ("Green Villain Space", "real"),
    ("Green Villain Space", "imag"),
    ("Compressibility Villain", "mean"),
    ("Current Time Villain", "mean"),
    ("Current Space Villain", "mean"),
    ("Stiffness Villain", "mean"),
    ("Stiffness Villain FT", "mean"),
    ("Current Villain FT", "real"),
    ("Current Villain FT", "imag"),
    ("Magnitization", "mean"),
  ]
  header = ["sample_num"]
  for ob in obs
    push!(header, string(ob[1], "_", ob[2]))
    push!(header, string(ob[1], "_", ob[2], "_sq"))
  end
  println(join(header, ","))
end


function print_raw_output(sim::Sim)
  """
  Print out the raw data to STDOUT in the CSV format
  The data aligns with the header from print_header function
  """
  obs = [
    ("Energy", "mean"),
    ("Complex Hall Conductivity", "real"),
    ("Complex Hall Conductivity", "imag"),
    ("Stiffness Loop All", "mean"),
    ("Current Time Loop", "mean"),
    ("Current Space Loop", "mean"),
    ("Current Loop FT", "real"),
    ("Current Loop FT", "imag"),
    ("Compressibility Loop", "mean"),
    ("Green Loop Time", "mean"),
    ("Green Loop Space", "mean"),
    ("ZRatio", "mean"),
    ("Green Villain Time", "real"),
    ("Green Villain Time", "imag"),
    ("Green Villain Space", "real"),
    ("Green Villain Space", "imag"),
    ("Compressibility Villain", "mean"),
    ("Current Time Villain", "mean"),
    ("Current Space Villain", "mean"),
    ("Stiffness Villain", "mean"),
    ("Stiffness Villain FT", "mean"),
    ("Current Villain FT", "real"),
    ("Current Villain FT", "imag"),
    ("Magnitization", "mean"),
  ]
  data = MCMC.get_measurements(sim.measurements)
  sample_num = sim.measurements.measurements[4].obs_data.c_num_of_measure
  print(json(json(sample_num)), ",")
  println(output_to_string(data, obs))
end
