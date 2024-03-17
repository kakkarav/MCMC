using JSON

# Constants for parameter types
const FLOAT_PARAMS = [
  "lambda1", "lambda2", "lambda3", "eta", "delta_theta",
  "chem_potential", "chem_disorder", "bond_disorder_villain",
  "bond_disorder_loop", "prob", "Lt_ratio", "Lt_power"
]
const INT_PARAMS = [
  "seed", "bond_type", "L", "num_of_thermal", "num_of_measure",
  "num_of_sweeps", "num_of_relative_sweeps", "dual_lambda",
  "nbin", "dual_L"
]

function output_to_string(meas::Dict{String,Tuple{Any,Any}}, obs::Vector{Tuple{String,String}})
  """
  Convert the raw measurements (dictionary) into strings.
  """
  data = [parse_obs(meas, pair...) for pair in obs]
  return join(data, ",")
end

function parse_obs(meas::Dict{String,Tuple{Any,Any}}, name::String, stat_type="mean")
  value, _ = get(meas, name, (nothing, nothing))
  if stat_type != "mean"
    value = apply_stat_function(value, stat_type)
  end
  return value_to_string(value)
end

function apply_stat_function(value, stat_type)
  try
    func = getfield(Main, Symbol(stat_type))
    return func(value)
  catch
    error("Statistical function $stat_type not found.")
  end
end

function value_to_string(value)
  if value isa Number
    return string(value)
  elseif value isa Tuple
    return join(json.(value), ",")
  else
    return join(string.(value), " ")
  end
end

function file_to_dict(file_path)
  """
  Parse the option file into a dictionary.
  """
  params = Dict{String,Real}()
  try
    open(file_path) do file
      for line in readlines(file)
        key, val = split(line, ",")
        params[key] = parse_param(key, val)
      end
    end
  catch e
    println("Error reading file $file_path: $e")
    rethrow(e)
  end

  post_process_params(params)
  return params
end

function parse_param(key, val)
  if key in FLOAT_PARAMS
    return parse(Float64, val)
  elseif key in INT_PARAMS
    return parse(Int64, val)
  else
    error("Unknown parameter type for key $key")
  end
end

function post_process_params(params)
  dual_keys = [("dual_lambda", "lambda2"), ("dual_L", "Lt")]
  for (dual_key, target_key) in dual_keys
    if params[dual_key] == 1
      params[target_key] = params[target_key[1:end-1]]
    end
  end

  params["Lt"] = params["dual_L"] == 1 ? params["L"] : calculate_lt(params)
  delete_keys!(params, ["Lt_ratio", "Lt_power", "dual_lambda", "dual_L"])
end

function calculate_lt(params)
  Int64(params["Lt_ratio"] * (params["L"]^params["Lt_power"]))
end

function delete_keys!(dict, keys)
  for key in keys
    pop!(dict, key, nothing)
  end
end

# Other functions (print_header, print_raw_output) remain unchanged.

