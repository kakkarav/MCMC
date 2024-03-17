using JSON

# Float parameters
const FLOATS = Set([
	"lambda1", "lambda2", "lambda3", "eta", "delta_theta", "chem_potential",
	"chem_disorder", "bond_disorder_villain", "bond_disorder_loop", "prob",
	"Lt_ratio", "Lt_power",
])

#Integer parameters
const INTS = Set([
	"seed", "bond_type", "L", "num_of_thermal", "num_of_measure",
	"num_of_sweeps", "num_of_relative_sweeps", "dual_lambda", "nbin", "dual_L",
])

# Parse an options file into a dictionary
function file_to_dict(file)
	params = Dict{String, Real}()
	open(file) do f
		for line in readlines(f)
			k, v = split(line, ",")
			try
				if k in FLOATS
					params[k] = parse(Float64, v)
				elseif k in INTS
					params[k] = parse(Int64, v)
				end
			catch e
				@error "Error parsing $k: $v"
				rethrow(e)
			end
		end
	end

	# Optional paramters controlling the lattice isotropy 
	params["lambda2"] = params["dual_lambda"] == 1 ? params["lambda1"] : params["lambda2"]
	params["Lt"] = params["dual_L"] == 1 ? params["L"] : Int64(params["Lt_ratio"] * (params["L"]^params["Lt_power"]))

	# Remove unused parameters
	for key in ["Lt_ratio", "Lt_power", "dual_lambda", "dual_L"]
		pop!(params, key, nothing)
	end

	return params
end

# Helper function to parse observations
function parseObs(meas::Dict{String, Tuple{Any, Any}}, name::String, type = "mean")
	value = meas[name]
	if type != "mean"
		value = getfield(Main, Symbol(type)).(value)
	end
	if eltype(value) <: Number
		return string(value)
	else
		return join(string.(value), " ")
	end
end

# Convert measurements to a comma-separated string
function output_to_string(meas::Dict{String, Tuple{Any, Any}}, obs::Vector{Tuple{String, String}})
	data = [parseObs(meas, pair...) for pair in obs]
	return join(data, ",")
end

# List of observables and the operator we want to act on the corresponding observables 
const OBSERVABLES = [
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

# Print out the headers for the raw data in CSV format
function print_header()
	"""
	Print out the headers for the raw data in the CSV format.
	"""
	header = ["sample_num"]
	for ob in OBSERVABLES
		push!(header, string(ob[1], "_", ob[2]))
		push!(header, string(ob[1], "_", ob[2], "_sq"))
	end
	println(join(header, ","))
end

# Print out the raw data to STDOUT in CSV format
function print_raw_output(sim::Sim)
	"""
	Print out the raw data to STDOUT in the CSV format.
	The data aligns with the header from the print_header function.
	"""
	data = MCMC.get_measurements(sim.measurements)
	sample_num = sim.measurements.measurements[4].obs_data.c_num_of_measure
	print(json(json(sample_num)), ",")
	println(output_to_string(data, OBSERVABLES))
end
