include("src/MCMC.jl")
using .MCMC

function main(params_dict::Dict{String,Real})
  params = MCMC.SimParameters.SimParams(params_dict)
  sim = MCMC.Sim(params)
  println("Thermalizing...")
  MCMC.thermalize!(sim)
  println("Running...")
  MCMC.run!(sim)

  #Print out output
  MCMC.print_header()
  MCMC.print_raw_output(sim)
  return MCMC.get_measurements(sim.measurements)
end

input_params = MCMC.file_to_dict(ARGS[1])
data = main(input_params)

