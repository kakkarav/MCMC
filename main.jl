include("src/MCMC.jl")
using .MCMC

function main(params_dict::Dict{String,Real})
  params = MCMC.SimParameters.SimParams(params_dict)
  println("Initializing...")
  sim = MCMC.Sim(params)
  println("Thermalizing...")
  MCMC.thermalize!(sim)
  println("Running...")
  MCMC.run!(sim)

  #Print raw measurements 
  MCMC.print_header()
  MCMC.print_raw_output(sim)
  return MCMC.get_measurements(sim.measurements)
end

params_dict = MCMC.file_to_dict(ARGS[1])
data = main(params_dict)

