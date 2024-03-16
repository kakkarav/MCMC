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
  MCMC.printHeader(params_dict)
  MCMC.printOutput(sim)
  return MCMC.get_measurements(sim.measurements)
end

input_params = MCMC.FileToDict(ARGS[1])
data = main(input_params)
