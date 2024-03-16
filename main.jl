include("src/MCMC.jl")

using .MCMC
using .MCMC: SimParameters

function main(params_dict::Dict{String,Real})
  params = SimParameters.SimParams(params_dict)
  sim = MCMC.Sim(params)
  MCMC.thermalize!(sim)
  MCMC.PrintHeader(params_dict)
  MCMC.run!(sim)
  return MCMC.get_measurements(sim.measurements)
end

input_params = MCMC.FileToDict(ARGS[1])
main(input_params)
