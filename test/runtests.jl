using Test
using TSMC
using TSMC: MCMC, SimParameters, Defs

sim_dict = Dict("lambda1" => 0.2,
    "lambda2" => 0.2,
    "lambda3" => 10000000000000.0,
    "eta" => 0.66,
    "delta_theta" => pi * 0.7,
    "chem_potential" => 0.0,
    "bond_type" => 1,
    "prob" => 0.0,
    "chem_disorder" => 0.65,
    "bond_disorder_villain" => 0.0,
    "bond_disorder_loop" => 0.0,
    "seed" => 2,
    "L" => 8,
    "Lt" => 8,
    "num_of_thermal" => 1000,
    "num_of_measure" => 100,
    "num_of_sweeps" => 10,
    "nbin" => 1,
)

sim_params = SimParameters.SimParams_new(sim_dict)
sim = MCMC.Sim(sim_params)
MCMC.thermalize!(sim)
MCMC.run!(sim)

include("Divergence.jl")
include("Energy.jl")
include("Curl.jl")
include("Angles.jl")
include("Stoke.jl")