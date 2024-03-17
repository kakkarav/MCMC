using MCMC
using Test
using MCMC: SimParameters, Defs

params_dict = Dict("lambda1" => 0.2,
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

params = SimParameters.SimParams(params_dict)
sim = MCMC.Sim(params)
MCMC.thermalize!(sim)
MCMC.run!(sim)

include("divergence.jl")
include("energy.jl")
include("curl.jl")
include("angles.jl")
include("stoke.jl")
