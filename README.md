# Markov Chain Monte Carlo simulation for compressible Bosonic Integer quantum hall 

This Julia package implemented the [worm](https://arxiv.org/abs/cond-mat/0103146) and Metropolis algorithm on the U(1)xU(1) model on a 3D cubic lattice.

This simulation is used to study the phase diagram of the compressible bosonic integer quantum hall state.

This simulation extends the beautiful work of [Scott D. Geraedts and Olexei I. Motrunich](https://arxiv.org/abs/1302.1436) by including a random chemical potential and is a part of my new paper.


## Installation

The dependencies can be installed by going into the root directory of the project and run the following command:

```julia
julia --project=.
]
instantiate
```

## Example

This example create a lattice of size 8x8x8 with 1000 measurements and 10000 thermalization sweeps.

The parameters are stored in a dictionary and the simulation is run as follows:
```julia
using MCMC

params_dict = Dict("lambda1" => 0.23,
  "lambda2" => 0.23,
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
  "L" => 12,
  "Lt" => 12,
  "num_of_thermal" => 100000,
  "num_of_measure" => 10000,
  "num_of_sweeps" => 1,
  "nbin" => 1,
)

params = MCMC.SimParameters.SimParams(params_dict)
sim = MCMC.Sim(params)
MCMC.thermalize!(sim)
MCMC.run!(sim)
# the raw data is return as a dictionary
println(MCMC.get_measurements(sim.measurements))
```

Alternatively, the code can be run by taking an option file with the help of the parser function:
```julia
julia --project=. main.jl test.opts
```

