# Markov Chain Monte Carlo simulation: Compressible Bosonic Integer Quantum Hall Effect 

This Julia package implemented the [worm](https://arxiv.org/abs/cond-mat/0103146) and Metropolis algorithm on the U(1)xU(1) model on a 3D cubic lattice.

This simulation is used to study the phase diagram of the compressible bosonic integer quantum hall state.

This simulation extends the beautiful work of [Scott D. Geraedts and Olexei I. Motrunich](https://arxiv.org/abs/1302.1436) by including a random chemical potential to the first Boson component and is a part of my [upcoming paper](https://kakkarav.com/publications/).

The U(1)xU(1) model is described by an effective action:
```math
$$ S = \sum_{\textbf{r}} \frac{1}{2 \lambda_1} \bigg|\textbf{J}(\textbf{r}) - \eta [ \nabla \times \textbf{p}] (\textbf{r}) \bigg|^2 
+ \sum_{\textbf{r}} \frac{\mu_0 + \mu(\textbf{r})}{\lambda_1} J_{\tau}(\textbf{r})
+ \sum_{\textbf{r}} \frac{1}{2\lambda_3} \bigg| \textbf{J}(\textbf{r}) \bigg|^2
+ \sum_\text{\textbf{R}} \frac{\lambda_2}{2} \bigg( [\nabla \theta] (\textbf{R})- 2 \pi \textbf{p} (\textbf{R})\bigg)^2$$
```
where $\textbf{J}$ is a loop bond variable, $\theta$ is an angle variable, and $\textbf{p}$ is an integer-valued gauge field.

## Installation

The dependencies can be installed by going into the root directory of the project and run the following command:

```julia
julia --project=.
```

In Julia REPL, run the following command to instantiate the environment:

```julia
> ]
> instantiate
```

## Parameters

-  `lambda1` (float): The coupling of the loop bond variable.
-  `lambda2` (float): The coupling of the angle variable.
-  `lambda3` (float): The additional self energy coupling 
- `eta` (float): The interspecies coupling.
- `delta_theta` (float): The angle update step size.
- `chem_potential` (float): The global chemical potential of the lop variable.
- `bond_type` (int): The type of bond disorder: binary or uniform distribution.
- `prob` (float): The probability of the missing bonds for the angle variable.
- `chem_disorder` (float): The strength of the chemical potential disorder (uniform distribution).
- `bond_disorder_villain` (float): The strength of the bond disorder.
- `bond_disorder_loop` (float): The strength of the bond disorder.
- `seed` (int): The seed for the random number generator.
- `L` (int): The spatial lattice dimension.
- `Lt` (int): The temporal lattice dimension.
- `num_of_thermal` (int): The number of thermalization sweeps.
- `num_of_measure` (int): The number of total measurements.
- `num_of_sweeps` (int): The number of Monte Carlo sweeps between measurement attempts.
- `nbin` (int): The number of bins for the measurements.

## Example

This example creates a lattice of size 12x12x12 with 10,000 measurements and 100,000 Monte Carlo sweeps.

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
# Warm up the simulation
MCMC.thermalize!(sim)
# Run the simulation
MCMC.run!(sim)
# the raw data is returned as a dictionary
println(MCMC.get_measurements(sim.measurements))
```

Alternatively, the script can also take an option file with the help of the parser function:
```julia
julia --project=. main.jl test.opts
```

