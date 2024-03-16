module SimParameters

export SimParams_new, SimParams

mutable struct SimParams
  # coupling for Villain model
  lambda1::Float64
  # coupling for loop-current binding model
  lambda2::Float64
  # coupling for loop-current model
  lambda3::Float64
  #coupling constant : 0-turn off binding, 1-turn on binding
  eta::Float64
  #chemical potential
  chem_potential::Float64
  #chemical potential disorder strength
  chem_disorder::Float64
  #type of disorder: 1 Box , 2 percolation
  bond_type::Int64
  #type of disorder: 1 Box , 2 percolation
  prob::Float64
  #bond disorder strength
  bond_disorder_villain::Float64
  #bond disorder strength
  bond_disorder_loop::Float64
  # Length
  L::Int64
  # Length in the temporal direction
  Lt::Int64
  #angle range for each step
  delta_theta::Float64
  # num of thermalization steps
  num_of_thermal::Int64
  # num of sweeps between measurements
  num_of_sweeps::Int64
  # num of measurements
  num_of_measure::Int64
  #number of relative sweeps between Villain and loop model
  num_of_relative_sweeps::Int64
  # rng seed
  seed::Int64
  # number of bins
  nbin::Int64
  # unintizlied constructor
  SimParams() = new()
end

function SimParams(params::Dict{String,Real})
  newSimParams = SimParams()
  for f in keys(params)
    setfield!(newSimParams, Symbol(f), params[f])
  end
  return newSimParams
end

end
