struct Disorder
  chemical_potential::Array{Float64,NDIMS - 1}
  chemical_weight::Array{Float64,NDIMS - 1}
  bond_disorder_villain::Array{Float64,NDIMS}
  bond_disorder_loop::Array{Float64,NDIMS}
  weights_inc::Array{Array{Float64,1},NDIMS}
  weights_dec::Array{Array{Float64,1},NDIMS}
  sum_villain_bond::Array{Float64,1}
  sum_loop_bond::Array{Float64,1}
  sum_loop_bond_inverse::Array{Float64,1}
end

function Disorder(sim_params::SimParams)
  rng = Random.Xoshiro(rand(UInt64))
  chemical_potential_disorder = rand(sim_params.L * ones(Int64, NDIMS - 1)...)
  chemical_potential =
    (chemical_potential_disorder * 2.0 .- 1) * sim_params.chem_disorder .+
    sim_params.chem_potential
  prob = sim_params.prob

  bond_struct = sim_params.bond_type

  if (bond_struct == 1)
    bond_disorder_loop =
      2 * rand(rng, sim_params.L * ones(Int64, NDIMS - 1)..., 2) .- 1
    bond_disorder_villain =
      (1.0 .+ bond_disorder_loop * sim_params.bond_disorder_loop) .*
      sim_params.lambda1
    bond_disorder_loop =
      (1.0 .+ bond_disorder_loop * sim_params.bond_disorder_loop) .*
      sim_params.lambda2
  elseif (bond_struct == 2)
    bond_disorder_villain = map(
      x -> (x < prob) ? sim_params.lambda1 ./ 3.0 : sim_params.lambda1,
      rand(rng, sim_params.L * ones(Int64, NDIMS - 1)..., 2),
    )
    bond_disorder_loop = map(
      x -> (x < prob) ? sim_params.lambda2 ./ 3.0 : sim_params.lambda2,
      rand(rng, sim_params.L * ones(Int64, NDIMS - 1)..., 2),
    )
  end

  chemical_weight = exp.(-chemical_potential ./ sim_params.lambda2)
  weights_inc, weights_dec = calc_ratios(sim_params, bond_disorder_loop)
  sum_villain_bond =
    [sum(bond_disorder_villain[:, :, 1]), sum(bond_disorder_villain[:, :, 2])]
  sum_loop_bond =
    [sum(1 ./ bond_disorder_loop[:, :, 1]), sum(1 ./ bond_disorder_loop[:, :, 2])]
  sum_loop_bond_inverse = [
    sum(1.0 ./ bond_disorder_loop[:, :, 1]),
    sum(1.0 ./ bond_disorder_loop[:, :, 2]),
  ]

  Disorder(
    chemical_potential,
    chemical_weight,
    bond_disorder_villain,
    bond_disorder_loop,
    weights_inc,
    weights_dec,
    sum_villain_bond,
    sum_loop_bond,
    sum_loop_bond_inverse,
  )
end

function calc_ratios(sim_params::SimParams, bond_loop::Array{Float64,NDIMS})
  """
  Calculates the ratios of the weights for the loop update
  """
  w_inc = Array{Array{Float64,1}}(undef, sim_params.L, sim_params.L, NDIMS - 1)
  w_dec = Array{Array{Float64,1}}(undef, sim_params.L, sim_params.L, NDIMS - 1)

  for site in eachindex(bond_loop)
    location = Tuple(CartesianIndices(bond_loop)[site])

    coupling = 1 / (2.0 * bond_loop[location...])

    energy_inc = coupling * [((n + 1)^2.0 .- n^2.0) for n = -MAX_BOND_J:MAX_BOND_J]
    weights_inc = exp.(-energy_inc)

    energy_dec = coupling * [((n - 1)^2.0 .- n^2.0) for n = -MAX_BOND_J:MAX_BOND_J]
    weights_dec = exp.(-energy_dec)

    w_inc[location...] = weights_inc
    w_dec[location...] = weights_dec
  end

  return w_inc, w_dec
end
