# Disorder stores information about spatially varying disorder in a lattice
struct Disorder
	chemical_potential::Array{Float64, NDIMS - 1}
	chemical_weight::Array{Float64, NDIMS - 1}
	bond_disorder_villain::Array{Float64, NDIMS}
	bond_disorder_loop::Array{Float64, NDIMS}
	weights_inc::Array{Array{Float64, 1}, NDIMS}
	weights_dec::Array{Array{Float64, 1}, NDIMS}
	sum_villain_bond::Array{Float64, 1}
	sum_loop_bond::Array{Float64, 1}
	sum_loop_bond_inverse::Array{Float64, 1}

	function Disorder(sim_params::SimParams)
		rng = Random.Xoshiro(rand(UInt64))
		dims = sim_params.L * ones(Int64, NDIMS - 1)
		chemical_potential_disorder = rand(dims...)
		chemical_potential = (chemical_potential_disorder * 2.0 .- 1) * sim_params.chem_disorder .+ sim_params.chem_potential
		bond_disorder_villain, bond_disorder_loop = generate_bond_disorder(sim_params, rng, dims)
		chemical_weight = exp.(-chemical_potential ./ sim_params.lambda2)
		weights_inc, weights_dec = calc_ratios(sim_params, bond_disorder_loop)
		sum_villain_bond = sum_bonds(bond_disorder_villain)
		sum_loop_bond = sum_bonds(1 ./ bond_disorder_loop)
		sum_loop_bond_inverse = sum_bonds(bond_disorder_loop)

		new(chemical_potential,
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
end

function generate_bond_disorder(sim_params::SimParams, rng, dims)
	bond_disorder_villain = Array{Float64, NDIMS}(undef, dims..., 2)
	bond_disorder_loop = Array{Float64, NDIMS}(undef, dims..., 2)

	if sim_params.bond_type == 1
		bond_disorder = 2 * rand(rng, dims..., 2) .- 1
		bond_disorder_villain .= (1.0 .+ bond_disorder * sim_params.bond_disorder_loop) .* sim_params.lambda1
		bond_disorder_loop .= (1.0 .+ bond_disorder * sim_params.bond_disorder_loop) .* sim_params.lambda2
	elseif sim_params.bond_type == 2
		random_values = rand(rng, dims..., 2)
		bond_disorder_villain .= map(x -> (x < sim_params.prob) ? sim_params.lambda1 / 3.0 : sim_params.lambda1, random_values)
		bond_disorder_loop .= map(x -> (x < sim_params.prob) ? sim_params.lambda2 / 3.0 : sim_params.lambda2, random_values)
	end

	return bond_disorder_villain, bond_disorder_loop
end

function sum_bonds(bond_array::Array{Float64, NDIMS})
	[sum(bond_array[:, :, i]) for i in 1:size(bond_array, 3)]
end

function calc_ratios(sim_params::SimParams, bond_loop::Array{Float64, NDIMS})
	w_inc = Array{Array{Float64, 1}}(undef, sim_params.L, sim_params.L, NDIMS - 1)
	w_dec = Array{Array{Float64, 1}}(undef, sim_params.L, sim_params.L, NDIMS - 1)

	for site in eachindex(bond_loop)
		location = Tuple(CartesianIndices(bond_loop)[site])
		coupling = 1 / (2.0 * bond_loop[location...])
		energy_inc = coupling .* [((n + 1)^2 - n^2) for n in -MAX_BOND_J:MAX_BOND_J]
		energy_dec = coupling .* [((n - 1)^2 - n^2) for n in -MAX_BOND_J:MAX_BOND_J]
		w_inc[location...] = exp.(-energy_inc)
		w_dec[location...] = exp.(-energy_dec)
	end

	return w_inc, w_dec
end
