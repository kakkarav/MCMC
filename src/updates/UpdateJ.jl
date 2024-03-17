struct WormMoveData
	weights_inc_with_p::Array{Array{Float64, 1}, 1}
	weights_dec_with_p::Array{Array{Float64, 1}, 1}
	weights_self_inc::Array{Float64, 1}
	weights_self_dec::Array{Float64, 1}
	rng::Random.Xoshiro
	function WormMoveData(sim_params::SimParams, rng::Random.Xoshiro)
		w_inc, w_dec, weights_self_inc, weights_self_dec = calc_ratios(sim_params)
		new(w_inc, w_dec, weights_self_inc, weights_self_dec, rng)
	end
end


function calc_ratios(sim_params::SimParams)
	"""
	Return the bond weight ratios for a given patameter up to order MAX_BOND
	"""

	weights_inc_with_p = Array{Array{Float64, 1}}(undef, 2 .* MAX_P + 1)
	weights_dec_with_p = Array{Array{Float64, 1}}(undef, 2 .* MAX_P + 1)
	for curl_p in -MAX_P:MAX_P
		location = curl_p + MAX_P + 1
		coupling1 = 1 / (2.0 * sim_params.lambda2)
		coupling2 = sim_params.eta * curl_p ./ sim_params.lambda2
		energy_inc =
			coupling1 * [((n + 1)^2.0 - n^2.0) for n in -MAX_BOND_J:MAX_BOND_J] .- coupling2
		weights_inc = exp.(-energy_inc)
		energy_dec =
			coupling1 * [((n - 1)^2.0 - n^2.0) for n in -MAX_BOND_J:MAX_BOND_J] .+ coupling2
		weights_dec = exp.(-energy_dec)
		weights_inc_with_p[location] = weights_inc
		weights_dec_with_p[location] = weights_dec
	end

	coupling2 = 1 / (2.0 * sim_params.lambda3)
	energy_inc = coupling2 * [((n + 1)^2.0 - n^2.0) for n in -MAX_BOND_J:MAX_BOND_J]
	weights_self_inc = exp.(-energy_inc)
	energy_dec = coupling2 * [((n - 1)^2.0 - n^2.0) for n in -MAX_BOND_J:MAX_BOND_J]
	weights_self_dec = exp.(-energy_dec)

	return weights_inc_with_p, weights_dec_with_p, weights_self_inc, weights_self_dec
end


function update_energy!(
	lat::Lattice,
	params::SimParams,
	disorder::Disorder,
	dir::Int64,
	lr::Int64,
	location::NTuple{NDIMS, Int64},
)
	"""
	update_energy:update the energy of the lattice after the worm moves
	"""
	nn = lr == 1 ? location : hop(location, dir, 2, size(lat.angle))
	bond_val = lat.current[nn..., dir] #current variable
	c_val = lat.curl_lattice[nn..., dir] #curl variable
	effective_bond = bond_val - params.eta * c_val
	dE = 0.0
	# Binding
	if lr == 1
		dE1 = (effective_bond + 1)^2.0 - (effective_bond)^2.0
	elseif lr == 2
		dE1 = (effective_bond - 1)^2.0 - (effective_bond)^2.0
	else
		println("Error")
	end

	if dir == NDIMS
		coupling = 1 / (2 * params.lambda2)
	elseif dir == 1 || dir == 2
		coupling = 1 / (2 * disorder.bond_disorder_loop[nn[1], nn[2], dir])
	else
		println("Error")
	end

	dE += dE1 * coupling

	# chemical potential
	if dir == NDIMS
		if lr == 1
			dE3 = disorder.chemical_potential[nn[1], nn[2]] / params.lambda2
		elseif lr == 2
			dE3 = -disorder.chemical_potential[nn[1], nn[2]] / params.lambda2
		else
			println("Error")
		end
		dE += dE3
	end

	return dE
end




function shift_worm!(
	worm_data::WormMoveData,
	lat::Lattice,
	params::SimParams,
	disorder::Disorder,
)
	"""
	single shift move of worm's head
	"""
	# draw next move
	move_rnd = rand(worm_data.rng, WORM_RANGE)
	# cartesian direction
	dir = fld1(move_rnd, 2)
	# left or right
	lr = mod1(move_rnd, 2)
	# nearest neighbor
	nn = hop(lat.head, dir, lr, size(lat.angle))

	weight = 0.0
	weight_final = 0.0
	#weight from the bond, no on-site weight
	if dir == NDIMS
		if lr == 1
			#current J
			bond_val = lat.current[lat.head..., dir]
			#curl of P
			c_val = lat.curl_lattice[lat.head..., dir]
			chem_weight = disorder.chemical_weight[lat.head[1], lat.head[2]]
			#switch worm_data to disorder
			loop_weight =
				worm_data.weights_inc_with_p[(c_val.+MAX_P.+1)][(bond_val.+MAX_BOND_J.+1)]
			# loop_weight = worm_data.weights_inc_with_p[1][1]
			weight = loop_weight * chem_weight

		elseif lr == 2
			bond_val = lat.current[nn..., dir]
			c_val = lat.curl_lattice[nn..., dir]
			chem_weight = disorder.chemical_weight[nn[1], nn[2]]
			loop_weight =
				worm_data.weights_dec_with_p[(c_val.+MAX_P.+1)][(bond_val.+MAX_BOND_J.+1)]
			# loop_weight = worm_data.weights_inc_with_p[1][1]
			weight = loop_weight ./ chem_weight
		end

		#chemical potential
		#time direction index = NDIMS
	elseif dir == 1 || dir == 2
		if lr == 1
			bond_val = lat.current[lat.head..., dir]
			c_val = lat.curl_lattice[lat.head..., dir]
			weight =
				worm_data.weights_inc_with_p[(c_val.+MAX_P.+1)][(bond_val.+MAX_BOND_J.+1)]
		elseif lr == 2
			bond_val = lat.current[nn..., dir]
			c_val = lat.curl_lattice[nn..., dir]
			weight =
				worm_data.weights_dec_with_p[(c_val.+MAX_P.+1)][(bond_val.+MAX_BOND_J.+1)]
		end
	end

	accept = (weight_final > rand(worm_data.rng))
	if (accept)
		# accept
		# E_before = lat.energy
		dE = update_energy!(lat, params, disorder, dir, lr, lat.head)
		lat.energy += dE


		if lr == 1
			#move to the right
			lat.current[lat.head..., dir] += 1
			lat.winding[dir] += 1
		elseif lr == 2
			# move to the left
			lat.current[nn..., dir] -= 1
			lat.winding[dir] -= 1
		end
		lat.head = nn
	end
	return nothing
end

function jump_worm!(worm_data::WormMoveData, lat::Lattice)
	"""
	jump of worm's head and tail
	"""
	# new_site=ind2sub(size(lat.angle),rand(worm_data.rng,1:length(lat.angle)))
	new_site = Tuple(CartesianIndices(lat.angle)[rand(worm_data.rng, 1:length(lat.angle))])
	if new_site != lat.head
		lat.head = new_site
		lat.tail = new_site
	end
	return nothing
end
