function metro_angle!(
	lattice::Lattice,
	params::SimParams,
	disorder::Disorder,
	rng::Random.Xoshiro,
	location::NTuple{NDIMS, Int64},
)
	angle_change = 2.0 * (rand(rng) - 0.5) * params.delta_theta
	dE =
		Energy_local_angle(lattice, params, disorder, location, angle_change) -
		Energy_local_angle(lattice, params, disorder, location, 0.0)
	if (dE <= 0.0 || rand(rng) <= exp.(-dE))
		lattice.angle[location...] = update_angle(lattice.angle[location...], angle_change)
		lattice.energy += dE
	end
	return nothing
end

function update_angle(angle::Float64, delta::Float64)
	new_angle = mod(angle + delta, 2 * pi)
	if new_angle > pi
		new_angle = new_angle - 2 * pi
	end
	return new_angle
end

function Energy_local_angle(
	lattice::Lattice,
	params::SimParams,
	disorder::Disorder,
	location::NTuple{NDIMS, Int64},
	delta_angle::Float64,
)
	new_angle = update_angle(lattice.angle[location...], delta_angle)
	#by changing an angle we change NDIMS*2 bonds states and energy, so the energy change calculation is looped over them
	E = 0.0
	fluataute = 0.0
	#For bonds in the XY plane
	#Disordered bonds are store in Disorder.
	for dir in 1:2
		for lc ∈ 1:2
			#we hop here in order to calculate \omega_{\mu} = phi(r+\mu) - phi(r)
			nn = hop(location, dir, lc, size(lattice.angle))
			if lc == 1
				fluctuate = lattice.fluctuate_lattice[dir] / params.L
				E +=
					disorder.bond_disorder_villain[location[1], location[2], dir] *
					(
						lattice.angle[nn...] - new_angle -
						2 * pi * lattice.p_lattice[location..., dir] - fluctuate
					)^2.0
			elseif lc == 2
				fluctuate = lattice.fluctuate_lattice[dir] / params.L
				E +=
					disorder.bond_disorder_villain[nn[1], nn[2], dir] *
					(
						new_angle - lattice.angle[nn...] -
						2 * pi * lattice.p_lattice[nn..., dir] - fluctuate
					)^2.0
			end
		end
	end

	#for bond in the time direction (dir==3)
	for lc ∈ 1:2
		#we hop here in order to calculate \omega_{\mu} = phi(r+\mu) - phi(r)
		nn = hop(location, NDIMS, lc, size(lattice.angle))
		if lc == 1
			fluctuate = lattice.fluctuate_lattice[3] / params.Lt
			E +=
				params.lambda1 *
				(
					lattice.angle[nn...] - new_angle -
					2 * pi * lattice.p_lattice[location..., NDIMS] - fluctuate
				)^2.0
		elseif lc == 2
			fluctuate = lattice.fluctuate_lattice[3] / params.Lt
			E +=
				params.lambda1 *
				(
					new_angle - lattice.angle[nn...] -
					2 * pi * lattice.p_lattice[nn..., NDIMS] - fluctuate
				)^2.0
		end
	end

	return E / 2.0
end
