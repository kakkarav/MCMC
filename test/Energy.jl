function energy_single_bond(
	lat::MCMC.Lattice,
	params::SimParameters.SimParams,
	disorder::MCMC.Disorder,
	location::NTuple{Defs.NDIMS, Int64},
	dir::Int64,
)
	nn = MCMC.hop(location, dir, 1, size(lat.angle))
	angle = lat.angle
	fluctuate =
		location[dir] == 3 ? lat.fluctuate_lattice[dir] / params.Lt :
		lat.fluctuate_lattice[dir] / params.L
	p_lattice = lat.p_lattice
	lambda1 =
		dir == Defs.NDIMS ? params.lambda1 :
		disorder.bond_disorder_villain[location[1], location[2], dir]
	lambda2 =
		dir == Defs.NDIMS ? params.lambda2 :
		disorder.bond_disorder_loop[location[1], location[2], dir]
	E_villain =
		(lambda1 / 2.0) *
		(
			angle[nn...] .- angle[location...] .- 2 * pi * (p_lattice[location..., dir]) .-
			fluctuate
		)^2
	E_loop =
		1 / (2.0 * lambda2) *
		(lat.current[location..., dir] - params.eta * MCMC.curl(lat, location, dir))^2
	E_disorder =
		dir == Defs.NDIMS ?
		disorder.chemical_potential[location[1], location[2]] *
		lat.current[location..., 3] / lambda2 : 0.0
	E_total = E_villain + E_loop + E_disorder
	return E_total
end


@testset "Energy" begin
	E = 0.0
	for j in eachindex(sim.lat.angle)
		j_site = Tuple(CartesianIndices(sim.lat.angle)[j])
		for dir in 1:Defs.NDIMS
			E += energy_single_bond(sim.lat, params, sim.disorder, j_site, dir)
		end
	end
	@test isapprox(E, sim.lat.energy, atol = 1e-10)
end
