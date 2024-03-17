function divergence(lattice::Array{Int64, 4}, params::SimParameters.SimParams, site::NTuple{Defs.NDIMS, Int64})
	total_current = 0
	for dir in 1:Defs.NDIMS
		for lr in 1:2
			current_site = lr == 1 ? site : MCMC.hop(site, dir, 2, size(lattice)[1:3])
			charge = (lr == 1 ? 1 : -1)
			total_current += charge * lattice[current_site..., dir]
		end
	end
	return total_current
end

@testset "Current Divergence" begin
	for site in eachindex(sim.lat.angle)
		site = Tuple(CartesianIndices(sim.lat.angle)[site])
		@test divergence(sim.lat.current, sim.sim_params, site) == 0
	end
end

@testset "Curl lattice Divergence" begin
	for site in eachindex(sim.lat.angle)
		site = Tuple(CartesianIndices(sim.lat.angle)[site])
		@test divergence(sim.lat.curl_lattice, sim.sim_params, site) == 0
	end
end
